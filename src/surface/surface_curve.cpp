#include "geometrycentral/surface/surface_curve.h"

// Public Utils

SurfaceCurve::SurfaceCurve(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry,
                           GeodesicAlgorithmExact* mmpSolver)
    : m_linearScaleCoefficient(1.0) {
  m_mesh.reset(mesh);
  m_geometry.reset(geometry);
  m_mmpSolver.reset(mmpSolver);
}

void SurfaceCurve::createFromPoints(std::vector<SurfacePoint>& curvePoints) {
  assert(curvePoints.size() > 1);

  m_points.clear();

  m_startPoint = curvePoints[0];

  std::vector<SurfacePoint> denseCurve;

  size_t N = curvePoints.size();

  SurfacePoint pt1, pt2;
  double dist;

  for (size_t i = 0; i < N - 1; i++) {
    pt1 = curvePoints[i];
    pt2 = curvePoints[i + 1];

    std::vector<SurfacePoint> path = connectPointsWithGeodesic(pt1, pt2, dist);
    std::vector<SurfacePoint> prunedPath = pruneApproxEqualEntries(path);

    denseCurve.insert(denseCurve.end(), prunedPath.begin(), prunedPath.end());
    denseCurve.pop_back(); // to avoid double counti
  }

  m_points = denseCurve;

  computeInitialDirection();
  computeAnglesAndDistances();
}

void SurfaceCurve::deformAngle(int pointIndex, std::complex<double> newDir) {
  assert(pointIndex < m_angles.size());
  m_angles[pointIndex] = newDir;
  recompute();
}

void SurfaceCurve::deformDistance(int pointIndex, double newDist) {
  assert(pointIndex < m_distances.size());
  m_distances[pointIndex] = newDist;
  recompute();
}

void SurfaceCurve::getPoints(std::vector<SurfacePoint>& points) { points = m_points; }

std::complex<double> SurfaceCurve::getAngleAtIndex(int index) {
  assert(index < m_angles.size());
  return m_angles[index];
}

double SurfaceCurve::getDistanceAtIndex(int index) {
  assert(index < m_distances.size());
  return m_distances[index];
}

double SurfaceCurve::getLinearScaleCoefficient() { return m_linearScaleCoefficient; }

void SurfaceCurve::getParameterized(std::vector<CurvePointParams>& parameterizedPoints) {
  size_t N = m_angles.size();

  std::vector<CurvePointParams> pCurve(N);

  // First axis point direction and all last point parameters will be ignored
  for (int i = 0; i < N; i++) {
    CurvePointParams intermediateParams;

    intermediateParams.nextPointDistance = m_distances[i];
    intermediateParams.nextPointDirection = m_angles[i];

    pCurve[i] = intermediateParams;
  }

  parameterizedPoints = pCurve;
}

SurfacePoint SurfaceCurve::getPointAtIndex(int index) {
  assert(index < m_points.size());
  return m_points[index];
}

int SurfaceCurve::getSize() { return m_points.size(); }

Vector2 SurfaceCurve::getStartDir() { return m_initDir; }

std::complex<double> SurfaceCurve::getTangentAtIndex(int index) {
  int nextDir = (index == m_points.size() - 1) ? -1 : 1;
  SurfacePoint pt1 = m_points[index];
  SurfacePoint pt2 = m_points[index + nextDir];
  std::complex<double> ref = localDir(pt1, pt2);
  return ref;
}

void SurfaceCurve::rotate(Vector2 newDir) {
  m_initDir = newDir;
  recompute();
}

void SurfaceCurve::saveStartParams() {
  m_savedStartPoint = m_startPoint;
  m_savedInitDir = m_initDir;
}

void SurfaceCurve::setLinearScaleCoefficient(double linearScaleCoefficient) {
  m_linearScaleCoefficient = linearScaleCoefficient;
  recompute();
}

void SurfaceCurve::setParameterized(std::vector<CurvePointParams>& parameterizedPoints) {
  assert(parameterizedPoints.size() == m_points.size());

  m_distances.clear();
  m_angles.clear();

  for (int i = 0; i < parameterizedPoints.size(); i++) {
    CurvePointParams pp = parameterizedPoints[i];

    m_distances.push_back(pp.nextPointDistance);
    m_angles.push_back(pp.nextPointDirection);
  }
}

void SurfaceCurve::setPoints(std::vector<SurfacePoint>& points) {
  assert(points.size() > 0);

  m_points.clear();

  m_points.insert(m_points.end(), points.begin(), points.end());

  m_startPoint = m_points[0];

  computeInitialDirection();
  computeAnglesAndDistances();
}

void SurfaceCurve::transfer(SurfaceCurve* target, const SurfacePoint& targetMeshStart,
                            const SurfacePoint& targetMeshDirEndpoint) {
  target->m_angles.clear();
  target->m_distances.clear();

  target->m_angles.insert(target->m_angles.end(), m_angles.begin(), m_angles.end());
  target->m_distances.insert(target->m_distances.end(), m_distances.begin(), m_distances.end());

  target->m_startPoint = targetMeshStart;
  target->m_initDir = target->localDir(targetMeshStart, targetMeshDirEndpoint);

  target->recompute();
}

void SurfaceCurve::translate(const SurfacePoint& newStartPoint) {
  assert(m_savedStartPoint != SurfacePoint());
  assert(m_savedInitDir != Vector2());

  m_initDir = parallelTransport(m_savedStartPoint, newStartPoint, m_savedInitDir.arg());
  m_startPoint = newStartPoint;

  recompute();
}

// Private Utils

bool SurfaceCurve::approxEqual(const SurfacePoint& pA, const SurfacePoint& pB) {

  if (pA.type != pB.type) return false;
  double eps = 1e-5;
  switch (pA.type) {
  case (SurfacePointType::Vertex):
    return pA.vertex == pB.vertex;
    break;
  case (SurfacePointType::Edge):
    return pA.edge == pB.edge && abs(pA.tEdge - pB.tEdge) < eps;
    break;
  case (SurfacePointType::Face):
    return pA.face == pB.face && (pA.faceCoords - pB.faceCoords).norm() < eps;
    break;
  }
  throw std::logic_error("bad switch"); // shouldn't get here
}

void SurfaceCurve::computeAnglesAndDistances() {
  size_t N = m_points.size();

  std::vector<std::complex<double>> angles(N);
  std::vector<double> distances(N);

  SurfacePoint currNode, prevNode, nextNode;
  std::complex<double> origDir, offsetDir, adjustedDir;
  size_t currIdx, nextIdx, prevIdx;

  for (size_t i = 1; i < N - 1; i++) {
    currIdx = i;
    nextIdx = currIdx + 1;
    prevIdx = currIdx - 1;
    currNode = m_points[currIdx];
    nextNode = m_points[nextIdx];
    prevNode = m_points[prevIdx];

    Vector3 a = currNode.interpolate(m_geometry->inputVertexPositions);
    Vector3 b = nextNode.interpolate(m_geometry->inputVertexPositions);

    distances[i] = (b - a).norm();

    origDir = localDir(currNode, nextNode);
    offsetDir = localDir(currNode, prevNode);

    // CCW angle of rotation to get to next tangent vector; for convex corners, rotation is CCW (and angle is
    // positive), CW for non-convex corners (angle is negative).
    adjustedDir = origDir / offsetDir;
    angles[i] = adjustedDir;
  }

  Vector3 startPos = m_points[0].interpolate(m_geometry->inputVertexPositions);
  Vector3 nextPos = m_points[1].interpolate(m_geometry->inputVertexPositions);

  distances[0] = (nextPos - startPos).norm();

  m_angles = angles;
  m_distances = distances;
}

void SurfaceCurve::computeInitialDirection() {
  m_geometry->requireVertexNormals();
  m_geometry->requireVertexTangentBasis();

  SurfacePoint firstPoint = m_points[0];
  SurfacePoint secondPoint = m_points[1];

  m_initDir = localDir(firstPoint, secondPoint);
  m_initDir /= m_initDir.norm();
}

std::vector<SurfacePoint> SurfaceCurve::connectPointsWithGeodesic(const SurfacePoint& pt1, const SurfacePoint& pt2,
                                                                  double& distance) {

  // If the source and target SurfacePoints already share a common face, no need to run the algorithm.
  if (sharedFace(pt1, pt2) != Face()) {
    VertexData<Vector3>& vertexPositions = m_geometry->vertexPositions;
    distance = (pt1.interpolate(vertexPositions) - pt2.interpolate(vertexPositions)).norm();
    return {pt1, pt2};
  }

  // Otherwise, use MMP to compute exact geodesic paths.
  SurfacePoint source = pt1;
  double max_propagation_distance = GEODESIC_INF;
  std::vector<SurfacePoint> stop_points = {pt2};
  m_mmpSolver->propagate(source, max_propagation_distance, stop_points);
  distance = m_mmpSolver->getDistance(pt2);
  std::vector<SurfacePoint> path = m_mmpSolver->traceBack(pt2); // gets path from query point to source
  std::reverse(path.begin(), path.end());

  return path;
}

Vector2 SurfaceCurve::inTangentBasis(const SurfacePoint& pA, const SurfacePoint& pB, const SurfacePoint& p) {

  Vector2 vecP;
  Vector3 posA = pA.interpolate(m_geometry->vertexPositions);
  Vector3 posB = pB.interpolate(m_geometry->vertexPositions);
  Vector3 vec = posB - posA;

  switch (p.type) {
  case (SurfacePointType::Vertex): {
    Vertex vP = p.vertex;
    Vector3 basisX = m_geometry->vertexTangentBasis[vP][0];
    Vector3 basisY = m_geometry->vertexTangentBasis[vP][1];
    vecP = {dot(vec, basisX), dot(vec, basisY)};
    break;
  }
  case (SurfacePointType::Edge): {
    Edge eP = p.edge;
    Face f = sharedFace(pA, pB); // face containing vec
    Vector3 fN = m_geometry->faceNormals[f];
    Vector3 basisX = m_geometry->halfedgeVector(eP.halfedge()).normalize();
    Vector3 basisY = cross(basisX, fN);
    vecP = {dot(vec, basisX), dot(vec, basisY)};
    break;
  }
  case (SurfacePointType::Face): {
    Face fP = p.face;
    Vector3 basisX = m_geometry->faceTangentBasis[fP][0];
    Vector3 basisY = m_geometry->faceTangentBasis[fP][1];
    vecP = {dot(vec, basisX), dot(vec, basisY)};
    break;
  }
  }
  return vecP.normalize();
}

Vector2 SurfaceCurve::localDir(const SurfacePoint& pt1, const SurfacePoint& pt2) {

  // Get the local basis vectors in global coordinates
  Halfedge heDir;
  Vector3 normal;
  Vector3 globalDir =
      pt2.interpolate(m_geometry->inputVertexPositions) - pt1.interpolate(m_geometry->inputVertexPositions);
  globalDir /= globalDir.norm();
  m_geometry->requireVertexNormals(); // GC doesn't seem to have immediate method for vertex normals
  double eps = 1e-8;

  if (pt1.type == SurfacePointType::Face) {
    heDir = pt1.face.halfedge();
    normal = m_geometry->faceNormal(pt1.face);
  } else if (pt1.type == SurfacePointType::Edge) {
    Edge eDir = pt1.edge;
    heDir = eDir.halfedge();
    Vector3 local_x = m_geometry->halfedgeVector(heDir).normalize();
    // If segment lies on the edge
    if (abs(1.0 - dot(local_x, globalDir)) < eps) {
      return {1.0, 0.0};
    }
    // If segment lies on the edge
    else if (abs(-1.0 - dot(local_x, globalDir)) < eps) {
      return {-1.0, 0.0};
    }
    // If segment lies on a face, find which face it lies in, and get its normal
    bool found = false;
    if (pt2.type == SurfacePointType::Face) {
      Face f = pt2.face;
      normal = m_geometry->faceNormal(f);
    } else if (pt2.type == SurfacePointType::Edge) {
      for (Face f : eDir.adjacentFaces()) {
        for (Edge e : f.adjacentEdges()) {
          if (pt2.edge == e) {
            normal = m_geometry->faceNormal(f);
            found = true;
            break;
          }
        }
        if (found) break;
      }
    } else {
      for (Face f : eDir.adjacentFaces()) {
        for (Vertex v : f.adjacentVertices()) {
          if (pt2.vertex == v) {
            normal = m_geometry->faceNormal(f);
            found = true;
            break;
          }
        }
        if (found) break;
      }
    }
  } else {
    heDir = pt1.vertex.halfedge();
    normal = m_geometry->vertexNormals[pt1.vertex]; // angle-weighted normal
  }

  Vector3 local_x = m_geometry->halfedgeVector(heDir).normalize();
  Vector3 local_y = cross(normal, local_x);

  Vector2 dir = {dot(globalDir, local_x), dot(globalDir, local_y)}; // not necessarily unit
  dir /= dir.norm();
  return dir;
}

Vector2 SurfaceCurve::parallelTransport(const SurfacePoint& startPoint, const SurfacePoint& endpoint,
                                        double initAngle) {
  double distance;

  std::vector<SurfacePoint> path = connectPointsWithGeodesic(startPoint, endpoint, distance);

  m_geometry->requireTransportVectorsAlongHalfedge();
  m_geometry->requireVertexTangentBasis();
  m_geometry->requireFaceTangentBasis();
  m_geometry->requireFaceNormals();

  Vector2 result = Vector2::fromAngle(initAngle);

  size_t nNodes = path.size();
  for (size_t i = 0; i < nNodes - 1; i++) {
    SurfacePoint pA = path[i];
    SurfacePoint pB = path[i + 1];
    assert(sharedFace(pA, pB) != Face());

    if (pA.type == SurfacePointType::Face && pB.type == SurfacePointType::Face && pA.face == pB.face) continue;
    if (pA.type == SurfacePointType::Edge && pB.type == SurfacePointType::Edge && pA.edge == pB.edge) continue;

    Vector2 rot;
    if (pA.type == SurfacePointType::Vertex && pB.type == SurfacePointType::Vertex) {
      // Determine the halfedge between pA.vertex -> pB.vertex
      Halfedge heAB = Halfedge();
      for (Halfedge he : pA.vertex.outgoingHalfedges()) {
        if (he.tipVertex() == pB.vertex) heAB = he;
      }
      assert(heAB != Halfedge());
      rot = m_geometry->transportVectorsAlongHalfedge[heAB];
    } else {
      Vector2 vecA = inTangentBasis(pA, pB, pA);
      Vector2 vecB = inTangentBasis(pA, pB, pB);
      rot = vecB / vecA;
    }
    result = rot * result;
  }

  return result;
}

std::vector<SurfacePoint> SurfaceCurve::pruneApproxEqualEntries(const std::vector<SurfacePoint>& source) {
  std::vector<SurfacePoint> result;

  if (source.size() == 2) {
    result = source;
  } else {
    result.push_back(source[0]);

    for (const SurfacePoint& pathPt : source) {
      if (!approxEqual(pathPt, result.back())) {
        result.push_back(pathPt);
      }
    }
  }

  return result;
}

void SurfaceCurve::recompute() {
  // Clear all existing point data
  m_points.clear();

  // Place first points
  m_points.push_back(m_startPoint);

  // Declare some variables
  size_t N = m_angles.size();

  SurfacePoint pathEndpoint;
  std::complex<double> endingDir;
  TraceGeodesicResult tracedGeodesic;
  TraceOptions traceOptions;
  traceOptions.includePath = true;

  // Trace out from first to second point
  m_initDir /= m_initDir.norm();
  tracedGeodesic =
      traceGeodesic(*(m_geometry), m_startPoint,
                    Vector2::fromComplex(m_initDir * m_linearScaleCoefficient * m_distances[0]), traceOptions);

  pathEndpoint = tracedGeodesic.endPoint;
  m_points.push_back(pathEndpoint);
  endingDir = -tracedGeodesic.endingDir;

  // Trace the rest
  for (size_t i = 1; i < N - 1; i++) {
    // Get the corresponding direction on S2
    endingDir /= std::abs(endingDir);
    std::complex<double> adjustedDirS2 = endingDir * m_angles[i];
    adjustedDirS2 /= std::abs(adjustedDirS2); // make sure direction is unit

    // Trace geodesic from last point and record where it ended up
    tracedGeodesic =
        traceGeodesic(*(m_geometry), m_points[m_points.size() - 1],
                      Vector2::fromComplex(adjustedDirS2 * m_linearScaleCoefficient * m_distances[i]), traceOptions);

    pathEndpoint = tracedGeodesic.endPoint;
    m_points.push_back(pathEndpoint);
    endingDir = -tracedGeodesic.endingDir;
  }
}
