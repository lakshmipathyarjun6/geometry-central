#include "geometrycentral/surface/surface_patch.h"

int mod(int a, int b) { return (b + (a % b)) % b; }

SurfacePatch::SurfacePatch(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry,
                           GeodesicAlgorithmExact* mmpSolver, VectorHeatMethodSolver* vectorHeatSolver)
    : m_patchSpreadCoefficient(1.0) {
  m_mesh.reset(mesh);
  m_geometry.reset(geometry);
  m_mmpSolver.reset(mmpSolver);
  m_vectorHeatSolver.reset(vectorHeatSolver);

  m_axis = new SurfaceCurve(mesh, geometry, mmpSolver);
}

void SurfacePatch::createCustomAxis(std::vector<SurfacePoint>& axisPoints) { m_axis->createFromPoints(axisPoints); }

/*
 * A quick and dirty method for automatically generating an initial axis. Heuristic: Want the axis to roughly be aligned
 * with the longest dimension of the shape's bounding object. Computing bbox is kind of hard; just estimate the pair of
 * points with the greatest distance between them, and draw a geodesic between them.
 *
 * This might be kind of slow. See if I can come up with a better method.
 */
void SurfacePatch::createDefaultAxis() {

  // TODO: To save time, the Laplacian should only be computed once during the lifetime of the whole program.
  SurfacePoint pt1, pt2;
  double maxSoFar = 0;
  m_geometry->requireCotanLaplacian();

  SparseMatrix<double> C = m_geometry->cotanLaplacian;
  shiftDiagonal(C, 1e-8);

  PositiveDefiniteSolver<double> solver(C);

  Vector<double> RHS = Vector<double>::Zero(m_mesh->nVertices());
  Vector<double> distToSource; // actually just a quantity that is inversely proportional to distance

  for (size_t i = 0; i < m_points.size(); i++) {
    SurfacePoint pt = m_points[i];
    RHS[pt.vertex.getIndex()] = 1; // set
    distToSource = solver.solve(RHS);

    for (size_t j = 0; j < m_points.size(); j++) {
      if (j != i) {
        double dist = 1.0 / distToSource[m_points[j].vertex.getIndex()];
        if (dist > maxSoFar) {
          maxSoFar = dist;
          pt1 = pt;
          pt2 = m_points[j];
        }
      }
    }
    RHS[pt.vertex.getIndex()] = 0; // unset
  }

  std::unique_ptr<FlipEdgeNetwork> edgeNetwork;
  edgeNetwork = FlipEdgeNetwork::constructFromDijkstraPath(*m_mesh, *m_geometry, pt1.vertex, pt2.vertex);
  edgeNetwork->iterativeShorten();

  std::vector<std::vector<SurfacePoint>> paths = edgeNetwork->getPathPolyline();

  m_axis->setPoints(paths[0]);
}

void SurfacePatch::deformAxis(int index, std::complex<double> newDir) { m_axis->deformAngle(index, newDir); }

void SurfacePatch::getAxis(std::vector<SurfacePoint>& axis) { m_axis->getPoints(axis); }

void SurfacePatch::getPoints(std::vector<SurfacePoint>& points) { points = m_points; }

void SurfacePatch::getParameterizedAxis(std::vector<CurvePointParams>& parameterizedAxis) {
  m_axis->getParameterized(parameterizedAxis);
}

void SurfacePatch::getParameterizedPoints(std::vector<PatchPointParams>& parameterizedPoints) {
  parameterizedPoints = m_parameterizedPoints;
}

std::complex<double> SurfacePatch::getDirAtAxisIndex(int index) { return m_axis->getAngleAtIndex(index); }

Vector2 SurfacePatch::getInitDir() { return m_axis->getStartDir(); }

double SurfacePatch::getPatchSpreadCoefficient() { return m_patchSpreadCoefficient; }

void SurfacePatch::parameterizePatch(bool useFastParameterization, std::map<size_t, size_t> customClosestPointBinding) {
  m_parameterizedPoints.clear();

  std::cout << m_points.size() << std::endl;

  // Do closest point interpolation.

  bool closestPointBindingProvided = customClosestPointBinding.size() > 0;

  std::vector<std::tuple<SurfacePoint, double>> zippedDistances; // distances along curve
  VertexData<double> closestPointData;

  int axisSize = m_axis->getSize();

  // Default behavior - ignore if custom binding is provided
  if (!closestPointBindingProvided) {
    SurfacePoint sp = m_axis->getPointAtIndex(0);

    zippedDistances.push_back(std::make_tuple(sp, 0));

    double totalDist = 0;

    for (size_t i = 0; i < axisSize - 1; i++) {
      SurfacePoint pt2 = m_axis->getPointAtIndex(i + 1);
      double dist = m_axis->getDistanceAtIndex(i);
      totalDist += dist;
      zippedDistances.push_back(std::make_tuple(pt2, totalDist));
    }

    closestPointData = m_vectorHeatSolver->extendScalar(zippedDistances);
  }

  std::vector<PatchPointParams> ptToParam(m_points.size());

  std::complex<double> dir;
  double distance;
  size_t cp;

  for (size_t i = 0; i < m_points.size(); i++) {
    SurfacePoint patchPoint = m_points[i];

    if (!closestPointBindingProvided) {
      double diffusedVal = patchPoint.interpolate(closestPointData);
      cp = indexOfClosestPointOnAxis(diffusedVal, zippedDistances);
    } else {
      cp = customClosestPointBinding[i];
    }

    SurfacePoint axisPoint = m_axis->getPointAtIndex(cp);

    // Fast = logamp via VHM (accuracy increases with denser sampling, but poor on coarse meshes)
    // Slow = MMP (good for coarse meshes, but does not scale to higher sampling density)
    if (useFastParameterization) {
      VertexData<Vector2> logMap = m_vectorHeatSolver->computeLogMap(axisPoint);
      dir = evaluateLogMap(logMap, patchPoint);
      distance = std::abs(dir);
    } else {
      std::vector<SurfacePoint> geodesic = connectPointsWithGeodesic(axisPoint, patchPoint, distance);
      std::vector<SurfacePoint> prunedGeodesic = pruneApproxEqualEntries(geodesic);
      dir = localDir(axisPoint, prunedGeodesic[1]);
    }

    // Edge case: if distance is 0, then point is on axis
    // In which case just assign an arbitrary direction and let the distance
    // 0 out to prevent divide by 0 complaints
    if (distance <= 0) {
      std::cout << "WARNING (DON'T PANIC): point found to lie on axis" << std::endl;
      dir = std::complex<double>{1.0, 0.0};
    }

    // Compute relative to the tangent direction of the axis
    dir /= std::abs(dir);

    std::complex<double> axisBasis;

    axisBasis = m_axis->getTangentAtIndex(cp);
    axisBasis /= std::abs(axisBasis);

    dir = dir / axisBasis;

    PatchPointParams prms = {cp, distance, dir};
    ptToParam[i] = prms;
  }

  m_parameterizedPoints.insert(m_parameterizedPoints.begin(), ptToParam.begin(), ptToParam.end());
}

void SurfacePatch::reconstructPatch() {
  m_points.clear();

  // Use distances/directions to reconstruct patches from sparse axis
  SurfacePoint pathEndpoint;
  TraceGeodesicResult tracedGeodesic;

  std::vector<SurfacePoint> constructedPatch(m_parameterizedPoints.size());

  for (size_t i = 0; i < m_parameterizedPoints.size(); i++) {
    PatchPointParams p = m_parameterizedPoints[i];
    std::complex<double> dir = p.axisPointDirection;
    double distance = m_patchSpreadCoefficient * p.axisPointDistance;
    SurfacePoint startPoint = m_axis->getPointAtIndex(p.closestAxisPoint);

    dir /= std::abs(dir);

    std::complex<double> axisBasis;

    axisBasis = m_axis->getTangentAtIndex(p.closestAxisPoint);
    axisBasis /= std::abs(axisBasis);

    dir *= axisBasis;

    // Handle distance = 0 edge cases
    if (distance <= 0) {
      pathEndpoint = startPoint;
    } else {
      tracedGeodesic = traceGeodesic(*(m_geometry), m_axis->getPointAtIndex(p.closestAxisPoint),
                                     Vector2::fromComplex(distance * dir));
      pathEndpoint = tracedGeodesic.endPoint;
    }

    constructedPatch[i] = pathEndpoint;
  }

  m_points.insert(m_points.begin(), constructedPatch.begin(), constructedPatch.end());
}

void SurfacePatch::rotateAxis(Vector2 newDir) { m_axis->rotate(newDir); }

void SurfacePatch::saveAxisStartPointAndDirection() { m_axis->saveStartParams(); }

void SurfacePatch::setPatchSpreadCoefficient(double patchSpreadCoefficient) {
  m_patchSpreadCoefficient = patchSpreadCoefficient;
  reconstructPatch();
}

void SurfacePatch::transfer(SurfacePatch* target, const SurfacePoint& targetMeshStart,
                            const SurfacePoint& targetMeshDirEndpoint) {
  // First transfer the axis

  m_axis->transfer(target->m_axis, targetMeshStart, targetMeshDirEndpoint);

  // Compute distances and directions on S1, then reconstruct contact on S2
  // Reconstruct patch only if patch has been parameterized on source domain

  if (m_parameterizedPoints.size() > 0) {
    target->m_parameterizedPoints.clear();

    for (int i = 0; i < m_parameterizedPoints.size(); i++) {
      PatchPointParams sourceParams = m_parameterizedPoints[i];

      std::complex<double> sourceDirection = sourceParams.axisPointDirection;

      std::complex<double> adjustedDir = {sourceDirection.real(), -sourceDirection.imag()};
      PatchPointParams targetParams = {sourceParams.closestAxisPoint, sourceParams.axisPointDistance, adjustedDir};
      target->m_parameterizedPoints.push_back(targetParams);
    }

    target->reconstructPatch();
  }
}

void SurfacePatch::translate(const SurfacePoint& newStartPoint) { m_axis->translate(newStartPoint); }

void SurfacePatch::setAxisPoints(std::vector<SurfacePoint>& axisPoints) { m_axis->setPoints(axisPoints); }

void SurfacePatch::setPatchPoints(const std::vector<SurfacePoint>& points) { m_points = points; }

void SurfacePatch::setParameterizedAxis(std::vector<CurvePointParams>& parameterizedAxisPoints) {
  m_axis->setParameterized(parameterizedAxisPoints);
}

void SurfacePatch::setParameterizedPatch(std::vector<PatchPointParams>& parameterizedPatchPoints) {
  m_parameterizedPoints.clear();
  m_parameterizedPoints.insert(m_parameterizedPoints.begin(), parameterizedPatchPoints.begin(),
                               parameterizedPatchPoints.end());
}

// Begin private utils

bool SurfacePatch::approxEqual(const SurfacePoint& pA, const SurfacePoint& pB) {

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

/*
 * Given two SurfacePoints, connect them with a geodesic if needed (if they aren't in the same face.)
 * Return a vector of SurfacePoints making up this geodesic, with the endpoints excluded.
 */

std::vector<SurfacePoint> SurfacePatch::connectPointsWithGeodesic(const SurfacePoint& pt1, const SurfacePoint& pt2,
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

Vector2 SurfacePatch::evaluateLogMap(const VertexData<Vector2>& logMap, const SurfacePoint& pt) {

  if (pt.type == SurfacePointType::Vertex) {
    return logMap[pt.vertex];
  }
  if (pt.type == SurfacePointType::Edge) {
    Edge e = pt.edge;
    double t = pt.tEdge;
    return (1.0 - t) * logMap[e.firstVertex()] + t * logMap[e.secondVertex()];
  }
  // face
  Face f = pt.face;
  Vector3 fCoords = pt.faceCoords;
  Vector2 result = {0, 0};
  size_t i = 0;
  for (Vertex v : f.adjacentVertices()) {
    result += fCoords[i] * logMap[v];
    i++;
  }
  return result;
}

/*
 * Given a value at a SurfacePoint and the initial values of points in the axis, determine the closest point on the
 * axis.
 */
size_t SurfacePatch::indexOfClosestPointOnAxis(double diffusedVal,
                                               const std::vector<std::tuple<SurfacePoint, double>>& zippedDistances) {

  // zippedDistances should be the same size as the axis
  for (size_t i = 0; i < zippedDistances.size() - 1; i++) {
    double dist0 = std::get<1>(zippedDistances[i]);
    double dist1 = std::get<1>(zippedDistances[i + 1]);
    if (diffusedVal >= dist0 && diffusedVal <= dist1) {
      if ((diffusedVal - dist0) > (dist1 - diffusedVal)) {
        return i;
      }
      return i + 1;
    }
  }
  assert(false); // should never get here

  return -1;
}

/*
 * Given two SurfacePoints assumed to lie within the same face, determine the direction from pt1 to pt2, in the local
 * tangent basis of pt1.
 */
Vector2 SurfacePatch::localDir(const SurfacePoint& pt1, const SurfacePoint& pt2) {

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

  // TODO: Need to handle boundary cases (if vertex is on the boundary, etc.)

  Vector3 local_x = m_geometry->halfedgeVector(heDir).normalize();
  Vector3 local_y = cross(normal, local_x);

  Vector2 dir = {dot(globalDir, local_x), dot(globalDir, local_y)}; // not necessarily unit
  dir /= dir.norm();
  return dir;
}

std::vector<SurfacePoint> SurfacePatch::pruneApproxEqualEntries(const std::vector<SurfacePoint>& source) {
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
