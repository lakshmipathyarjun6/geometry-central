#include "geometrycentral/surface/surface_patch.h"
#include "geometrycentral/surface/surface_point.h"

int mod(int a, int b) { return (b + (a % b)) % b; }

SurfacePatch::SurfacePatch(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry,
                           HeatMethodDistanceSolver* distanceHeatSolver, VectorHeatMethodSolver* vectorHeatSolver) {
  m_mesh.reset(mesh);
  m_geometry.reset(geometry);
  m_distanceHeatSolver.reset(distanceHeatSolver);
  m_vectorHeatSolver.reset(vectorHeatSolver);
}

/*
 * To support smooth translation, parallel transport the base initial direction to all vertices.
 * Table can be referenced during drag, but should be recomputed at start.
 */
void SurfacePatch::computeGlobalAxisTransport() {
  Vertex startVertex = m_startPoint.nearestVertex();
  m_axisDirectionTable = m_vectorHeatSolver->transportTangentVector(startVertex, m_initDir);
}

/*
 * The start direction should be determined using the global direction between the first two points
 * but in the tangent basis of the first point.  It can subsequently be overridden by user input via rotation.
 */
void SurfacePatch::computeInitialAxisDirection() {
  m_geometry->requireVertexNormals();
  m_geometry->requireVertexTangentBasis();

  SurfacePoint firstPoint = m_patchAxisSparse[0];
  SurfacePoint secondPoint = m_patchAxisSparse[1];

  m_initDir = localDir(firstPoint, secondPoint);
  m_initDir /= m_initDir.norm();
}

void SurfacePatch::createCustomAxis(std::vector<Vertex>& axisPoits) {
  m_patchAxisSparse.clear();

  std::unique_ptr<FlipEdgeNetwork> edgeNetwork;
  edgeNetwork = FlipEdgeNetwork::constructFromPiecewiseDijkstraPath(*m_mesh, *m_geometry, axisPoits);

  std::vector<std::vector<SurfacePoint>> paths = edgeNetwork->getPathPolyline();

  m_patchAxisSparse = paths[0];
  m_startPoint = m_patchAxisSparse[0];
  computeInitialAxisDirection();
  constructDenselySampledAxis();
  computeAxisAnglesAndDistances();
}

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

  for (size_t i = 0; i < m_patchBoundary.size(); i++) {
    SurfacePoint pt = m_patchBoundary[i];
    RHS[pt.vertex.getIndex()] = 1; // set
    distToSource = solver.solve(RHS);

    for (size_t j = 0; j < m_patchBoundary.size(); j++) {
      if (j != i) {
        double dist = 1.0 / distToSource[m_patchBoundary[j].vertex.getIndex()];
        if (dist > maxSoFar) {
          maxSoFar = dist;
          pt1 = pt;
          pt2 = m_patchBoundary[j];
        }
      }
    }
    RHS[pt.vertex.getIndex()] = 0; // unset
  }

  std::unique_ptr<FlipEdgeNetwork> edgeNetwork;
  edgeNetwork = FlipEdgeNetwork::constructFromDijkstraPath(*m_mesh, *m_geometry, pt1.vertex, pt2.vertex);

  std::vector<std::vector<SurfacePoint>> paths = edgeNetwork->getPathPolyline();

  m_patchAxisSparse = paths[0];
  m_startPoint = m_patchAxisSparse[0];
  computeInitialAxisDirection();
  constructDenselySampledAxis();
  computeAxisAnglesAndDistances();
}

void SurfacePatch::get(std::vector<SurfacePoint>& axis, std::vector<SurfacePoint>& boundary) {
  axis = m_patchAxisSparse;
  boundary = m_patchBoundary;
}

std::vector<std::string> SurfacePatch::getAxisSerialized() {

  std::vector<std::string> result;

  int numElements = m_patchAxisSparse.size();

  for (int i = 0; i < numElements; i++) {
    SurfacePoint p = m_patchAxisSparse[i];
    std::complex<double> dir = m_patchAxisSparseAngles[i];
    double dist = m_patchAxisSparseDistances[i];

    SurfacePointType t = p.type;

    std::string serializedExtraParams =
        std::to_string(dir.real()) + " " + std::to_string(dir.imag()) + " " + std::to_string(dist);

    std::string serialized;

    if (t == SurfacePointType::Edge) {
      Edge e = p.edge;
      int idx = e.getIndex();
      double tEdge = p.tEdge;

      serialized = "e " + std::to_string(idx) + " " + std::to_string(tEdge) + " " + serializedExtraParams;
    } else if (t == SurfacePointType::Face) {
      Face f = p.face;
      int idx = f.getIndex();
      Vector3 faceCoords = p.faceCoords;

      serialized = "f " + std::to_string(idx) + " " + std::to_string(faceCoords.x) + " " +
                   std::to_string(faceCoords.y) + " " + std::to_string(faceCoords.z) + " " + serializedExtraParams;
    } else if (t == SurfacePointType::Vertex) {
      Vertex v = p.vertex;
      int idx = v.getIndex();

      serialized = "v " + std::to_string(idx) + " " + serializedExtraParams;
    } else {
      std::cout << "Unknown surface type found. Killing" << std::endl;
      serialized = "NULL";
    }

    result.push_back(serialized);
  }

  return result;
}

std::vector<std::string> SurfacePatch::getBoundarySerialized() { return getSerializedSurfacePoints(m_patchBoundary); }

Vector2 SurfacePatch::getInitDir() { return m_initDir; }

void SurfacePatch::invertAxisOrder() {
  reverse(m_patchAxisSparse.begin(), m_patchAxisSparse.end());

  m_startPoint = m_patchAxisSparse[0];
  computeInitialAxisDirection();
  constructDenselySampledAxis();
  computeAxisAnglesAndDistances();
}

void SurfacePatch::invertBoundaryOrder() { reverse(m_parameterizedBoundary.begin(), m_parameterizedBoundary.end()); }

void SurfacePatch::leftShiftOrder() {
  params firstElem = m_parameterizedBoundary[0];
  m_parameterizedBoundary.push_back(firstElem);
  m_parameterizedBoundary.erase(m_parameterizedBoundary.begin());
}

void SurfacePatch::linkPatch(std::string childName, SurfacePatch* child) {
  m_children[childName] = child;
  child->m_parent = this;

  Vertex parentStartpoint = m_startPoint.nearestVertex();
  Vertex childStartpoint = child->m_startPoint.nearestVertex();
  Vertex childSecondpoint = child->m_patchAxisSparse[1].nearestVertex();

  TraceGeodesicResult tracedGeodesic;
  TraceOptions traceOptions;
  traceOptions.includePath = true;

  std::unique_ptr<FlipEdgeNetwork> edgeNetwork;
  edgeNetwork = FlipEdgeNetwork::constructFromDijkstraPath(*m_mesh, *m_geometry, parentStartpoint, childStartpoint);

  std::vector<std::vector<SurfacePoint>> paths = edgeNetwork->getPathPolyline();
  std::vector<SurfacePoint> result = paths[0];

  double totalDist1 = 0;

  for (size_t i = 0; i < result.size() - 1; i++) {
    SurfacePoint pt1 = result[i];
    SurfacePoint pt2 = result[i + 1];
    Vector3 a = pt1.interpolate(m_geometry->inputVertexPositions);
    Vector3 b = pt2.interpolate(m_geometry->inputVertexPositions);
    totalDist1 += (a - b).norm();
  }

  edgeNetwork = FlipEdgeNetwork::constructFromDijkstraPath(*m_mesh, *m_geometry, parentStartpoint, childSecondpoint);

  paths = edgeNetwork->getPathPolyline();
  result = paths[0];

  double totalDist2 = 0;

  for (size_t i = 0; i < result.size() - 1; i++) {
    SurfacePoint pt1 = result[i];
    SurfacePoint pt2 = result[i + 1];
    Vector3 a = pt1.interpolate(m_geometry->inputVertexPositions);
    Vector3 b = pt2.interpolate(m_geometry->inputVertexPositions);
    totalDist2 += (a - b).norm();
  }

  Vector2 offsetDepartDir = m_initDir;
  Vector2 targetDepartDir1 = localDir(parentStartpoint, childStartpoint);
  Vector2 targetDepartDir2 = localDir(parentStartpoint, child->m_patchAxisSparse[1].nearestVertex());

  std::complex<double> localAngleDepart1 = targetDepartDir1 / offsetDepartDir;
  std::complex<double> localAngleDepart2 = targetDepartDir2 / offsetDepartDir;

  std::tuple<std::complex<double>, std::complex<double>, double, double> childTraceParams =
      std::make_tuple(localAngleDepart1, localAngleDepart2, totalDist1, totalDist2);

  m_childTraceParams[childName] = childTraceParams;
}

void SurfacePatch::propagateChildUpdates() {
  SurfacePoint pathEndpoint;

  TraceGeodesicResult tracedGeodesic;
  TraceOptions traceOptions;
  traceOptions.includePath = true;

  for (auto const& child : m_children) {
    std::string childName = child.first;
    SurfacePatch* childPatch = child.second;

    std::tuple<std::complex<double>, std::complex<double>, double, double> childTraceParams =
        m_childTraceParams[childName];

    std::complex<double> localAngleDepart1 = std::get<0>(childTraceParams);
    std::complex<double> localAngleDepart2 = std::get<1>(childTraceParams);
    double childLocalDist1 = std::get<2>(childTraceParams);
    double childLocalDist2 = std::get<3>(childTraceParams);

    std::complex<double> offsetDepartDir = m_initDir;
    std::complex<double> targetDepartDir1 = offsetDepartDir * localAngleDepart1;

    tracedGeodesic = traceGeodesic(*(m_geometry), m_startPoint,
                                   Vector2::fromComplex(targetDepartDir1 * childLocalDist1), traceOptions);
    pathEndpoint = tracedGeodesic.endPoint;

    childPatch->m_patchAxisSparse[0] = pathEndpoint;
    childPatch->m_startPoint = childPatch->m_patchAxisSparse[0];

    std::complex<double> targetDepartDir2 = offsetDepartDir * localAngleDepart2;

    tracedGeodesic = traceGeodesic(*(m_geometry), m_startPoint,
                                   Vector2::fromComplex(targetDepartDir2 * childLocalDist2), traceOptions);
    pathEndpoint = tracedGeodesic.endPoint;
    childPatch->m_patchAxisSparse[1] = pathEndpoint;

    childPatch->m_initDir = localDir(childPatch->m_patchAxisSparse[0], childPatch->m_patchAxisSparse[1]);

    childPatch->traceAxis();

    childPatch->propagateChildUpdates();
  }
}

void SurfacePatch::reconstructBoundary() {
  m_patchBoundary.clear();

  // Use distances/directions to reconstruct patches from sparse axis
  SurfacePoint pathEndpoint;
  TraceGeodesicResult tracedGeodesic;

  std::vector<SurfacePoint> constructedBoundary(m_parameterizedBoundary.size());

  for (size_t i = 0; i < m_parameterizedBoundary.size(); i++) {
    params p = m_parameterizedBoundary[i];
    std::complex<double> dir = p.dir;
    std::complex<double> axisBasis = axisTangent(m_patchAxisSparseDenseIdx[p.cp], m_patchAxisDense);
    dir /= std::abs(dir);
    axisBasis /= std::abs(axisBasis);
    dir *= axisBasis;
    tracedGeodesic = traceGeodesic(*(m_geometry), m_patchAxisSparse[p.cp], Vector2::fromComplex(p.dist * dir));
    pathEndpoint = tracedGeodesic.endPoint;
    constructedBoundary[i] = pathEndpoint;
  }

  m_patchBoundary.insert(m_patchBoundary.begin(), constructedBoundary.begin(), constructedBoundary.end());
}

void SurfacePatch::reparameterizeBoundary() {
  m_parameterizedBoundary.clear();

  VertexData<double> distToSource = m_distanceHeatSolver->computeDistance(m_patchAxisSparse);

  std::cout << m_patchBoundary.size() << std::endl;

  // Do closest point interpolation.

  std::vector<std::tuple<SurfacePoint, double>> zippedDistances; // distances along curve
  zippedDistances.push_back(std::make_tuple(m_startPoint, 0));
  double totalDist = 0;
  for (size_t i = 0; i < m_patchAxisSparse.size() - 1; i++) {
    SurfacePoint pt2 = m_patchAxisSparse[i + 1];
    double dist = m_patchAxisSparseDistances[i];
    totalDist += dist;
    zippedDistances.push_back(std::make_tuple(pt2, totalDist));
  }
  VertexData<double> closestPoint = m_vectorHeatSolver->extendScalar(zippedDistances);

  std::vector<params> bdyPtToParam(m_patchBoundary.size());

  for (size_t i = 0; i < m_patchBoundary.size(); i++) {
    SurfacePoint bdyPoint = m_patchBoundary[i];
    double diffusedVal = evaluateVertexDataAtPoint(closestPoint, bdyPoint);
    size_t cp = indexOfClosestPointOnAxis(diffusedVal, zippedDistances);
    SurfacePoint axisPoint = m_patchAxisSparse[cp];
    VertexData<Vector2> logMap = m_vectorHeatSolver->computeLogMap(axisPoint);
    std::complex<double> dir = evaluateLogMap(logMap, bdyPoint);
    double logMapDist = std::abs(dir);

    // Edge case: if distance is 0, then point is on axis
    // In which case just assign an arbitrary direction and let the distance
    // 0 out to prevent divide by 0 complaints
    if (logMapDist <= 0) {
      std::cout << "WARNING (DON'T PANIC): point found to lie on axis" << std::endl;
      dir = std::complex<double>{1.0, 0.0};
    }

    // Compute relative to the tangent direction of the axis
    dir /= std::abs(dir);
    std::complex<double> axisBasis = axisTangent(m_patchAxisSparseDenseIdx[cp], m_patchAxisDense);
    axisBasis /= std::abs(axisBasis);
    dir = dir / axisBasis;
    double heatDist = evaluateVertexDataAtPoint(distToSource, bdyPoint);
    // TODO: Also try log map dist
    params prms = {cp, heatDist, -dir};
    bdyPtToParam[i] = prms;
  }

  m_parameterizedBoundary.insert(m_parameterizedBoundary.begin(), bdyPtToParam.begin(), bdyPtToParam.end());
}

void SurfacePatch::rightShiftOrder() {
  params lastElem = m_parameterizedBoundary[m_parameterizedBoundary.size() - 1];
  m_parameterizedBoundary.insert(m_parameterizedBoundary.begin(), lastElem);
  m_parameterizedBoundary.pop_back();
}

void SurfacePatch::rotateAxis(Vector2 newDir) {
  m_initDir = newDir;
  traceAxis();
  propagateChildUpdates();
}

void SurfacePatch::setBulkTransferParams(SurfacePatch* sourcePatch, std::string sourcePatchName,
                                         std::string destinationPatchName) {
  std::cout << "Overriding child trace params of " << destinationPatchName << " with " << sourcePatchName << std::endl;

  std::tuple<std::complex<double>, std::complex<double>, double, double> sourceChildTraceParams =
      sourcePatch->m_childTraceParams[sourcePatchName];

  std::complex<double> localAngleDepart1 = std::get<0>(sourceChildTraceParams);
  std::complex<double> localAngleDepart2 = std::get<1>(sourceChildTraceParams);
  double sourcedLocalDist1 = std::get<2>(sourceChildTraceParams);
  double sourceLocalDist2 = std::get<3>(sourceChildTraceParams);

  std::complex<double> invertedLocalAngleDepart1 =
      std::complex<double>{localAngleDepart1.real(), -localAngleDepart1.imag()};
  std::complex<double> invertedLocalAngleDepart2 =
      std::complex<double>{localAngleDepart2.real(), -localAngleDepart2.imag()};

  m_childTraceParams[destinationPatchName] =
      std::make_tuple(invertedLocalAngleDepart1, invertedLocalAngleDepart2, sourcedLocalDist1, sourceLocalDist2);

  propagateChildUpdates();
}

void SurfacePatch::transfer(SurfacePatch* target, const Vertex& targetMeshStart, const Vertex& targetMeshDirEndpoit) {
  // All angles and distances should be the same, save for the first one (which is ignored)
  target->m_patchAxisSparseAngles = m_patchAxisSparseAngles;
  target->m_patchAxisSparseDistances.insert(target->m_patchAxisSparseDistances.end(),
                                            m_patchAxisSparseDistances.begin(), m_patchAxisSparseDistances.end());

  SurfacePoint tmStart = SurfacePoint(targetMeshStart);
  SurfacePoint tmEnd = SurfacePoint(targetMeshDirEndpoit);

  target->m_startPoint = tmStart;
  target->m_initDir = target->localDir(tmStart, tmEnd);
  target->traceAxis();

  // Compute distances and directions on S1, then reconstruct contact on S2

  // Need to also invert boundary order (unintuitive)
  std::vector<params> targetParameterizedBoundary(m_parameterizedBoundary.size());

  for (int i = 0; i < m_parameterizedBoundary.size(); i++) {
    params sourceParams = m_parameterizedBoundary[i];
    params targetParams = {sourceParams.cp, sourceParams.dist, -sourceParams.dir};
    targetParameterizedBoundary[i] = targetParams;
  }

  target->m_parameterizedBoundary = targetParameterizedBoundary;
  target->reconstructBoundary();
}

void SurfacePatch::translate(const Vertex& newStartVertex) {
  m_startPoint = SurfacePoint(newStartVertex);
  m_initDir = m_axisDirectionTable[newStartVertex];
  traceAxis();
  propagateChildUpdates();
}

void SurfacePatch::setPatchAxis(const std::vector<SurfacePoint>& axis, const std::vector<std::complex<double>>& dirs,
                                const std::vector<double>& dists) {
  m_patchAxisSparse = axis;
  m_startPoint = m_patchAxisSparse[0];
  m_patchAxisSparseAngles = dirs;
  m_patchAxisSparseDistances = dists;

  computeInitialAxisDirection();
  traceAxis();
}

void SurfacePatch::setPatchBoundary(const std::vector<SurfacePoint>& boundary) { m_patchBoundary = boundary; }

void SurfacePatch::unlinkAllPatches() {
  std::vector<std::string> allChildren;

  for (auto const& child : m_children) {
    std::string childName = child.first;
    allChildren.push_back(childName);
  }

  for (std::string childName : allChildren) {
    SurfacePatch* childPatch = m_children[childName];
    childPatch->unlinkAllPatches();
    unlinkPatch(childName);
  }
}

void SurfacePatch::unlinkPatch(std::string childName) {
  if (m_childTraceParams.count(childName) <= 0) {
    std::cout << "Child does not exist - doing nothing" << std::endl;
    return;
  }

  SurfacePatch* child = m_children[childName];
  child->m_parent = NULL;

  m_children.erase(childName);
  m_childTraceParams.erase(childName);
}

// Begin private utils

/*
 * Compute the direction of the basis of the axis curve, at the given point.
 */
std::complex<double> SurfacePatch::axisTangent(size_t idx, const std::vector<SurfacePoint>& axis) {

  int nextDir = (idx == axis.size() - 1) ? -1 : 1;
  SurfacePoint pt1 = axis[idx];
  SurfacePoint pt2 = axis[idx + nextDir];
  std::complex<double> ref = localDir(pt1, pt2);
  return ref;
}

/*
 * Given a series of SurfacePoints (representing a curve) on a mesh, compute a complex number z representing the
 * "angle" at each point of the polygon. The argument std::arg(z) represents the rotation from the negative of the
 * tangent vector of the current edge to the tangent vector of the next edge; the reason for doing this this way is that
 * the local directions need to be computed at a common point (in the same local tangent space), i.e. at each point we
 * consider the two *outgoing* vectors along its adjacent edges. (The turning angle is sgn(angle)*PI - angle.) Also for
 * each point, computes the distance to the next point.
 *
 * <sparseCurve> is the original curve. The angles and distances returned will correspond to this curve.
 * <denseCurve> is assumed to be densely sampled, s.t. angles & distances are easily computed without using the log map.
 * <idxIntoDense> indicates how sparse points get mapped into the dense vector.
 *
 * Warning: If the curve is open (like if the curve represents an axis), then the angles at the endpoints will be
 * garbage (as well as the distance value at the last point) -- this isn't a problem since we never use those values
 * later on, but just to be aware.
 */
void SurfacePatch::computeAxisAnglesAndDistances() {

  size_t N = m_patchAxisSparse.size(); // # of vertices in the polygon
  size_t M = m_patchAxisDense.size();

  std::vector<std::complex<double>> angles(N);
  std::vector<double> distances(N);

  SurfacePoint currNode, prevNode, nextNode;
  std::complex<double> origDir, offsetDir, adjustedDir;
  Vector3 currPos, nextPos;
  size_t currIdx, nextIdx, prevIdx;
  for (size_t i = 0; i < N; i++) {
    currIdx = m_patchAxisSparseDenseIdx[i];
    nextIdx = mod(currIdx + 1, M);
    prevIdx = mod(currIdx - 1, M);
    currNode = m_patchAxisDense[currIdx];
    nextNode = m_patchAxisDense[nextIdx];
    prevNode = m_patchAxisDense[prevIdx];

    double totalDist = 0;
    SurfacePoint pt1, pt2;

    for (size_t j = currIdx; j < m_patchAxisSparseDenseIdx[mod(i + 1, N)]; j++) {
      pt1 = m_patchAxisDense[j];
      pt2 = m_patchAxisDense[j + 1];
      Vector3 a = pt1.interpolate(m_geometry->inputVertexPositions);
      Vector3 b = pt2.interpolate(m_geometry->inputVertexPositions);
      totalDist += (a - b).norm();
    }
    distances[i] = totalDist;

    origDir = localDir(currNode, nextNode);
    offsetDir = localDir(currNode, prevNode);

    // CCW angle of rotation to get to next tangent vector; for convex corners, rotation is CCW (and angle is
    // positive), CW for non-convex corners (angle is negative).
    adjustedDir = origDir / offsetDir;
    angles[i] = adjustedDir;
  }

  m_patchAxisSparseAngles = angles;
  m_patchAxisSparseDistances = distances;
}

/*
 * Given two SurfacePoints, connect them with a geodesic if needed (if they aren't in the same face.)
 * Return a vector of SurfacePoints making up this geodesic, with the endpoints excluded.
 */
std::vector<SurfacePoint> SurfacePatch::connectPointsWithGeodesic(const SurfacePoint& pt1, const SurfacePoint& pt2) {

  if (sharedFace(pt1, pt2) != Face()) return std::vector<SurfacePoint>();
  assert(pt1.type == SurfacePointType::Vertex);
  assert(pt2.type == SurfacePointType::Vertex);
  std::unique_ptr<FlipEdgeNetwork> edgeNetwork;
  edgeNetwork = FlipEdgeNetwork::constructFromDijkstraPath(*m_mesh, *m_geometry, pt1.vertex, pt2.vertex);
  edgeNetwork->iterativeShorten();
  std::vector<std::vector<SurfacePoint>> paths = edgeNetwork->getPathPolyline();
  assert(paths.size() == 1);
  std::vector<SurfacePoint> fullPath = paths[0];
  assert(fullPath.size() > 2);
  std::vector<SurfacePoint> path(fullPath.begin() + 1, fullPath.begin() + fullPath.size() - 1);
  return path;
}

/*
 * Given a curve (possibly closed), fill in any gaps with geodesics.
 * Returns a vector of SurfacePoints with no gaps (first argument of the tuple.)
 * Also returns a map <idxIntoDense> (2nd argument of the tuple) where the ith entry is the index of the ith point in
 * the original curve in the dense curve vector.
 */
void SurfacePatch::constructDenselySampledAxis() {
  m_patchAxisDense.clear();
  m_patchAxisSparseDenseIdx.clear();

  std::vector<SurfacePoint> denseCurve;
  std::vector<size_t> idxIntoDense(m_patchAxisSparse.size());

  size_t currIdx = 0;
  size_t N = m_patchAxisSparse.size();

  SurfacePoint pt1, pt2;

  for (size_t i = 0; i < N; i++) {
    pt1 = m_patchAxisSparse[i];
    pt2 = m_patchAxisSparse[mod(i + 1, N)];
    denseCurve.push_back(pt1);
    idxIntoDense[i] = currIdx;
    currIdx++;

    std::vector<SurfacePoint> path = connectPointsWithGeodesic(pt1, pt2);
    denseCurve.insert(denseCurve.end(), path.begin(), path.end());
    currIdx += path.size();
  }

  m_patchAxisDense = denseCurve;
  m_patchAxisSparseDenseIdx = idxIntoDense;
}

/*
 * Evaluate a log map (given at vertices) at any point.
 * Just linearly interpolate.
 */
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
 * Given a function defined as data on vertices, evaluate this function at the given SurfacePoint.
 */
double SurfacePatch::evaluateVertexDataAtPoint(const VertexData<double>& u, const SurfacePoint& pt) {

  if (pt.type == SurfacePointType::Vertex) {
    return u[pt.vertex];
  }
  if (pt.type == SurfacePointType::Edge) {
    Edge e = pt.edge;
    double a = u[e.firstVertex()];
    double b = u[e.secondVertex()];
    double t = pt.tEdge;
    return (1.0 - t) * a + t * b;
  }
  if (pt.type == SurfacePointType::Face) {
    std::vector<double> vals;
    Vector3 faceCoords = pt.faceCoords;
    Face f = pt.face;
    for (Vertex v : f.adjacentVertices()) {
      vals.push_back(u[v]);
    }
    double val = 0;
    for (size_t i = 0; i < vals.size(); i++) {
      val += faceCoords[i] * vals[i];
    }
    return val;
  }
  return -1; // shouldn't get here
}

std::vector<std::string> SurfacePatch::getSerializedSurfacePoints(std::vector<SurfacePoint>& source) {
  std::vector<std::string> result;

  for (SurfacePoint p : source) {
    SurfacePointType t = p.type;

    std::string serialized;

    if (t == SurfacePointType::Edge) {
      Edge e = p.edge;
      int idx = e.getIndex();
      double tEdge = p.tEdge;

      serialized = "e " + std::to_string(idx) + " " + std::to_string(tEdge);
    } else if (t == SurfacePointType::Face) {
      Face f = p.face;
      int idx = f.getIndex();
      Vector3 faceCoords = p.faceCoords;

      serialized = "f " + std::to_string(idx) + " " + std::to_string(faceCoords.x) + " " +
                   std::to_string(faceCoords.y) + " " + std::to_string(faceCoords.z);
    } else if (t == SurfacePointType::Vertex) {
      Vertex v = p.vertex;
      int idx = v.getIndex();

      serialized = "v " + std::to_string(idx);
    } else {
      std::cout << "Unknown surface type found. Killing" << std::endl;
      serialized = "NULL";
    }

    result.push_back(serialized);
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

void SurfacePatch::traceAxis() {
  // Clear all existing data
  m_patchAxisSparse.clear();
  m_patchAxisSparseDenseIdx.clear();
  m_patchAxisDense.clear();

  // Place first points
  m_patchAxisSparse.emplace_back(m_startPoint);
  m_patchAxisDense.emplace_back(m_startPoint);
  m_patchAxisSparseDenseIdx.push_back(0);

  // Declare some variables
  size_t N = m_patchAxisSparseAngles.size();
  SurfacePoint pathEndpoint;
  std::complex<double> endingDir;
  TraceGeodesicResult tracedGeodesic;
  TraceOptions traceOptions;
  traceOptions.includePath = true;

  // Trace out from first to second point
  m_initDir /= m_initDir.norm();
  tracedGeodesic = traceGeodesic(*(m_geometry), m_startPoint,
                                 Vector2::fromComplex(m_initDir * m_patchAxisSparseDistances[0]), traceOptions);

  pathEndpoint = tracedGeodesic.endPoint;
  m_patchAxisDense.insert(m_patchAxisDense.end(), tracedGeodesic.pathPoints.begin() + 1,
                          tracedGeodesic.pathPoints.end());
  m_patchAxisSparse.push_back(pathEndpoint);
  m_patchAxisSparseDenseIdx.push_back(m_patchAxisDense.size() - 1);
  endingDir = -tracedGeodesic.endingDir;

  // Trace the rest
  for (size_t i = 1; i < N - 1; i++) {
    // Get the corresponding direction on S2
    endingDir /= std::abs(endingDir);
    std::complex<double> adjustedDirS2 = endingDir * m_patchAxisSparseAngles[i];
    adjustedDirS2 /= std::abs(adjustedDirS2); // make sure direction is unit

    // Trace geodesic from last point and record where it ended up
    double dist = m_patchAxisSparseDistances[i];
    tracedGeodesic = traceGeodesic(*(m_geometry), m_patchAxisSparse[m_patchAxisSparse.size() - 1],
                                   Vector2::fromComplex(adjustedDirS2 * dist), traceOptions);

    pathEndpoint = tracedGeodesic.endPoint;
    m_patchAxisDense.insert(m_patchAxisDense.end(), tracedGeodesic.pathPoints.begin() + 1,
                            tracedGeodesic.pathPoints.end());
    m_patchAxisSparse.push_back(pathEndpoint);
    m_patchAxisSparseDenseIdx.push_back(m_patchAxisDense.size() - 1);
    endingDir = -tracedGeodesic.endingDir;
  }
}
