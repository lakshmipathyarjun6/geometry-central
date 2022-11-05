#include "geometrycentral/surface/surface_patch.h"

int mod(int a, int b) { return (b + (a % b)) % b; }

SurfacePatch::SurfacePatch(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry,
                           GeodesicAlgorithmExact* mmpSolver, HeatMethodDistanceSolver* distanceHeatSolver,
                           VectorHeatMethodSolver* vectorHeatSolver) {
  m_mesh.reset(mesh);
  m_geometry.reset(geometry);
  m_distanceHeatSolver.reset(distanceHeatSolver);
  m_mmpSolver.reset(mmpSolver);
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

void SurfacePatch::createCustomAxis(std::vector<Vertex>& axisPoints) {
  m_patchAxisSparse.clear();

  std::unique_ptr<FlipEdgeNetwork> edgeNetwork;
  edgeNetwork = FlipEdgeNetwork::constructFromPiecewiseDijkstraPath(*m_mesh, *m_geometry, axisPoints);
  edgeNetwork->iterativeShorten();

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

  for (size_t i = 0; i < m_patchPoints.size(); i++) {
    SurfacePoint pt = m_patchPoints[i];
    RHS[pt.vertex.getIndex()] = 1; // set
    distToSource = solver.solve(RHS);

    for (size_t j = 0; j < m_patchPoints.size(); j++) {
      if (j != i) {
        double dist = 1.0 / distToSource[m_patchPoints[j].vertex.getIndex()];
        if (dist > maxSoFar) {
          maxSoFar = dist;
          pt1 = pt;
          pt2 = m_patchPoints[j];
        }
      }
    }
    RHS[pt.vertex.getIndex()] = 0; // unset
  }

  std::unique_ptr<FlipEdgeNetwork> edgeNetwork;
  edgeNetwork = FlipEdgeNetwork::constructFromDijkstraPath(*m_mesh, *m_geometry, pt1.vertex, pt2.vertex);
  edgeNetwork->iterativeShorten();

  std::vector<std::vector<SurfacePoint>> paths = edgeNetwork->getPathPolyline();

  m_patchAxisSparse = paths[0];
  m_startPoint = m_patchAxisSparse[0];
  computeInitialAxisDirection();
  constructDenselySampledAxis();
  computeAxisAnglesAndDistances();
}

void SurfacePatch::deformAxis(int index, std::complex<double> newDir) {
  assert(index < m_patchAxisSparseAngles.size());
  m_patchAxisSparseAngles[index] = newDir;
  traceAxis();
}

void SurfacePatch::get(std::vector<SurfacePoint>& axis, std::vector<SurfacePoint>& patch) {
  axis = m_patchAxisSparse;
  patch = m_patchPoints;
}

void SurfacePatch::getParameterized(std::vector<AxisPointParams>& parameterizedAxis,
                                    std::vector<PatchPointParams>& parameterizedPatch) {
  parameterizedPatch = m_parameterizedPatchPoints;

  size_t N = m_patchAxisSparseAngles.size();

  std::vector<AxisPointParams> pAxis(N);

  // First axis point direction and all last point parameters will be ignored
  for (int i = 0; i < N; i++) {
    AxisPointParams intermediateParams;

    intermediateParams.nextPointDistance = m_patchAxisSparseDistances[i];
    intermediateParams.nextPointDirection = m_patchAxisSparseAngles[i];

    pAxis[i] = intermediateParams;
  }

  parameterizedAxis = pAxis;
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

std::complex<double> SurfacePatch::getDirAtAxisIndex(int index) {
  assert(index < m_patchAxisSparseAngles.size());
  return m_patchAxisSparseAngles[index];
}

std::vector<std::string> SurfacePatch::getPatchSerialized() { return getSerializedSurfacePoints(m_patchPoints); }

Vector2 SurfacePatch::getInitDir() { return m_initDir; }

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

void SurfacePatch::parameterizePatch() {
  m_parameterizedPatchPoints.clear();

  std::cout << m_patchPoints.size() << std::endl;

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

  std::vector<PatchPointParams> ptToParam(m_patchPoints.size());

  double distance;

  for (size_t i = 0; i < m_patchPoints.size(); i++) {
    SurfacePoint patchPoint = m_patchPoints[i];
    double diffusedVal = patchPoint.interpolate(closestPoint);
    size_t cp = indexOfClosestPointOnAxis(diffusedVal, zippedDistances);
    SurfacePoint axisPoint = m_patchAxisSparse[cp];

    std::vector<SurfacePoint> geodesic = connectPointsWithGeodesicMMP(axisPoint, patchPoint, distance);
    std::vector<SurfacePoint> prunedGeodesic = pruneApproxEqualEntries(geodesic);

    // Compute direction.
    std::complex<double> dir = localDir(axisPoint, prunedGeodesic[1]);

    // Edge case: if distance is 0, then point is on axis
    // In which case just assign an arbitrary direction and let the distance
    // 0 out to prevent divide by 0 complaints
    if (distance <= 0) {
      std::cout << "WARNING (DON'T PANIC): point found to lie on axis" << std::endl;
      dir = std::complex<double>{1.0, 0.0};
    }

    // Compute relative to the tangent direction of the axis
    dir /= std::abs(dir);
    std::complex<double> axisBasis = axisTangent(m_patchAxisSparseDenseIdx[cp], m_patchAxisDense);
    axisBasis /= std::abs(axisBasis);
    dir = dir / axisBasis;

    PatchPointParams prms = {cp, distance, dir};
    ptToParam[i] = prms;
  }

  m_parameterizedPatchPoints.insert(m_parameterizedPatchPoints.begin(), ptToParam.begin(), ptToParam.end());
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

void SurfacePatch::reconstructPatch() {
  m_patchPoints.clear();

  // Use distances/directions to reconstruct patches from sparse axis
  SurfacePoint pathEndpoint;
  TraceGeodesicResult tracedGeodesic;

  std::vector<SurfacePoint> constructedPatch(m_parameterizedPatchPoints.size());

  for (size_t i = 0; i < m_parameterizedPatchPoints.size(); i++) {
    PatchPointParams p = m_parameterizedPatchPoints[i];
    std::complex<double> dir = p.axisPointDirection;
    double distance = p.axisPointDistance;
    SurfacePoint startPoint = m_patchAxisSparse[p.closestAxisPoint];

    std::complex<double> axisBasis = axisTangent(m_patchAxisSparseDenseIdx[p.closestAxisPoint], m_patchAxisDense);
    dir /= std::abs(dir);
    axisBasis /= std::abs(axisBasis);
    dir *= axisBasis;

    // Handle distance = 0 edge cases
    if (distance <= 0) {
      pathEndpoint = startPoint;
    } else {
      tracedGeodesic =
          traceGeodesic(*(m_geometry), m_patchAxisSparse[p.closestAxisPoint], Vector2::fromComplex(distance * dir));
      pathEndpoint = tracedGeodesic.endPoint;
    }

    constructedPatch[i] = pathEndpoint;
  }

  m_patchPoints.insert(m_patchPoints.begin(), constructedPatch.begin(), constructedPatch.end());
}

void SurfacePatch::rotateAxis(Vector2 newDir) {
  m_initDir = newDir;
  traceAxis();
  propagateChildUpdates();
}

void SurfacePatch::saveAxisStartPointAndDirection() {
  m_savedStartPoint = m_startPoint;
  m_savedInitDir = m_initDir;
}

void SurfacePatch::setBulkTransferParams(SurfacePatch* sourcePatch, std::string sourcePatchName,
                                         std::string destinationPatchName) {
  std::cout << "Overriding child trace PatchPointParams of " << destinationPatchName << " with " << sourcePatchName
            << std::endl;

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

void SurfacePatch::transfer(SurfacePatch* target, const SurfacePoint& targetMeshStart,
                            const SurfacePoint& targetMeshDirEndpoit) {
  // All angles and distances should be the same, save for the first one (which is ignored)

  // Note that we deliberately adjust the phasor multiplier since we want the MIRROR image
  // to be generated on the target domain

  target->m_patchAxisSparseAngles.clear();
  target->m_patchAxisSparseDistances.clear();

  target->m_patchAxisSparseDistances.insert(target->m_patchAxisSparseDistances.end(),
                                            m_patchAxisSparseDistances.begin(), m_patchAxisSparseDistances.end());

  target->m_startPoint = targetMeshStart;
  target->m_initDir = target->localDir(targetMeshStart, targetMeshDirEndpoit);

  for (int i = 0; i < m_patchAxisSparseAngles.size(); i++) {
    std::complex<double> adjustedDir = {m_patchAxisSparseAngles[i].real(), -m_patchAxisSparseAngles[i].imag()};
    target->m_patchAxisSparseAngles.push_back(adjustedDir);
  }

  target->traceAxis();

  // Compute distances and directions on S1, then reconstruct contact on S2
  // Reconstruct patch only if patch has been parameterized on source domain

  if (m_parameterizedPatchPoints.size() > 0) {
    target->m_parameterizedPatchPoints.clear();

    for (int i = 0; i < m_parameterizedPatchPoints.size(); i++) {
      PatchPointParams sourceParams = m_parameterizedPatchPoints[i];

      std::complex<double> sourceDirection = sourceParams.axisPointDirection;

      std::complex<double> adjustedDir = {sourceDirection.real(), -sourceDirection.imag()};
      PatchPointParams targetParams = {sourceParams.closestAxisPoint, sourceParams.axisPointDistance, adjustedDir};
      target->m_parameterizedPatchPoints.push_back(targetParams);
    }

    target->reconstructPatch();
  }
}

void SurfacePatch::transferAxisOnly(SurfacePatch* target, const SurfacePoint& targetMeshStart,
                                    const SurfacePoint& targetMeshDirEndpoint) {
  // All angles and distances should be the same, save for the first one (which is ignored)
  target->m_patchAxisSparseAngles = m_patchAxisSparseAngles;
  target->m_patchAxisSparseDistances.insert(target->m_patchAxisSparseDistances.end(),
                                            m_patchAxisSparseDistances.begin(), m_patchAxisSparseDistances.end());

  target->m_startPoint = targetMeshStart;
  target->m_initDir = target->localDir(targetMeshStart, targetMeshDirEndpoint);
  target->traceAxis();
}

void SurfacePatch::transferContactPointsOnly(SurfacePatch* target) {
  // Compute distances and directions on S1, then reconstruct contact on S2
  // Recreate boundary only if boundary has been parameterized on source domain
  // Invert if correspondence should result in mirror image

  if (m_parameterizedPatchPoints.size() > 0) {
    target->m_parameterizedPatchPoints.clear();

    target->m_parameterizedPatchPoints.insert(target->m_parameterizedPatchPoints.begin(),
                                              m_parameterizedPatchPoints.begin(), m_parameterizedPatchPoints.end());

    target->reconstructPatch();
  }
}

void SurfacePatch::translate(const SurfacePoint& newStartPoint) {
  assert(m_savedStartPoint != SurfacePoint());
  assert(m_savedInitDir != Vector2());

  m_initDir = parallelTransport(m_savedStartPoint, newStartPoint, m_savedInitDir.arg());
  m_startPoint = newStartPoint;

  traceAxis();
  propagateChildUpdates();
}

void SurfacePatch::setAxisPoints(std::vector<SurfacePoint>& axisPoints) {
  m_patchAxisSparse.clear();

  m_patchAxisSparse.insert(m_patchAxisSparse.begin(), axisPoints.begin(), axisPoints.end());

  m_startPoint = m_patchAxisSparse[0];

  computeInitialAxisDirection();
  constructDenselySampledAxis();
  computeAxisAnglesAndDistances();
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

void SurfacePatch::setPatchPoints(const std::vector<SurfacePoint>& points) { m_patchPoints = points; }

void SurfacePatch::setParameterizedAxis(std::vector<AxisPointParams>& parameterizedAxisPoints) {
  assert(parameterizedAxisPoints.size() == m_patchAxisSparse.size());
  assert(m_patchAxisSparse.size() == m_patchAxisDense.size());

  m_patchAxisSparseDistances.clear();
  m_patchAxisSparseAngles.clear();

  for (int i = 0; i < parameterizedAxisPoints.size(); i++) {
    AxisPointParams app = parameterizedAxisPoints[i];

    m_patchAxisSparseDistances.push_back(app.nextPointDistance);
    m_patchAxisSparseAngles.push_back(app.nextPointDirection);
  }
}

void SurfacePatch::setParameterizedPatch(std::vector<PatchPointParams>& parameterizedPatchPoints) {
  m_parameterizedPatchPoints.clear();
  m_parameterizedPatchPoints.insert(m_parameterizedPatchPoints.begin(), parameterizedPatchPoints.begin(),
                                    parameterizedPatchPoints.end());
}

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

std::vector<SurfacePoint> SurfacePatch::connectPointsWithGeodesicMMP(const SurfacePoint& pt1, const SurfacePoint& pt2,
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
  double dist;

  // NOTE: Tossing out dense axis - will be same as sparse
  // Workaround until cutting to SurfaceCurve

  for (size_t i = 0; i < N; i++) {
    pt1 = m_patchAxisSparse[i];
    pt2 = m_patchAxisSparse[mod(i + 1, N)];
    denseCurve.push_back(pt1);
    idxIntoDense[i] = currIdx;
    currIdx++;
  }

  m_patchAxisDense = denseCurve;
  m_patchAxisSparseDenseIdx = idxIntoDense;
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

  return -1;
}

Vector2 SurfacePatch::inTangentBasis(const SurfacePoint& pA, const SurfacePoint& pB, const SurfacePoint& p) {

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

Vector2 SurfacePatch::parallelTransport(const SurfacePoint& startPoint, const SurfacePoint& endpoint,
                                        double initAngle) {
  double distance;

  std::vector<SurfacePoint> path = connectPointsWithGeodesicMMP(startPoint, endpoint, distance);

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
