#pragma once

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/exact_geodesics.h"
#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/vector2.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

struct AxisPointParams {
  double nextPointDistance;                // distance to next point
  std::complex<double> nextPointDirection; // outgoing direction relative to incoming direction
};

struct PatchPointParams {
  size_t closestAxisPoint;                 // index of closest point on axis
  double axisPointDistance;                // distance away from axis point
  std::complex<double> axisPointDirection; // outgoing direction from axis point
};

class SurfacePatch {

public:
  // constructors
  SurfacePatch(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry, GeodesicAlgorithmExact* mmpSolver,
               VectorHeatMethodSolver* vectorHeatSolver);

  void computeInitialAxisDirection();

  void createCustomAxis(std::vector<Vertex>& axisPoints);

  void createDefaultAxis();

  void deformAxis(int index, std::complex<double> newDir);

  void get(std::vector<SurfacePoint>& axis, std::vector<SurfacePoint>& patch);

  void getParameterized(std::vector<AxisPointParams>& parameterizedAxis,
                        std::vector<PatchPointParams>& parameterizedPatch);

  std::complex<double> getDirAtAxisIndex(int index);

  Vector2 getInitDir();

  double getPatchSpreadCoefficient();

  void parameterizePatch();

  void reconstructPatch();

  void rotateAxis(Vector2 newDir);

  void saveAxisStartPointAndDirection();

  void setPatchSpreadCoefficient(double patchSpreadCoefficient);

  void transfer(SurfacePatch* target, const SurfacePoint& targetMeshStart, const SurfacePoint& targetMeshDirEndpoint);

  void translate(const SurfacePoint& newStartPoint);

  void setAxisPoints(std::vector<SurfacePoint>& axisPoints);

  void setPatchAxis(const std::vector<SurfacePoint>& axis, const std::vector<std::complex<double>>& dirs,
                    const std::vector<double>& dists);

  void setPatchPoints(const std::vector<SurfacePoint>& points);

  void setParameterizedAxis(std::vector<AxisPointParams>& parameterizedAxisPoints);

  void setParameterizedPatch(std::vector<PatchPointParams>& parameterizedPatchPoints);

private:
  bool approxEqual(const SurfacePoint& pA, const SurfacePoint& pB);

  std::complex<double> axisTangent(size_t idx, const std::vector<SurfacePoint>& axis);

  void computeAxisAnglesAndDistances();

  std::vector<SurfacePoint> connectPointsWithGeodesic(const SurfacePoint& pt1, const SurfacePoint& pt2,
                                                      double& distance);

  void constructDenselySampledAxis();

  size_t indexOfClosestPointOnAxis(double diffusedVal,
                                   const std::vector<std::tuple<SurfacePoint, double>>& zippedDistances);

  Vector2 inTangentBasis(const SurfacePoint& pA, const SurfacePoint& pB, const SurfacePoint& p);

  Vector2 localDir(const SurfacePoint& pt1, const SurfacePoint& pt2);

  Vector2 parallelTransport(const SurfacePoint& startPoint, const SurfacePoint& endpoint, double initAngle);

  std::vector<SurfacePoint> pruneApproxEqualEntries(const std::vector<SurfacePoint>& source);

  void traceAxis();

  std::unique_ptr<ManifoldSurfaceMesh> m_mesh;
  std::unique_ptr<VertexPositionGeometry> m_geometry;
  std::unique_ptr<VectorHeatMethodSolver> m_vectorHeatSolver;
  std::unique_ptr<GeodesicAlgorithmExact> m_mmpSolver;

  std::vector<SurfacePoint> m_patchAxisSparse;
  std::vector<std::complex<double>> m_patchAxisSparseAngles;
  std::vector<double> m_patchAxisSparseDistances;
  std::vector<size_t> m_patchAxisSparseDenseIdx;

  std::vector<SurfacePoint> m_patchAxisDense;

  std::vector<PatchPointParams> m_parameterizedPatchPoints;
  std::vector<SurfacePoint> m_patchPoints;

  SurfacePoint m_startPoint;
  Vector2 m_initDir;

  SurfacePoint m_savedStartPoint;
  Vector2 m_savedInitDir;

  double m_patchSpreadCoefficient;
};
