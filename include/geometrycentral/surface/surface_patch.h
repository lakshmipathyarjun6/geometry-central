#pragma once

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/exact_geodesics.h"
#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/heat_method_distance.h"
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
               HeatMethodDistanceSolver* distanceHeatSolver, VectorHeatMethodSolver* vectorHeatSolver);

  void computeGlobalAxisTransport();

  void computeInitialAxisDirection();

  void createCustomAxis(std::vector<Vertex>& axisPoints);

  void createDefaultAxis();

  void deformAxis(int index, std::complex<double> newDir);

  void get(std::vector<SurfacePoint>& axis, std::vector<SurfacePoint>& patch);

  void getParameterized(std::vector<AxisPointParams>& parameterizedAxis,
                        std::vector<PatchPointParams>& parameterizedPatch);

  std::vector<std::string> getAxisSerialized();

  std::complex<double> getDirAtAxisIndex(int index);

  std::vector<std::string> getPatchSerialized();

  Vector2 getInitDir();

  void linkPatch(std::string childName, SurfacePatch* child);

  void parameterizePatch();

  void propagateChildUpdates();

  void reconstructPatch();

  void rotateAxis(Vector2 newDir);

  void saveAxisStartPointAndDirection();

  void setBulkTransferParams(SurfacePatch* sourcePatch, std::string sourcePatchName, std::string destinationPatchName);

  void transfer(SurfacePatch* target, const SurfacePoint& targetMeshStart, const SurfacePoint& targetMeshDirEndpoint);

  void transferAxisOnly(SurfacePatch* target, const SurfacePoint& targetMeshStart,
                        const SurfacePoint& targetMeshDirEndpoint);

  void transferContactPointsOnly(SurfacePatch* target);

  void translate(const SurfacePoint& newStartPoint);

  void setAxisPoints(std::vector<SurfacePoint>& axisPoints);

  void setPatchAxis(const std::vector<SurfacePoint>& axis, const std::vector<std::complex<double>>& dirs,
                    const std::vector<double>& dists);

  void setPatchPoints(const std::vector<SurfacePoint>& points);

  void setParameterizedAxis(std::vector<AxisPointParams>& parameterizedAxisPoints);

  void setParameterizedPatch(std::vector<PatchPointParams>& parameterizedPatchPoints);

  void unlinkAllPatches();

  void unlinkPatch(std::string childName);

private:
  bool approxEqual(const SurfacePoint& pA, const SurfacePoint& pB);

  std::complex<double> axisTangent(size_t idx, const std::vector<SurfacePoint>& axis);

  void computeAxisAnglesAndDistances();

  std::vector<SurfacePoint> connectPointsWithGeodesic(const SurfacePoint& pt1, const SurfacePoint& pt2);

  std::vector<SurfacePoint> connectPointsWithGeodesicMMP(const SurfacePoint& pt1, const SurfacePoint& pt2,
                                                         double& distance);

  void constructDenselySampledAxis();

  double evaluateVertexDataAtPoint(const VertexData<double>& u, const SurfacePoint& pt);

  std::vector<std::string> getSerializedSurfacePoints(std::vector<SurfacePoint>& source);

  size_t indexOfClosestPointOnAxis(double diffusedVal,
                                   const std::vector<std::tuple<SurfacePoint, double>>& zippedDistances);

  Vector2 inTangentBasis(const SurfacePoint& pA, const SurfacePoint& pB, const SurfacePoint& p);

  Vector2 localDir(const SurfacePoint& pt1, const SurfacePoint& pt2);

  Vector2 parallelTransport(const SurfacePoint& startPoint, const SurfacePoint& endpoint, double initAngle);

  std::vector<SurfacePoint> pruneApproxEqualEntries(const std::vector<SurfacePoint>& source);

  void traceAxis();

  std::unique_ptr<ManifoldSurfaceMesh> m_mesh;
  std::unique_ptr<VertexPositionGeometry> m_geometry;
  std::unique_ptr<HeatMethodDistanceSolver> m_distanceHeatSolver;
  std::unique_ptr<VectorHeatMethodSolver> m_vectorHeatSolver;
  std::unique_ptr<GeodesicAlgorithmExact> m_mmpSolver;

  std::vector<SurfacePoint> m_patchAxisSparse;
  std::vector<std::complex<double>> m_patchAxisSparseAngles;
  std::vector<double> m_patchAxisSparseDistances;
  std::vector<size_t> m_patchAxisSparseDenseIdx;

  VertexData<Vector2> m_axisDirectionTable;
  std::vector<SurfacePoint> m_patchAxisDense;

  std::vector<PatchPointParams> m_parameterizedPatchPoints;
  std::vector<SurfacePoint> m_patchPoints;

  // To support hierarchal organization
  SurfacePatch* m_parent;
  std::map<std::string, SurfacePatch*> m_children;
  std::map<std::string, std::tuple<std::complex<double>, std::complex<double>, double, double>> m_childTraceParams;

  SurfacePoint m_startPoint;
  Vector2 m_initDir;

  SurfacePoint m_savedStartPoint;
  Vector2 m_savedInitDir;
};
