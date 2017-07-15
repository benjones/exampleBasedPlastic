#pragma once

#include "eigenTypedefs.h"
#include <vector>
#include "Quat.h"

#include "cl.hpp"

#include <LinearMath/btTransform.h>
#include <LinearMath/btDefaultMotionState.h>
#include <BulletCollision/Gimpact/btGImpactShape.h>
#include <BulletDynamics/Dynamics/btRigidBody.h>
#include <BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h>
#include <BulletCollision/CollisionShapes/btCompoundShape.h>
#include <BulletCollision/CollisionShapes/btTetrahedronShape.h>
#include <BulletCollision/CollisionShapes/btSphereShape.h>

using TetrahedronShape = btBU_Simplex1to4;


class ExampleGraph;
class PlasticBody;
class World;

class PlasticPiece{
public:

  //move only type because of the unique_ptrs
  //also, we should never actually need to copy them

  
  PlasticPiece(){}
  ~PlasticPiece() = default;
  PlasticPiece(PlasticPiece&&) = default;
  PlasticPiece& operator=(PlasticPiece&&) = default;
  PlasticPiece(const PlasticPiece&) = delete;
  PlasticPiece& operator=(const PlasticPiece&) = delete;
  
  void initialize(const std::string& directory,
	  const PlasticBody& parent,
	  const RMMatrix3d& verticesIn,
	  World& world,
	  int pieceNumber);

  void initializeNoDynamics(const std::string& directory, 
	  const PlasticBody& parent,
	  const RMMatrix3d& verticesIn);

  void computeMassesAndVolume(double density);

  void updateBulletProperties();

  int getNearestVertex(const Eigen::Vector3d& localPoint) const;
  size_t getNearestSphere(const btVector3& localPoint) const;

  
  void saveBulletSnapshot();
  void restoreBulletSnapshot();

  void skinMeshVaryingBarycentricCoords(const Eigen::MatrixXd& boneWeights,
										const std::vector<int>& boneIndices,
										const ExampleGraph& exampleGraph);
  
  void skinMeshOpenCL(World& world, PlasticBody& parent, cl::Kernel& clKernel);

  void skinSpheres(World& world);

  
  RMMatrix3i computeAllTriangles(const std::vector<size_t>& tetIndices) const;
  RMMatrix3i computeTriangleFaces(const std::vector<size_t>& tetIndices) const;

  //neighbors within this piece.  
  //We can just union with vertices
  //in other pieces starting at vertices with the same index
  void computeVertexNeighbors();

  //propogate out from the seed points
  //there can be more than one if the real seed is in a different piece
  void geodesicDistancesPropogate(const std::vector<size_t>& seeds,
	double radius);
  

  struct BulletSnapshot{
	btVector3 linearVelocity, angularMomentum;
	btTransform comTransform;
  };
  BulletSnapshot bulletSnapshot;
  
  
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>  barycentricCoordinates; //per vertex
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> deltaBarycentricCoordinates;


  //std::vector<Eigen::Vector3d> perVertexTranslations; 
  //std::vector<Quat> perVertexRotations;
  //should call world.removeRigidBody....
  std::unique_ptr<btRigidBody> bulletBody;

  std::unique_ptr<btCompoundShape> bulletShape;
  std::vector<TetrahedronShape> bulletTets;
  std::unique_ptr<btDefaultMotionState> motionState;

  //for the planar fracture mesh
  std::unique_ptr<btTriangleIndexVertexArray> btTriMesh; //only used for trimesh shapes
  std::unique_ptr<btGImpactMeshShape> fractureMesh;

  //as an alternative to the tetmesh approach
  bool useVolumetricCollisions;
  std::unique_ptr<btTriangleIndexVertexArray> bulkTriMesh;
  std::unique_ptr<btGImpactMeshShape> bulkMesh;

  size_t numSpheres;
  std::vector<double> sphereMasses;
  std::vector<btSphereShape> bulletSpheres;
  //which tetmesh vertex is closest to the sphere center?
  std::vector<size_t> sphereToVertexMap;
  std::vector<Eigen::Vector3d> originalSpherePositions, skinnedSpherePositions;
  std::vector<double> skinnedSphereRadii;



  
  std::vector<std::vector<Eigen::Vector3d> > deformedSpherePositions;
  std::vector<std::vector<double> > deformedSphereRadii;
  
  void computeSphereToVertexMap();
  void loadDeformedSpheres(std::ifstream& ins);
  
  RMMatrix3d tetmeshVertices;
  RMMatrix4i tetmeshTets;
  RMMatrix3i tetmeshTriangles; 
  Eigen::VectorXd tetmeshVertexMasses;  //for computing inertia tensor, etc

  int numPhysicsVertices; //include the ones that aren't actually part of any
  //tets in this piece



  RMMatrix3d currentBulletVertexPositions;
  RMMatrix3f localOffsets;



  
  double mass, volume, scaleFactor;

  Eigen::Matrix3d inertiaTensor;
  btTransform inertiaAligningTransform, worldTransform;

  std::vector<std::vector<size_t>> vertexNeighbors;
  std::vector<std::vector<double>> neighborDistances;
  std::vector<double> geodesicDistances; //to a given point

  void computeActiveVertices();
  std::vector<size_t> activeVertices; //vertices that are actually part of tets
  void computeCollisionTetIndices();
  std::vector<size_t> cutTets; //read from slicer output file
  std::vector<size_t> collisionTetIndices; 

  //indices of the vertices in the containing tet,
  //and the corresponding barycentric coordinates
  void computeClippingPlaneTetInfo();
  std::vector<std::tuple<Eigen::Vector4i,Eigen::Vector4d>> clippingPlaneTetInfo;
  


  RMMatrix3d splittingPlaneVertices;
  RMMatrix3d currentBulletSplittingPlaneVertices;
  RMMatrix3i splittingPlaneTriangles;


  
  void dumpPly(const std::string& filename) const;
  //void dumpImpulses(const std::string& filename) const; //skip for now
  void dumpBcc(const std::string& filename) const;
  void dumpSpheres(const std::string& filename) const;
  
  void dumpColoredShapeOnly(const std::string& filename) const;

  void updateAabbs();

  //maybe should go in plasticBody, but whatever
  struct Plane{
	Plane(Eigen::Vector3d n, double o)
	  :normal{std::move(n)}, offset{o}
	{}

	Eigen::Vector3d normal;
	double offset;
  };
  std::vector<Plane> planes;
  Eigen::Vector3d intersectPlane(size_t ind1, size_t ind2, const Plane& plane) const;

  void computeCutTetGeometry();
  RMMatrix3d cutGeomVertices;
  RMMatrix3i cutGeomIndices;

  //Eigen::VectorXf texU, texV;

  //CL stuff

  std::vector<float> hostBarycentricCoordinates;
  std::vector<float> hostSkinnedPositions;
  //  std::vector<float> hostPerVertexTranslations;
  //  std::vector<float> hostPerVertexRotations;

  cl::Buffer deviceBarycentricCoordinates, deviceSkinnedPositions;
  //devicePerVertexTranslations, devicePerVertexRotations;

  int framesToSkin;

  bool plinkoObject;

  void computeAndPrintExtents() const;
  
};
