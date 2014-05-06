#pragma once

#include <memory>
#include <exampleGraph.h>
#include "egTraverser.h"
#include <Skeleton.cpp> //templates...
#include <Bone.h>

#include <LinearMath/btDefaultMotionState.h>
#include <BulletCollision/Gimpact/btGImpactShape.h>
#include <BulletCollision/CollisionShapes/btCompoundShape.h>
#include <BulletDynamics/Dynamics/btRigidBody.h>
#include <BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h>


using RMMatrix3d = Eigen::Matrix<double,Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix3i = Eigen::Matrix<int,   Eigen::Dynamic, 3, Eigen::RowMajor>;
using RMMatrix4i = Eigen::Matrix<int,   Eigen::Dynamic, 4, Eigen::RowMajor>;

class PlasticObject{
public:
  PlasticObject(){}
  ~PlasticObject() = default;
  PlasticObject(PlasticObject&&) = default;
  PlasticObject(const PlasticObject&) = delete;

  void writeOutputFiles(std::string filename);
  
  void loadFromFiles(std::string directory);
  
  void computeMassesAndVolume();
  
  void updateBulletProperties(const RMMatrix3d& vertices,
							  const RMMatrix4i& tets);
  

  void deformBasedOnImpulses(btPersistentManifold* man, bool isObject0);

  Eigen::Vector3d getBarycentricCoordinates(btVector3 localPoint, int triangleIndex);
  int getNearestVertex(Eigen::Vector3d localPoint);
  
  
  void skinMesh();

  void dump(std::string filename);

  ExampleGraph exampleGraph;
  EGTraverser egTraverser;
  std::unique_ptr<btRigidBody> bulletBody;
  std::unique_ptr<btGImpactMeshShape> bulletShape;
  std::unique_ptr<btDefaultMotionState> motionState;
  
  
  std::unique_ptr<btTriangleIndexVertexArray> btTriMesh; //only used for trimesh shapes
  //std::vector<double> trimeshVertices;
  //std::vector<int> trimeshIndices;
  
  RMMatrix3d trimeshVertices, tetmeshVertices;
  RMMatrix3i trimeshTriangles, tetmeshTriangles;
  RMMatrix4i tetmeshTets;
  
  RMMatrix3d currentBulletVertexPositions;


  double density;
  double mass;
  double volume;
  Eigen::Matrix3d inertiaTensor;
  
  Eigen::VectorXd tetmeshVertexMasses;
  
  Skeleton<Bone> skeleton;
  std::vector<Bone*> boneRoots;
  
  Eigen::MatrixXd boneWeights;
  Eigen::MatrixXd currentTransformMatrix;
  
  //an extra term to add to the translation of each handle based on
  //collision impacts
  Eigen::MatrixXd collisionHandleAdjustments;
  //only deform once things cross this point
  double plasticityImpulseYield;
  //scale*(mag(impulse) - yield) gets distributed to handles
  double plasticityImpulseScale;



  Eigen::MatrixXd LBSMatrix;

  btTransform inertiaAligningTransform; //align the skinned position so that it's 
  //centered at its COM and aligned with the principal axes of inertia
  
  btTransform worldTransform; //the transform bullet gets is the product of these 2

};
