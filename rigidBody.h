#pragma once

//#include <BulletCollision/CollisionShapes/btCollisionShape.h>
//#include <BulletDynamics/Dynamics/btRigidBody.h>
#include <memory> //for unique_ptr
//#include "couplingConstraint.h"
//#include "couplingConstraintSolver.h"
#include <LinearMath/btVector3.h>
#include <vector>
#include "Eigen/Dense"
#include <BulletCollision/CollisionShapes/btStaticPlaneShape.h>
#include <BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h>

class btCollisionShape;

class btRigidBody;
class btMotionState;


class RigidBody{
  public:
  RigidBody() {} //empty constructor
  
  void dump(std::string filename);

  //no copying or moving allowed
  RigidBody(const RigidBody& rb) = delete;
  RigidBody(RigidBody&& rb) = default; //move only, no copying

  std::unique_ptr<btRigidBody> bulletBody;
  std::unique_ptr<btMotionState> motionState;
  std::unique_ptr<btCollisionShape> shape;

  
  std::unique_ptr<btTriangleIndexVertexArray> triMesh; //only used for trimesh shapes
  void loadTrimesh(std::string filename, double scale);
  std::vector<double> meshVertices;
  std::vector<int> meshIndices;


  //CouplingConstraintSolver couplingSolver; //there will be one per RB, keep it here to avoid allocations
  //std::vector<CouplingConstraint> constraints;

  enum RBType{
    RB_BOX,
    RB_PLANE,
	RB_TRIMESH,
	RB_SPHERE
  };
  
  RBType rbType;
  int constantForceFrames;
  btVector3 constantForce; //apply this force each frame
  
};


inline Eigen::Matrix<float, 4, 3, Eigen::RowMajor>
getPlaneVertices(const btStaticPlaneShape* planeShape){
	auto normal = planeShape->getPlaneNormal();
	auto localPoint = -planeShape->getPlaneConstant()*normal;
	
	auto span1 = normal.cross(btVector3{1, 0, 0});
	if(span1.length2() < 1e-3){ //normal is nearly parallel to 1, 0, 0, so don't use it
	  span1 = normal.cross(btVector3{0, 0, 1});
	}
	auto span2 = span1.cross(normal);

	Eigen::Matrix<float, 4, 3, Eigen::RowMajor> vertices;
	auto tmp = localPoint - 100*span1 - 100*span2;
	//printPoint(tmp);
	vertices(0,0) = tmp[0];
	vertices(0,1) = tmp[1];
	vertices(0,2) = tmp[2];
	tmp = localPoint + 100*span1 - 100*span2;
	//printPoint(tmp);
	vertices(1,0) = tmp[0];
	vertices(1,1) = tmp[1];
	vertices(1,2) = tmp[2];
	tmp = localPoint + 100*span1 + 100*span2;
	//printPoint(tmp);
	vertices(2,0) = tmp[0];
	vertices(2,1) = tmp[1];
	vertices(2,2) = tmp[2];
	tmp = localPoint - 100*span1 + 100*span2;
	//printPoint(tmp);
	vertices(3,0) = tmp[0];
	vertices(3,1) = tmp[1];
	vertices(3,2) = tmp[2];


	return vertices;
}
