#pragma once

//#include <BulletCollision/CollisionShapes/btCollisionShape.h>
//#include <BulletDynamics/Dynamics/btRigidBody.h>
#include <memory> //for unique_ptr
#include "couplingConstraint.h"
#include "couplingConstraintSolver.h"
#include <vector>

class btCollisionShape;
class btTriangleIndexVertexArray;
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


  CouplingConstraintSolver couplingSolver; //there will be one per RB, keep it here to avoid allocations
  std::vector<CouplingConstraint> constraints;

  enum RBType{
    RB_BOX,
    RB_PLANE,
	RB_TRIMESH
  };
  
  RBType rbType;

  btVector3 constantForce; //apply this force each frame
};
