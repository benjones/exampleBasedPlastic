#pragma once

#include <BulletCollision/CollisionShapes/btCollisionShape.h>
#include <BulletDynamics/Dynamics/btRigidBody.h>
#include <memory> //for unique_ptr


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

  enum RBType{
    RB_BOX
  };
  
  RBType rbType;
};
