#pragma once
#include <vector>
#include <list>
#include <BulletCollision/BroadphaseCollision/btBroadphaseInterface.h>
#include <BulletCollision/CollisionDispatch/btDefaultCollisionConfiguration.h>
#include <BulletCollision/BroadphaseCollision/btDispatcher.h>
#include <BulletDynamics/ConstraintSolver/btSequentialImpulseConstraintSolver.h>
#include <BulletDynamics/Dynamics/btDiscreteDynamicsWorld.h>

#include "rigidBody.h"
#include "fem.H"



class World{
  public:
  
  World(std::string filename);
  World(const World& other) = delete;  //move only, no copying
  World(World&& other) = default;


  void timeStep();
  void dumpFrame();

  void computeFemVelocities();
  void updateFemPositions();


  std::unique_ptr<btBroadphaseInterface> broadphaseInterface;  
  std::unique_ptr<btCollisionConfiguration> collisionConfiguration;
  std::unique_ptr<btDispatcher> dispatcher;

  //maybe using the MLCP solver makes sense, since it's supposed to give higher fidelity results
  std::unique_ptr<btSequentialImpulseConstraintSolver> bulletSolver;

  btDiscreteDynamicsWorld bulletWorld;

  bool ground;
  double groundHeight;
  double dt;
  btVector3 gravity;
  double friction;

  std::vector<RigidBody> rigidBodies;
  std::vector<FemObject> femObjects;

  int currentFrame;

};
