#pragma once
#include <vector>
#include <list>
#include <BulletCollision/BroadphaseCollision/btBroadphaseInterface.h>
#include <BulletCollision/CollisionDispatch/btDefaultCollisionConfiguration.h>
#include <BulletCollision/BroadphaseCollision/btDispatcher.h>
#include <BulletDynamics/ConstraintSolver/btSequentialImpulseConstraintSolver.h>
#include <BulletDynamics/Dynamics/btDiscreteDynamicsWorld.h>

#include "rigidBody.h"
#include "plasticBody.h"

#include "profiler.hpp"
#include <tbb/tbb.h>

namespace Json{ class Value;}

class World{
  public:

  //  World() {} 
  
  World(std::string filename);
  World(const World& other) = delete;  //no copying
  World(World&& other) = delete;//moving only

  World& operator=(World&& other) = delete;
  World& operator=(const World& other) = delete;


  //void loadPlasticObjects(const Json::Value& root);

  void loadPlasticBodies(const Json::Value& root);

  void timeStepDynamicSprites();
  // void timeStepDynamicSpritesNoDouble(); //hold off on this for now

  
  void collectImpulses();


  void dumpFrame();



  const int pyramidSize = 4;
  int getNumBarrels();
  void makeBarrelPyramid();


  std::unique_ptr<btBroadphaseInterface> broadphaseInterface;  
  std::unique_ptr<btCollisionConfiguration> collisionConfiguration;
  std::unique_ptr<btCollisionDispatcher> dispatcher;

  //maybe using the MLCP solver makes sense, since it's supposed to give higher fidelity results
  std::unique_ptr<btSequentialImpulseConstraintSolver> bulletSolver;

  std::unique_ptr<btDiscreteDynamicsWorld> bulletWorld;

  bool ground;
  double groundHeight;
  double dt, duration;
  btVector3 gravity;
  double friction;
  
  bool makeBarrels;
  bool singleStep;

  std::vector<RigidBody> rigidBodies;
  std::vector<PlasticBody> plasticBodies;


  int currentFrame;
 
  benlib::Profiler profiler;


  //openCL stuff

  void initCL();

  cl::Device device;
  cl::Context context;
  cl::CommandQueue queue;

  cl::Program clProgram;
  
  tbb::enumerable_thread_specific<cl::Kernel> clKernels;


};
