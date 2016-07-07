
#include "world.h"
#include "rigidBody.h"

#include <fstream>
#include <iostream>
#include <BulletCollision/BroadphaseCollision/btDbvtBroadphase.h>
#include <BulletCollision/CollisionShapes/btStaticPlaneShape.h>
#include <LinearMath/btVector3.h>

#include <BulletCollision/CollisionShapes/btBoxShape.h>
#include <BulletCollision/CollisionShapes/btSphereShape.h>
#include <BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h>

#include "json/json.h"

#include "utils.h"
#include <numeric>
#include <limits>

#include <tbb/tbb.h>

#include "range.hpp"
#include "enumerate.hpp"
using benlib::range;
using benlib::enumerate;




World::World(std::string filename)
  :broadphaseInterface{std::unique_ptr<btBroadphaseInterface>{
	new btDbvtBroadphase()}},
   collisionConfiguration{std::unique_ptr<btCollisionConfiguration>{
	   new btDefaultCollisionConfiguration()}},
   dispatcher{std::unique_ptr<btCollisionDispatcher>{
	   new btCollisionDispatcher{collisionConfiguration.get()}}},
   bulletSolver{std::unique_ptr<btSequentialImpulseConstraintSolver>{
	   new btSequentialImpulseConstraintSolver()}},
   bulletWorld{std::unique_ptr<btDiscreteDynamicsWorld>{
	   new btDiscreteDynamicsWorld{
	   dispatcher.get(), 
		 broadphaseInterface.get(), 
		 bulletSolver.get(), 
		 collisionConfiguration.get()}}},
   currentFrame{0}
{

  initCL();



  btGImpactCollisionAlgorithm::registerAlgorithm(dispatcher.get());

  std::ifstream ins(filename.c_str());
  Json::Value root;
  Json::Reader reader;

  bool jsonSuccess = reader.parse(ins, root);
  ins.close();

  if(!jsonSuccess){
    std::cout << "couldn't read input file: " << filename << std::endl 
	      << reader.getFormattedErrorMessages();
    exit(1);
  }

  makeBarrels = root.get("makeBarrels", false).asBool();
  singleStep = root.get("singleStep", true).asBool();
  auto& bulletObjectsIn = root["bulletObjects"];


  ground = root.get("ground", true).asBool();

  if(ground){ //default to 0
    groundHeight = root.get("groundHeight", 0.0).asDouble();
  
    //set up the ground RB
    rigidBodies.emplace_back();
    auto& RB = rigidBodies.back();
	auto normal = btVector3{0,1,0};
	normal.normalize();
    RB.shape = 
      std::unique_ptr<btCollisionShape>(new btStaticPlaneShape(normal, groundHeight));
    RB.rbType = RigidBody::RBType::RB_PLANE;
    RB.motionState = 
      std::unique_ptr<btMotionState>(new btDefaultMotionState(btTransform(btQuaternion::getIdentity(),
									  btVector3(0,0,0))));
    RB.bulletBody = std::unique_ptr<btRigidBody>(new btRigidBody(0, 
								 RB.motionState.get(),
								 RB.shape.get(),
								 btVector3(0,0,0)));
								 
	RB.bulletBody->setUserIndex(-1);
	RB.bulletBody->setRestitution(0.5);
	//RB.bulletBody->setFriction(0.8);
	if(!root["groundFriction"].isNull()){
	  RB.bulletBody->setFriction(root.get("groundFriction", 0.99).asDouble());
	}
	if(!root["groundRestitution"].isNull()){
	  RB.bulletBody->setRestitution(root["groundRestitution"].asDouble());
	}
    bulletWorld->addRigidBody(RB.bulletBody.get());

  }

  dt = root.get("dt", 1.0/60.0).asDouble();
  duration = root.get("duration", 5.0).asDouble();
  friction = root.get("friction", 0.5).asDouble();

  auto gravityIn = root["gravity"];
  if(!gravityIn.isNull() && gravityIn.isArray()){
    gravity = btVector3(gravityIn[0].asDouble(),
			gravityIn[1].asDouble(),
			gravityIn[2].asDouble());
  } else if (!gravityIn.isNull() && gravityIn.isNumeric()){
    gravity = btVector3(0, gravityIn.asDouble(), 0);
  } else {
    gravity = btVector3(0, -9.81, 0);
  }


  bulletWorld->setGravity(gravity);
  for(auto i : range(bulletObjectsIn.size())){
    std::cout << "adding rb: " << i << " "<<bulletObjectsIn[i]["mass"] <<std::endl;
    auto bo = bulletObjectsIn[i];
    rigidBodies.emplace_back();
    auto& RB = rigidBodies.back();

    if(!bo["mass"].isNumeric()){
      std::cout << "rigid body mass must be provided" << std::endl;
      exit(1);
    }
    double mass = bo["mass"].asDouble();

	if(bo["constantForce"].isNull()){
	  RB.constantForce = btVector3{0,0,0};
	} else {
	  RB.constantForce = btVector3{ bo["constantForce"][0].asDouble(),
									bo["constantForce"][1].asDouble(),
									bo["constantForce"][2].asDouble()};
	}
	RB.constantForceFrames = bo.get("constantForceFrames", -1).asInt();
    btVector3 offset(0.0,0.0,0.0); //default to no offset
    auto offsetIn = bo["offset"];
    if(offsetIn.isArray()){
      assert(offsetIn.size() == 3);
      offset = btVector3(offsetIn[0].asDouble(),
			 offsetIn[1].asDouble(),
			 offsetIn[2].asDouble());
    }
    
    auto rotation = btQuaternion::getIdentity();
    auto rotationIn = bo["rotation"];
    if(!rotationIn.isNull() && rotationIn.isArray()){
      assert(rotationIn.size() == 4);
      rotation = btQuaternion(rotationIn[0].asDouble(),
			      rotationIn[1].asDouble(),
			      rotationIn[2].asDouble(),
			      rotationIn[3].asDouble());
    }

    RB.motionState = 
      std::unique_ptr<btMotionState>(new btDefaultMotionState(btTransform(rotation, offset)));



	//LOAD SHAPE SPECIFIC STUFF
    auto shapeTypeIn = bo["shape"];
    if(shapeTypeIn.asString() == "box"){
      RB.rbType = RigidBody::RBType::RB_BOX;
      auto extentsIn = bo["extents"];
      if(!extentsIn.isArray()){
		std::cout << "boxest must have extents, half widths of the box" << std::endl;
		exit(1);
      }
      btVector3 extents{extentsIn[0].asDouble(),
		  extentsIn[1].asDouble(),
		  extentsIn[2].asDouble()};
      
      RB.shape = std::unique_ptr<btCollisionShape>(new btBoxShape(extents));
      
    } else if(shapeTypeIn.asString() == "mesh"){
	  RB.rbType = RigidBody::RBType::RB_TRIMESH;
	  RB.loadTrimesh(bo["filename"].asString(), bo.get("scale", 1.0).asDouble());
	  std::cout<<"loading "<<bo["filename"]<<std::endl;
	}else if (shapeTypeIn.asString() == "plane"){

	  RB.rbType = RigidBody::RBType::RB_PLANE;
	  btVector3 normal;
	  for(auto k : range(3)){
		normal[k] = bo["normal"][k].asDouble();
	  }
	  normal.normalize();
	  double offset = bo.get("offset", 0.0).asDouble();
	  RB.shape = std::unique_ptr<btCollisionShape>(
		  new btStaticPlaneShape(normal, offset));
	} else if(shapeTypeIn.asString() == "sphere"){
	  RB.rbType = RigidBody::RBType::RB_SPHERE;
	  double sphereRadius = bo["radius"].asDouble();
	  RB.shape = std::unique_ptr<btCollisionShape>(
		  new btSphereShape(sphereRadius));
	  
	} else { //add other shape types here
      std::cout << "unknown shape: " << bo["shape"].asString() << std::endl;
      exit(1);
    }



	
    btVector3 inertiaTensor(0,0,0);
	
    RB.shape->calculateLocalInertia(mass, inertiaTensor);
    RB.bulletBody = std::unique_ptr<btRigidBody> {
	  new btRigidBody{mass, 
					  RB.motionState.get(), 
					  RB.shape.get(),
					  inertiaTensor}};
	
	RB.bulletBody->setUserIndex(-1);
	//RB.bulletBody->setFriction(0.8);
    bulletWorld->addRigidBody(RB.bulletBody.get());
    
  }

  loadPlasticBodies(root);
  makeBarrelPyramid();  


}

void World::dumpFrame(){

  char framestring[80];
  sprintf(framestring, "frames/foo-%%03i.%%04i.ply");
  
  std::cout << "writing frame: " << currentFrame << std::endl;

  int objectCount = 0;
  char fname[80];

  for(auto& rigidBody : rigidBodies){

	if(currentFrame == 0 || rigidBody.bulletBody->getInvMass() != 0){
	  sprintf(fname, framestring, objectCount, currentFrame);
	  rigidBody.dump(fname);
	}
	objectCount++;
  }
  
  for(auto& plasticBody : plasticBodies){
	objectCount = plasticBody.dump(currentFrame, objectCount);
  }
	  
  currentFrame++;

  profiler.dumpPercentages(std::cout);
}


void World::loadPlasticBodies(const Json::Value& root){
  
  plasticBodies.clear();
  auto& plasticBodiesIn = root["plasticBodies"];
  plasticBodies.reserve(plasticBodiesIn.size() + getNumBarrels());

  //way easier than before...
  for(auto& pbi : plasticBodiesIn){
	plasticBodies.emplace_back();
	auto& pb = plasticBodies.back();
	pb.loadFromJson(pbi, *this, plasticBodies.size() -1 );
  }
}




/*void World::timeStepDynamicSpritesNoDouble(){


  bulletWorld.stepSimulation(dt, 10, dt);
  for(auto& po : plasticObjects){
	po.saveBulletSnapshot();
  }
  
  collectImpulses();

  //deformBasedOnImpulses();
  for(auto& po : plasticObjects){
	
	//po.projectImpulsesOntoExampleManifold();
	po.projectImpulsesOntoExampleManifoldLocally();

	po.skinMesh();//bulletWorld);
	po.updateBulletProperties(po.currentBulletVertexPositions, 
							  po.tetmeshTets);
	//po.updateCompoundShape();
	//	po.restoreBulletSnapshot();
  }
  //  bulletWorld.stepSimulation(dt, 10, dt);

  }*/

void World::timeStepDynamicSprites(){
  /*
  bulletWorld.stepSimulation(dt, 10, dt);
  for(auto& po: plasticBodies){
	for(auto& pp : po.plasticPieces){
	  //if(!singleStep){
	  //  pp.bulletBody->setRestitution(po.restitution);
	  //}
	  if(po.hasConstantVelocity && currentFrame < po.constantVelocityFrames){
		pp.bulletBody->setLinearVelocity(po.constantVelocity);
	  }
	}
  }
  return;*/


  
  if(currentFrame > 5){
	for(auto & po : plasticBodies){
	  for(auto& constraint: po.constraints){
		std::get<3>(constraint)->setBreakingImpulseThreshold(
			po.breakingThreshold);
	  }
	}
  }
  {
	auto timer = profiler.timeName("bullet step 1");
	bulletWorld->stepSimulation(dt, 10, dt);
  }
  {
	auto time = profiler.timeName("collect impulses");
	collectImpulses();
  }

  {
	auto timer = profiler.timeName("deform and skin");
	//another layer of pararalellism
	tbb::parallel_for(tbb::blocked_range<size_t>(0, plasticBodies.size()),
		[&](const tbb::blocked_range<size_t>& r){
		  
		  bool kernelExists;
		  cl::Kernel& clKernel = clKernels.local(kernelExists);
		  if(!kernelExists){
			clKernel = cl::Kernel(clProgram, "skinVertexVarying");
		  }
		  //for(auto& po : plasticBodies){
		  for(auto i = r.begin(); i != r.end(); ++i){
			
			auto& po = plasticBodies[i];
			po.projectImpulsesOntoExampleManifoldLocally(dt);
			
			po.skinAndUpdate(); //skin pieces and update bullet props
			//po.skinAndUpdateCL(*this, clKernel); //same, but use the opencl skinning code

			//po.updateConstraints(); //make sure the point2point constraints are right
			
			if(po.hasConstantVelocity && currentFrame < po.constantVelocityFrames){
			  for(auto& pp : po.plasticPieces){
				pp.bulletBody->setLinearVelocity(po.constantVelocity);
			  }
			}
			
			for(auto& pp : po.plasticPieces){
			  double avgBcNorm = 
				pp.deltaBarycentricCoordinates.norm()/pp.numPhysicsVertices;
			  
			  //auto f1 = [](double a){ return 1 - 100*a;};
			  auto f2 = [&po](double a){ return exp(-po.restitutionExponent*a);};
			  
			  double restitutionScale = std::max(0.0, f2(avgBcNorm));
			  pp.bulletBody->setRestitution(
				  std::min(po.restitution, 
					  std::min(pp.bulletBody->getRestitution() + dt*po.restitutionAddition , 
						  restitutionScale*po.restitution)));
			  //std::cout << "restitution: " << pp.bulletBody->getRestitution() << std::endl;
			}
			//no breaking in this step
			for(auto& constraint : po.constraints){
			  std::get<3>(constraint)->
				setBreakingImpulseThreshold(std::numeric_limits<double>::infinity());
			}
		  }
		});
  }

  for(auto& rb : rigidBodies){
	if(rb.rbType == RigidBody::RBType::RB_PLANE) continue;
	if (rb.constantForceFrames < 0 || rb.constantForceFrames >  currentFrame) {
	  //rb.bulletBody->applyCentralForce(rb.constantForce);
	  rb.bulletBody->setLinearVelocity(rb.constantForce);
	  std::cout<<rb.constantForce<<std::endl;
	}
  }

  {
	auto timer = profiler.timeName("restitution and constant velocity");
	for(auto& po: plasticBodies){
	  for(auto& pp : po.plasticPieces){
		//if(!singleStep){
		//  pp.bulletBody->setRestitution(po.restitution);
		//}
		if(po.hasConstantVelocity && currentFrame < po.constantVelocityFrames){
		  pp.bulletBody->setLinearVelocity(po.constantVelocity);
		}
		if(pp.plinkoObject && pp.bulletBody->getCenterOfMassPosition().y() < -5){
		  pp.bulletBody->translate(btVector3(0,8, 0));
		}
	  }
	}
  }
}


void World::collectImpulses(){
  for(auto& po : plasticBodies){
	po.manifoldPoints.clear();
  }
  
  for(auto i : range(dispatcher->getNumManifolds())){
	auto* man = dispatcher->getManifoldByIndexInternal(i);

	auto* rb1 = btRigidBody::upcast(man->getBody0());
	auto* rb2 = btRigidBody::upcast(man->getBody1());

	for(auto j : range(man->getNumContacts())){
	  auto& manPoint = man->getContactPoint(j);
	  if(rb1 && rb1->getUserIndex() >= 0){
		plasticBodies[rb1->getUserIndex()].manifoldPoints.emplace_back(
			rb1, manPoint, true);
		
	  }
	  if(rb2 && rb2->getUserIndex() >= 0){
		plasticBodies[rb2->getUserIndex()].manifoldPoints.emplace_back(
			rb2, manPoint, false);
	  }
	}
  }
}

int World::getNumBarrels(){
  int barrelCount = 0;
  for(auto i : range(pyramidSize+1)){
	barrelCount += i*i;
  }
  return barrelCount;
}

void World::makeBarrelPyramid(){
  if(!makeBarrels)
	return;
  int barrelCount = getNumBarrels();
  
  std::cout << "making " << barrelCount << " barrels" << std::endl;
  
  double startHeight = 0.39;
  double currentHeight = startHeight;
  double deltaHeight = 0.77;

  double barrelRadius = 0.34;
  double startX = -3.0, startZ = 0;

  //set up most stuff here...
  std::ifstream ins("inputFiles/barrelStackParams.json");
  Json::Value pbj;
  Json::Reader reader;
  if(!reader.parse(ins, pbj)){
	std::cout << "couldn't read barrel params: " << std::endl
			  << reader.getFormatedErrorMessages();
	exit(1);
  }

  pbj["offset"].resize(3);
  
  for(auto l : range(pyramidSize)){
	auto level = pyramidSize - l;
	for(auto row : range(level)){
	  for(auto col : range(level)){
		
		plasticBodies.emplace_back();
		auto& po = plasticBodies.back();
		
		btVector3 offset{startX + col*2*barrelRadius,
			currentHeight,
			startZ + row*2*barrelRadius};
		
		pbj["offset"][0] = offset.x();
		pbj["offset"][1] = offset.y();
		pbj["offset"][2] = offset.z();
		
		po.loadFromJson(pbj, *this, plasticBodies.size() -1);

	  }
	}
	currentHeight += deltaHeight;
	startX += barrelRadius;
	startZ += barrelRadius;
  }
}


void World::initCL(){

  std::vector<cl::Platform> platforms;
  cl::Platform::get(&platforms);
  std::cout << platforms.size() << std::endl;
 
  std::vector<cl::Device> devices;

  for(auto& plat : platforms){
	std::string s;
	plat.getInfo(CL_PLATFORM_NAME, &s);
	std::cout << "Platform: " << s << std::endl;
	
	plat.getInfo(CL_PLATFORM_VENDOR, &s);
	std::cout << "\tVendor:  " << s << std::endl;
	
	plat.getInfo(CL_PLATFORM_VERSION, &s);
	std::cout << "\tVersion: " << s << std::endl;
	
	// Discover number of devices

	plat.getDevices(CL_DEVICE_TYPE_GPU, &devices);
	std::cout << "\n\tNumber of devices: " << devices.size() << std::endl;

	for (auto& dev : devices){
	  dev.getInfo(CL_DEVICE_NAME, &s);
	  std::cout << "\t\tName: " << s << std::endl;
	    
	  dev.getInfo(CL_DEVICE_OPENCL_C_VERSION, &s);
	  std::cout << "\t\tVersion: " << s << std::endl;
	    
	  int i;
	  dev.getInfo(CL_DEVICE_MAX_COMPUTE_UNITS, &i);
	  std::cout << "\t\tMax. Compute Units: " << i << std::endl;
	  size_t size;
	  dev.getInfo(CL_DEVICE_LOCAL_MEM_SIZE, &size);
	  std::cout << "\t\tLocal Memory Size: " << size/1024 << " KB" << std::endl;
	    
	  dev.getInfo(CL_DEVICE_GLOBAL_MEM_SIZE, &size);
	  std::cout << "\t\tGlobal Memory Size: " << size/(1024*1024) << " MB" << std::endl;
	    
	  dev.getInfo(CL_DEVICE_MAX_MEM_ALLOC_SIZE, &size);
	  std::cout << "\t\tMax Alloc Size: " << size/(1024*1024) << " MB" << std::endl;
	    
	  dev.getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &size);
	  std::cout << "\t\tMax Work-group Total Size: " << size << std::endl;
	    
	  std::vector<size_t> d;
	  dev.getInfo(CL_DEVICE_MAX_WORK_ITEM_SIZES, &d);
	  std::cout << "\t\tMax Work-group Dims: (";
	  for (auto&& st : d)
		std::cout << st << " ";
	  std::cout << "\x08)" << std::endl;

	}

  }

  assert(!devices.empty());
  device = devices[0];
  context = cl::Context(device); //Use the first one
  queue = cl::CommandQueue(context);
  
  std::string programSource = readFile("skinning.cl");
  clProgram = cl::Program(context, programSource, true);
  

}
