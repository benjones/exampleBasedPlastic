
#include "world.h"
#include "rigidBody.h"

#include <fstream>
#include <iostream>
#include <BulletCollision/BroadphaseCollision/btDbvtBroadphase.h>
#include <BulletCollision/CollisionShapes/btStaticPlaneShape.h>
#include <LinearMath/btVector3.h>

#include <BulletCollision/CollisionShapes/btBoxShape.h>
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
  :broadphaseInterface
   (std::unique_ptr<btBroadphaseInterface>(new btDbvtBroadphase())),
   collisionConfiguration
   (std::unique_ptr<btCollisionConfiguration>(new btDefaultCollisionConfiguration())),
   dispatcher
   (std::unique_ptr<btCollisionDispatcher>(new btCollisionDispatcher(collisionConfiguration.get()))),
   bulletSolver
   (std::unique_ptr<btSequentialImpulseConstraintSolver>
    (new btSequentialImpulseConstraintSolver())),
   bulletWorld(dispatcher.get(), 
	       broadphaseInterface.get(), 
	       bulletSolver.get(), 
	       collisionConfiguration.get()),
  currentFrame(0)
{

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

  auto bulletObjectsIn = root["bulletObjects"];


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
	
    bulletWorld.addRigidBody(RB.bulletBody.get());

  }

  dt = root.get("dt", 1.0/60.0).asDouble();
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


  bulletWorld.setGravity(gravity);
  for(auto i : range(bulletObjectsIn.size())){
    std::cout << "adding rb: " << i << std::endl;
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
	  RB.loadTrimesh(bo["filename"].asString());

	}else { //add other shape types here
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
    bulletWorld.addRigidBody(RB.bulletBody.get());
    
  }

  loadPlasticBodies(root);
  makeBarrelPyramid();  

  btGImpactCollisionAlgorithm::registerAlgorithm(dispatcher.get());
}

void World::dumpFrame(){

  char framestring[80];
  sprintf(framestring, "frames/foo-%%03i.%%04i.ply");
  
  std::cout << "writing frame: " << currentFrame << std::endl;

  int objectCount = 0;
  char fname[80];

  for(auto& rigidBody : rigidBodies){
    sprintf(fname, framestring, objectCount, currentFrame);
    rigidBody.dump(fname);
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
	pb.loadFromJson(pbi, bulletWorld, plasticBodies.size() -1 );
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
  
  {
	auto timer = profiler.timeName("save snapshots");
	for(auto& po : plasticBodies){
	  po.saveBulletSnapshots();
	  if(currentFrame > 5){
		for(auto& constraint: po.constraints){
		  std::get<3>(constraint)->setBreakingImpulseThreshold(
			  po.breakingThreshold);
		}
	  }
	}
  }
  {
	auto timer = profiler.timeName("bullet step 1");
	bulletWorld.stepSimulation(dt, 10, dt);
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
		  //  for(auto& po : plasticBodies){
		  for(auto i = r.begin(); i != r.end(); ++i){
			
			auto& po = plasticBodies[i];
			po.projectImpulsesOntoExampleManifoldLocally(dt);
			
			po.skinAndUpdate(); //skin pieces and update bullet props
			po.updateConstraints(); //make sure the point2point constraints are right
			
			
			po.restoreBulletSnapshots();
			
			if(po.hasConstantVelocity){
			  for(auto& pp : po.plasticPieces){
				pp.bulletBody->setLinearVelocity(po.constantVelocity);
			  }
			}
			
			for(auto& pp : po.plasticPieces){
			  double avgBcNorm = 
				pp.deltaBarycentricCoordinates.norm()/pp.activeVertices.size();
			  
			  //auto f1 = [](double a){ return 1 - 100*a;};
			  auto f2 = [](double a){ return exp(-300*a);};
			  
			  double restitutionScale = std::max(0.0, f2(avgBcNorm));
			  pp.bulletBody->setRestitution(
				  std::min(pp.bulletBody->getRestitution(), 
					  restitutionScale*po.restitution));
			  
			}
			//no breaking in this step
			for(auto& constraint : po.constraints){
			  std::get<3>(constraint)->
				setBreakingImpulseThreshold(std::numeric_limits<double>::infinity());
			}
		  }
		});
  }
  {
	auto timer = profiler.timeName("bullet step 2");
	bulletWorld.stepSimulation(dt, 10, dt);
  }
  {
	auto timer = profiler.timeName("restitution and constant velocity");
	for(auto& po: plasticBodies){
	  for(auto& pp : po.plasticPieces){
		pp.bulletBody->setRestitution(po.restitution);
		if(po.hasConstantVelocity){
		  pp.bulletBody->setLinearVelocity(po.constantVelocity);
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
  return;
  int barrelCount = getNumBarrels();
  
  std::cout << "making " << barrelCount << " barrels" << std::endl;
  
  double startHeight = 0.235;
  double currentHeight = startHeight;
  double deltaHeight = 0.77;

  double barrelRadius = 0.34;
  double startX = 0, startZ = 0;

  //set up most stuff here...
  Json::Value pbj;
  pbj["directory"] = "inputFiles/barrel4Examples";
  pbj["density"] = 1000.0;
  pbj["restitution"] = 0.8;
  pbj["plasticityImpulseYield"] = 3e-4;
  pbj["plasticityImpulseScale"] = 200;
  pbj["plasticityKernelScale"] = 0.02;
  pbj["plasticityRate"] = 0.3;
  pbj["jacobianAlpha"] = 1.0;
  pbj["scaleFactor"] = 0.01;
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
		
		po.loadFromJson(pbj, bulletWorld, plasticBodies.size() -1);

	  }
	}
	currentHeight += deltaHeight;
	startX += barrelRadius;
	startZ += barrelRadius;
  }
}

#if 0
//no longer used
void World::loadPlasticObjects(const Json::Value& root){

  plasticObjects.clear();
  auto plasticObjectsIn = root["plasticObjects"];

  std::cout << "slots reserved: " << plasticObjectsIn.size() + getNumBarrels() << std::endl;
  plasticObjects.reserve(plasticObjectsIn.size() + getNumBarrels());
  for(auto i : range(plasticObjectsIn.size())){
	auto& poi = plasticObjectsIn[i];
	
	plasticObjects.emplace_back();
	auto& po = plasticObjects.back();
	//load meshes
	po.loadFromFiles(poi["directory"].asString());


	po.dt = dt;
	po.density = poi.get("density", 1000).asDouble();
	
	po.plasticityImpulseYield = poi.get("plasticityImpulseYield", 1e6).asDouble();
	po.plasticityImpulseScale = poi.get("plasticityImpulseScale", 0).asDouble();
	po.plasticityKernelScale = poi.get("plasticityKernelScale", 1.0).asDouble();
	po.plasticityRate = poi.get("plasticityRate", 1.0).asDouble();
	po.jacobianAlpha = poi.get("jacobianAlpha", 1.0).asDouble();
	

	po.localPlasticityImpulseYield = poi.get("localPlasticityImpulseYield", 1e6).asDouble();
	po.localPlasticityImpulseScale = poi.get("localPlasticityImpulseScale", 0).asDouble();

	auto& pesIn = poi["perExampleScale"];
	if(!pesIn.isNull()){
	  for(auto j : range(pesIn.size())){
		po.perExampleScale[j] = pesIn[j].asDouble();
	  }
	}


	po.scaleFactor = poi.get("scaleFactor", 1).asDouble();
	

	po.currentBulletVertexPositions = po.scaleFactor*po.tetmeshVertices;

	//po.localImpulseBasedOffsets = RMMatrix3d::Zero(po.numPhysicsVertices, 3);
	
	//	po.currentBulletVertexPositions.resize( po.numPhysicsVertices,
	//											Eigen::NoChange);
	
	//make space for the vertex array
	po.btTriMesh = std::unique_ptr<btTriangleIndexVertexArray>{
	  new btTriangleIndexVertexArray{
		static_cast<int>(po.tetmeshTriangles.rows()),
		po.tetmeshTriangles.data(),
		3*sizeof(int), //UGH, THIS IS TURRRRRIBLE
		static_cast<int>(po.currentBulletVertexPositions.rows()),
		po.currentBulletVertexPositions.data(),
		3*sizeof(double), //UGH THIS IS ALSO TURRRRIBLE
	  }
	};
	
	std::cout << "nverts: " << po.currentBulletVertexPositions.rows() << "  ntris: " 
			  << po.tetmeshTriangles.rows() << std::endl;

	//	std::cout << "first 3 vertices: " << 
	//	  po.currentBulletVertexPositions.block(0,0, 3, 3) << std::endl;
	/*for(auto row : range(po.tetmeshTriangles.rows())){
	  std::cout  << "row: " << row << " ";
	  for(auto col : range(po.tetmeshTriangles.cols())){
		
		std::cout << po.tetmeshTriangles.data()[row*3 + col] << ' ';
	  }
	  std::cout << std::endl;
	  }*/

	
	



	//apply COM offset stuff
	btVector3 offset(0.0,0.0,0.0); //default to no offset
	auto offsetIn = poi["offset"];
	if(!offsetIn.isNull() && offsetIn.isArray()){
	  assert(offsetIn.size() == 3);
	  offset = btVector3{offsetIn[0].asDouble(),
						 offsetIn[1].asDouble(),
						 offsetIn[2].asDouble()};
	}
	
	auto rotation = btQuaternion::getIdentity();
	auto rotationIn = poi["rotation"];
	if(!rotationIn.isNull()){
	  btVector3 axis{rotationIn["axis"][0].asDouble(),
		  rotationIn["axis"][1].asDouble(),
		  rotationIn["axis"][2].asDouble()};
	  rotation = btQuaternion{axis,
							  rotationIn["angle"].asDouble()};
	  
	  std::cout << "rotation: " << rotation << std::endl;
	}

	auto constantVelocityIn = poi["constantVelocity"];
	if(!constantVelocityIn.isNull() && 
	   constantVelocityIn.isArray() && 
	   constantVelocityIn.size() == 3){
	  
	  po.hasConstantVelocity = true;
	  po.constantVelocity = btVector3{constantVelocityIn[0].asDouble(),
									  constantVelocityIn[1].asDouble(),
									  constantVelocityIn[2].asDouble()};
	  

	}

	//place it where it should go
	//this isn't mass weighted, but is probably close if the mesh
	//close to uniform

	//actually, just jam this into the worldTransform 
	/*Eigen::Vector3d centerOfMass = 
	  po.tetmeshVertices.colwise().sum().transpose()/
	  po.numPhysicsVertices;
	

	for(auto i : range(po.currentBulletVertexPositions.rows())){
	  po.currentBulletVertexPositions.row(i) =
		(rotation*(po.tetmeshVertices.row(i).transpose() -
				   centerOfMass) +
		 offset).transpose();
	}
	*/




	po.bulletShape = std::unique_ptr<btGImpactMeshShape>{
	  new btGImpactMeshShape{po.btTriMesh.get()}};
	

	//po.compoundShape = std::unique_ptr<btCompoundShape>{
	//	  new btCompoundShape{}};

	po.motionState = 
	  std::unique_ptr<btDefaultMotionState>{
	  new btDefaultMotionState{}};
	

	
	po.mass = 1; //save a bit of grief in rb construction
	po.bulletBody = std::unique_ptr<btRigidBody>{
	  new btRigidBody{po.mass,
					  po.motionState.get(),
					  po.bulletShape.get()}
					  //po.compoundShape.get()}
	};

	//compute node masses and volume
	po.computeMassesAndVolume();

	Eigen::Vector3d com = 
	  (po.tetmeshVertexMasses.asDiagonal()*
	   po.currentBulletVertexPositions).colwise().sum()/
	  po.mass;
	  

	offset -= quatRotate(rotation, eigenToBullet(com));
	


	//apply the transform now
	po.inertiaAligningTransform.setIdentity();
	po.worldTransform = btTransform{rotation, offset};
	po.bulletBody->setCenterOfMassTransform(po.worldTransform);
	
	po.saveBulletSnapshot(); //to get the COM for updateBulletProperties
	
	//aliasing issues?
	po.updateBulletProperties(po.currentBulletVertexPositions,
							  po.tetmeshTets);

	//once we know the inertia tensor and stuff
	po.saveBulletSnapshot(); //to save a snapshot with correct inertia

	//po.updateCompoundShape();
	

	po.bulletShape->updateBound();
	//	po.compoundShape->updateBound();
	
	
	std::cout << "po mass: " << 1.0/po.bulletBody->getInvMass() << std::endl;
	std::cout << "po inertia inverse: " << po.bulletBody->getInvInertiaTensorWorld() << std::endl;

	po.restitution = poi.get("restitution", 0.8).asDouble();
	po.minRestitution = po.restitution;
	po.bulletBody->setRestitution(po.restitution);

	bulletWorld.addRigidBody(po.bulletBody.get());
	//user index holds index into the plasticObjects array
	po.bulletBody->setUserIndex(plasticObjects.size() -1);
	
	//	po.egTraverser.restPosition = EGPosition{0, {0,1, 0}};
	//	po.egTraverser.currentPosition = EGPosition{0, {1,0, 0}};
	//	po.egTraverser.restSpringStrength = 0.00;//1;

  }
  
}
#endif
