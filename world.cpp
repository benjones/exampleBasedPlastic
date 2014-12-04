
#include "world.h"
#include "rigidBody.h"
//#include "fem.H"

#include <fstream>
#include <iostream>
#include <BulletCollision/BroadphaseCollision/btDbvtBroadphase.h>
#include <BulletCollision/CollisionShapes/btStaticPlaneShape.h>
#include <LinearMath/btVector3.h>

#include <BulletCollision/CollisionShapes/btBoxShape.h>
#include <BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h>

#include "json/json.h"
//#include "cppitertools/range.hpp"
//#include "cppitertools/enumerate.hpp"
#include "utils.h"
#include <numeric>
#include <limits>

/*
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
extern "C" {
#include <cblas.h>
}
#endif
*/


//#include "tminres.hpp"
//#include "SimpleVector.hpp"

//using iter::range;
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
  currentFrame(0)//,
  //	FEMS(false),
  //	RIGIDS(false)
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
  //  auto femObjectsIn = root["femObjects"];
  //  if (bulletObjectsIn.size() > 0) RIGIDS = true;
  //  if (femObjectsIn.size() > 0) FEMS = true;

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


  //  regularizerAlpha = root.get("regularizerAlpha", 1.0).asDouble();

  /*
  femObjects.reserve(femObjectsIn.size());
  for(auto i : range(femObjectsIn.size())){
		auto& femObjectIn = femObjectsIn[i];
		MaterialProperties mp;
		mp.density = femObjectIn.get("density", 1000.0).asDouble();
		mp.lambda = femObjectIn.get("lambda", 5e4).asDouble();
		mp.mu = femObjectIn.get("mu", 1e5).asDouble();
		mp.alphaCap = femObjectIn.get("alphaCap", 0.0).asDouble();
		//		mp.scale = femObjectIn.get("scale", 1e-2).asDouble();
		mp.scale = femObjectIn.get("dampingscale", 6e-1).asDouble();
		mp.yieldStress = femObjectIn.get("yieldStress", DBL_MAX).asDouble();
		mp.plasticModulus = femObjectIn.get("plasticModulus", 0.0).asDouble();
		mp.flowrate = femObjectIn.get("flowrate", 0.0).asDouble();
		mp.toughness = femObjectIn.get("toughness", DBL_MAX).asDouble();
    auto fnameIn = femObjectIn["filename"];
    femObjects.emplace_back();
		auto& femObj = femObjects.back();
		femObj.load(fnameIn.asString().c_str(), mp, &bulletWorld);
		if (RIGIDS) femObj.setupBulletParticles();

		auto centerOfMass = SlVector3{0, 0, 0};
		double totalMass = 0.0;
		for(auto vInd : range(femObj.nv)){
			centerOfMass += femObj.mass[vInd]*femObj.pos[vInd];
			totalMass += femObj.mass[vInd];
		}
		centerOfMass /= totalMass;
		//center the object about its COM to make rotation make sense
		femObj.bbMin -= centerOfMass;
		femObj.bbMax -= centerOfMass;
		for(auto  vInd : range(femObj.nv)){
			femObj.pos[vInd] -= centerOfMass;
		}
		std::cout << femObj.bbMin << femObj.bbMax << std::endl;
		std::cout<<mp<<std::endl;

	if(!femObjectIn["rotation"].isNull()){
	  SlVector3 axis;
	  axis(0) = femObjectIn["rotation"]["axis"][0].asDouble();
	  axis(1) = femObjectIn["rotation"]["axis"][1].asDouble();
	  axis(2) = femObjectIn["rotation"]["axis"][2].asDouble();
	  normalize(axis);
	  double angle = femObjectIn["rotation"]["angle"].asDouble();
	  
	  std::cout << "rotating about " << axis << " by " 
				<< angle << " radians" << std::endl;

	  SlMatrix3x3 rot;
	  SlMomentToPosMatrix(angle*axis, rot);
	  
	  //recompute bounding boxes
	  femObj.bbMin.set(std::numeric_limits<double>::max());
	  femObj.bbMax.set(std::numeric_limits<double>::lowest());
	  for(auto vInd : range(femObj.nv)){
			femObj.pos[vInd] = rot*femObj.pos[vInd];
			for(auto j : range(3)){
				femObj.bbMin[j] = std::min(femObj.bbMin[j], femObj.pos[vInd][j]);
				femObj.bbMax[j] = std::max(femObj.bbMax[j], femObj.pos[vInd][j]);
			}
	  }
		std::cout<<"new bounding box "<<femObj.bbMin<<"x"<<femObj.bbMax<<std::endl;
	}	  

	

    if(!femObjectIn["offset"].isNull()){
      assert(femObjectIn["offset"].isArray() &&
	     femObjectIn["offset"].size() == 3);
      auto offset = SlVector3{femObjectIn["offset"][0].asDouble(),
							  femObjectIn["offset"][1].asDouble(),
							  femObjectIn["offset"][2].asDouble()};
	  
			double scale = femObjectIn.get("scale", 1.0).asDouble();

      femObj.bbMin += offset;
      femObj.bbMax += offset;

      femObj.bbMin *= scale;
      femObj.bbMax *= scale;

			std::cout<<"offsetting by: "<<offset<<std::endl;
			std::cout<<"scaling by: "<<scale<<std::endl;
      for(auto  vInd : range(femObj.nv)){
				femObj.pos[vInd] += offset;
				femObj.pos[vInd] *= scale;
      }
      
      if(i == 0){
				femObj.firstNodeIndex = 0;
      } else {
				femObj.firstNodeIndex = femObjects[i-1].firstNodeIndex +
					femObjects[i -1].nv;
      }

			std::cout<<"new bounding box "<<femObj.bbMin<<"x"<<femObj.bbMax<<std::endl;
		}
	 
		femObj.initializeRestState();


		if(!femObjectIn["velocity"].isNull()){
			assert(femObjectIn["velocity"].isArray() &&
						 femObjectIn["velocity"].size() == 3);
			auto velocity = SlVector3{femObjectIn["velocity"][0].asDouble(),
																femObjectIn["velocity"][1].asDouble(),
																femObjectIn["velocity"][2].asDouble()};
			
	      for(auto  vInd : range(femObj.nv)){
					femObj.vel[vInd] = velocity;
					//std::cout<<"setting velocity "<<velocity<<std::endl;
				}
		}
		
  }
  */
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
  std::cout << "read in " << rigidBodies.size() << " rbs" << std::endl;


  loadPlasticObjects(root);
  makeBarrelPyramid();  

  //  computeConstraints();
  //  countConstraints();
  /*  
  std::function<void(btBroadphasePair&,
					 btCollisionDispatcher&, 
					 const btDispatcherInfo&)> collisionCallback = 
	[this](btBroadphasePair& collisionPair,
		   btCollisionDispatcher& _dispatcher,
		   const btDispatcherInfo& dispatcherInfo){
	constraintCullingCallback(collisionPair, _dispatcher, dispatcherInfo);
  };

  dispatcher.get()->setNearCallback(collisionCallback);
  */
  btGImpactCollisionAlgorithm::registerAlgorithm(dispatcher.get());
}

void World::dumpFrame(){

  char framestring[80];
  char bcsstring[80];
  char impulsestring[80];
  sprintf(framestring, "frames/foo-%%03i.%%04i.ply");
  sprintf(bcsstring, "frames/foo-%%03i.%%04i.bcs");
  sprintf(impulsestring, "frames/foo-%%03i.%%04i.impulses");
  std::cout << "writing frame: " << currentFrame << std::endl;

  int objectCount = 0;
  char fname[80];
  /*  for(auto& femObject : femObjects){
    sprintf(fname, framestring, objectCount, currentFrame);
    femObject.dumpObj(fname);
    objectCount++;
	}*/
  for(auto& rigidBody : rigidBodies){
    sprintf(fname, framestring, objectCount, currentFrame);
    rigidBody.dump(fname);
	objectCount++;
  }
  
  for(auto& plasticObject : plasticObjects){

	/*sprintf(fname, framestring, objectCount, currentFrame);
	objectCount++;
	std::ofstream outs(fname);
	
	btTransform trans = plasticObject.bulletBody->getCenterOfMassTransform();
	btVector3 abmin, abmax;
	plasticObject.bulletShape->getAabb(trans, abmin, abmax);
	outs << "v " << abmin.x() << ' ' << abmin.y() << ' ' << abmin.z() << '\n'
		 << "v " << abmax.x() << ' ' << abmin.y() << ' ' << abmin.z() << '\n'
		 << "v " << abmin.x() << ' ' << abmax.y() << ' ' << abmin.z() << '\n'
		 << "v " << abmax.x() << ' ' << abmax.y() << ' ' << abmin.z() << '\n'
		 << "v " << abmin.x() << ' ' << abmin.y() << ' ' << abmax.z() << '\n'
		 << "v " << abmax.x() << ' ' << abmin.y() << ' ' << abmax.z() << '\n'
		 << "v " << abmin.x() << ' ' << abmax.y() << ' ' << abmax.z() << '\n'
		 << "v " << abmax.x() << ' ' << abmax.y() << ' ' << abmax.z() << std::endl;
	
	trans.setIdentity();
	plasticObject.bulletShape->getAabb(trans, abmin, abmax);
	for(auto i : range(plasticObject.currentBulletVertexPositions.rows())){
	  
	  if(plasticObject.currentBulletVertexPositions(i, 0) < abmin.x() ||
		 plasticObject.currentBulletVertexPositions(i, 1) < abmin.y() ||
		 plasticObject.currentBulletVertexPositions(i, 2) < abmin.z() ||
		 plasticObject.currentBulletVertexPositions(i, 0) > abmax.x() ||
		 plasticObject.currentBulletVertexPositions(i, 1) > abmax.y() ||
		 plasticObject.currentBulletVertexPositions(i, 2) > abmax.z()){
		
		std::cout << "vertex "<< i << " outside of bounding box :( " << std::endl;
		std::cout << abmin << "\n"
				  << abmax << "\n"
				  << plasticObject.currentBulletVertexPositions.row(i) << std::endl;


		Eigen::Vector3d realMin = plasticObject.currentBulletVertexPositions.colwise().minCoeff();
		Eigen::Vector3d realMax = plasticObject.currentBulletVertexPositions.colwise().maxCoeff();
		
		std::cout << "eigen: " << realMin << std::endl << realMax << std::endl;

		std::cout << "bb size: " << abmax - abmin << std::endl;
		std::cout << "eigen: " << realMax - realMin << std::endl;


		break;
		//exit(1);
	  }

	  }*/
	/*if(currentFrame == 0){
	  sprintf(fname, framestring, objectCount, currentFrame);
	  plasticObject.dump(fname);
	  }*/
	sprintf(fname, framestring, objectCount, currentFrame);
	plasticObject.dump(fname);
	sprintf(fname, bcsstring, objectCount, currentFrame);
	plasticObject.dumpBarycentricCoords(fname);
	sprintf(fname, impulsestring, objectCount, currentFrame);
	plasticObject.dumpImpulses(fname);
	objectCount++;
	  
  }

  /*{
	sprintf(fname, framestring, objectCount, currentFrame);
	std::ofstream outs(fname);
	
	
	objectCount++;
	
	for(auto i : range(dispatcher->getNumManifolds())){
	  auto* man = dispatcher->getManifoldByIndexInternal(i);
	  for(auto j : range(man->getNumContacts())){
		auto& manPoint = man->getContactPoint(j);
		auto& va = manPoint.getPositionWorldOnA();
		auto& vb = manPoint.getPositionWorldOnB();
		outs << "v " << va.x() << ' ' << va.y() << ' ' << va.z() << std::endl;
		outs << "v " << vb.x() << ' ' << vb.y() << ' ' << vb.z() << std::endl;
	  }
	}
	}*/
  currentFrame++;

}



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

	po.localImpulseBasedOffsets = RMMatrix3d::Zero(po.numPhysicsVertices, 3);
	
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

void World::timeStepDynamicSpritesNoDouble(){


  bulletWorld.stepSimulation(dt, 10, dt);
  for(auto& po : plasticObjects){
	po.saveBulletSnapshot();
  }
  
  collectImpulses();

  deformBasedOnImpulses();
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

}

void World::timeStepDynamicSprites(){

  for(auto& po : plasticObjects){
	po.saveBulletSnapshot();
	po.deformedThisFrame = currentFrame == 1 || (po.deltaBarycentricCoordinates.norm() > 0);
  }
  //  std::cout << "do fist step" << std::endl;
  bulletWorld.stepSimulation(dt, 10, dt);

  //  std::cout << "deform" << std::endl;  
  collectImpulses();

  deformBasedOnImpulses();

  for(auto& po : plasticObjects){
	
	//	if(po.projectImpulsesOntoExampleManifold()){
	if(po.projectImpulsesOntoExampleManifoldLocally()){
	  po.deformedThisFrame = true;
	}

	if(po.deformedThisFrame){
	  po.skinMesh();//bulletWorld);
	  po.updateBulletProperties(po.currentBulletVertexPositions, 
								po.tetmeshTets);
	}
	//po.updateCompoundShape();
	po.restoreBulletSnapshot();
	
	if(po.hasConstantVelocity){
	  po.bulletBody->setLinearVelocity(btVector3{po.constantVelocity.x(),
			po.bulletBody->getLinearVelocity().y(),
			po.constantVelocity.z()});
	}
	
	double avgBcNorm = po.deltaBarycentricCoordinates.norm()/po.numPhysicsVertices;
	
	auto f1 = [](double a){ return 1 - 100*a;};
	auto f2 = [](double a){ return exp(-300*a);};

	double restitutionScale = std::max(0.0, f2(avgBcNorm));
	//	if(restitutionScale != 1){ std::cout << "restitution scale: " << restitutionScale << '\n';}
	po.minRestitution = std::min(po.minRestitution, po.restitution*restitutionScale); 
	po.bulletBody->setRestitution(po.minRestitution);

  }
  
  auto minIt = 
	std::min_element(plasticObjects.begin(),
					 plasticObjects.end(),
					 [](const PlasticObject& a,
						const PlasticObject& b){
					   return a.bulletBody->getRestitution() <
					   b.bulletBody->getRestitution();
					 });
  
  std::cout << "min restitution scale: " 
	//<< "total : " << minIt->bulletBody->getRestitution()
	//		<< "ratio: "
			<< minIt->bulletBody->getRestitution()/minIt->restitution
			<< '\n';


  //std::cout << "do second step: " << std::endl;
  bulletWorld.stepSimulation(dt, 10, dt);
  //std::cout << "step done" << std::endl;

  for(auto& po: plasticObjects){
	po.bulletBody->setRestitution(po.restitution);
	if(po.hasConstantVelocity){
	  po.bulletBody->setLinearVelocity(btVector3{po.constantVelocity.x(),
			po.bulletBody->getLinearVelocity().y(),
			po.constantVelocity.z()});
	}
  }

  //deformBasedOnImpulses();

  

  //  for(auto& po : plasticObjects){
	


	//po.egTraverser.traverse();
	//po.skinMesh();
	//po.currentBulletVertexPositions = po.tetmeshVertices;
	//po.updateBulletProperties(po.currentBulletVertexPositions,
	//							  po.tetmeshTets);
  //  }
}


void World::deformBasedOnImpulses(){
  //  std::cout << "num contacts this frame: " << 
  //	dispatcher->getNumManifolds() << std::endl;
  for(auto i : range(dispatcher->getNumManifolds())){
	auto* man = dispatcher->getManifoldByIndexInternal(i);
	
	auto* rb1 = btRigidBody::upcast(man->getBody0());
	auto* rb2 = btRigidBody::upcast(man->getBody1());
	if(rb1 && rb1->getUserIndex() >= 0){
	  //	  std::cout << "first" << std::endl;
	  //std::cout << "index: " << rb1->getUserIndex() << std::endl;
	  if(plasticObjects[rb1->getUserIndex()].deformBasedOnImpulseLocal(man, true)){
		plasticObjects[rb1->getUserIndex()].deformedThisFrame = true;
	  }
	}
	if(rb2 && rb2->getUserIndex() >= 0){
	  //	  std::cout << "second" << std::endl;
	  //std::cout << "index: " << rb2->getUserIndex() << std::endl;
	  if(plasticObjects[rb2->getUserIndex()].deformBasedOnImpulseLocal(man, false)){
		plasticObjects[rb2->getUserIndex()].deformedThisFrame = true;
	  }
	}
  }
}

void World::collectImpulses(){
  for(auto& po : plasticObjects){
	po.manifoldPoints.clear();
  }

  for(auto i : range(dispatcher->getNumManifolds())){
	auto* man = dispatcher->getManifoldByIndexInternal(i);

	auto* rb1 = btRigidBody::upcast(man->getBody0());
	auto* rb2 = btRigidBody::upcast(man->getBody1());


	for(auto j : range(man->getNumContacts())){
	  auto& manPoint = man->getContactPoint(j);
	  if(rb1 && rb1->getUserIndex() >= 0){
		//std::cout << "index: " << rb1->getUserIndex() << std::endl;
		plasticObjects[rb1->getUserIndex()].manifoldPoints.push_back(std::make_pair(manPoint,
																					true));
	  }
	  if(rb2 && rb2->getUserIndex() >= 0){
		//std::cout << "index: " << rb2->getUserIndex() << std::endl;
		plasticObjects[rb2->getUserIndex()].manifoldPoints.push_back(std::make_pair(manPoint,
																					false));
	  }
	}
	//man->clearManifold();
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

  int barrelCount = getNumBarrels();
  
  std::cout << "making " << barrelCount << " barrels" << std::endl;
  
  //plasticObjects.clear();
  //plasticObjects.reserve(barrelCount);

  double startHeight = 0.235;
  double currentHeight = startHeight;
  double deltaHeight = 0.77;

  double barrelRadius = 0.34;
  double startX = 0, startZ = 0;

  for(auto l : range(pyramidSize)){
	auto level = pyramidSize - l;
	for(auto row : range(level)){
	  for(auto col : range(level)){
		
		plasticObjects.emplace_back();
		auto& po = plasticObjects.back();
		
		po.loadFromFiles("inputFiles/barrel4Examples");
		po.dt = dt;
		po.density = 1000;
		
		po.plasticityImpulseYield = 3e-4;//0.0004;//0.001;
		po.plasticityImpulseScale = 200;//80;
		po.plasticityKernelScale = 0.02;
		po.plasticityRate = 0.3;
		po.jacobianAlpha = 1.0;


		po.perExampleScale(0) = 1;
		po.perExampleScale(1) = 1;
		po.perExampleScale(2) = 1;
		po.perExampleScale(3) = 1;

		po.localPlasticityImpulseScale = 0;//3;
		po.localPlasticityImpulseYield = 0.0001;
		


		po.scaleFactor = 0.01;
		
		po.currentBulletVertexPositions = po.scaleFactor*po.tetmeshVertices;

		po.localImpulseBasedOffsets = RMMatrix3d::Zero(po.numPhysicsVertices, 3);

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
		
		btVector3 offset{startX + col*2*barrelRadius,
			currentHeight,
			startZ + row*2*barrelRadius};

		
		po.bulletShape = std::unique_ptr<btGImpactMeshShape>{
		  new btGImpactMeshShape{po.btTriMesh.get()}};
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


		//apply the transform now
		po.inertiaAligningTransform.setIdentity();
		po.worldTransform = btTransform{btQuaternion::getIdentity(), 
										offset};
		po.bulletBody->setCenterOfMassTransform(po.worldTransform);
		
		po.saveBulletSnapshot();
	
		po.restitution = 0.8;
		po.minRestitution = po.restitution;
		po.bulletBody->setRestitution(po.restitution);

		//aliasing issues?
		po.updateBulletProperties(po.currentBulletVertexPositions,
								  po.tetmeshTets);

		//po.updateCompoundShape();
	

		po.bulletShape->updateBound();

		
		bulletWorld.addRigidBody(po.bulletBody.get());
		//user index holds index into the plasticObjects array
		po.bulletBody->setUserIndex(plasticObjects.size() -1);

		//		po.egTraverser.restPosition = 
		//		  EGPosition{0, std::vector<double>(po.exampleGraph.simplices[0].size(), 0.0)};
		//		po.egTraverser.restPosition.coords[0] = 1;
		//		po.egTraverser.currentPosition = po.egTraverser.restPosition;
		//		po.egTraverser.restSpringStrength = 0.0;
		

	  }
	}

	
	currentHeight += deltaHeight;
	startX += barrelRadius;
	startZ += barrelRadius;
  }


}
