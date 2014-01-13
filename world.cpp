
#include "world.h"
#include "rigidBody.h"
#include "fem.H"

#include <fstream>
#include <iostream>
#include <BulletCollision/BroadphaseCollision/btDbvtBroadphase.h>
#include <BulletCollision/CollisionShapes/btStaticPlaneShape.h>
#include <LinearMath/btVector3.h>
#include <LinearMath/btDefaultMotionState.h>
#include <BulletCollision/CollisionShapes/btBoxShape.h>

#include "json/json.h"
#include "cppitertools/range.hpp"
#include "cppitertools/enumerate.hpp"
#include "utils.h"
#include <numeric>
#include <limits>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#elif linux
extern "C" {
#include <cblas.h>
}
#endif

#include "tminres.hpp"
#include "SimpleVector.hpp"

using iter::range;
using iter::enumerate;

void callbackWrapper(btBroadphasePair& collisionPair, 
					 btCollisionDispatcher& _dispatcher, 
					 btDispatcherInfo& dispatcherInfo){
  
}

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


  regularizerAlpha = root.get("regularizerAlpha", 1.0).asDouble();


  auto femObjectsIn = root["femObjects"];
  femObjects.reserve(femObjectsIn.size());
  for(auto i : range(femObjectsIn.size())){
    auto& femObjectIn = femObjectsIn[i];
    auto fnameIn = femObjectIn["filename"];
    femObjects.emplace_back();
	auto& femObj = femObjects.back();
    femObj.load(fnameIn.asString().c_str(), &bulletWorld);

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
		  femObj.bbMin(j) = std::min(femObj.bbMin(j), femObj.pos[vInd](j));
		  femObj.bbMax(j) = std::max(femObj.bbMax(j), femObj.pos[vInd](j));
		}
	  }
	  
	}
	

    if(!femObjectIn["offset"].isNull()){
      assert(femObjectIn["offset"].isArray() &&
	     femObjectIn["offset"].size() == 3);
      auto offset = SlVector3{femObjectIn["offset"][0].asDouble(),
							  femObjectIn["offset"][1].asDouble(),
							  femObjectIn["offset"][2].asDouble()};
	  
      femObj.bbMin += offset;
      femObj.bbMax += offset;
      for(auto  vInd : range(femObj.nv)){
				femObj.pos[vInd] += offset;
      }
      
      if(i == 0){
				femObj.firstNodeIndex = 0;
      } else {
				femObj.firstNodeIndex = femObjects[i-1].firstNodeIndex +
					femObjects[i -1].nv;
      }

    }
    
  }

  auto bulletObjectsIn = root["bulletObjects"];
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

    
    auto shapeTypeIn = bo["shape"];
    if(shapeTypeIn.asString() == "box"){
      RB.rbType = RigidBody::RBType::RB_BOX;
      auto extentsIn = bo["extents"];
      if(!extentsIn.isArray()){
	std::cout << "boxest must have extents, half widths of the box" << std::endl;
	exit(1);
      }
      btVector3 extents = btVector3(extentsIn[0].asDouble(),
				    extentsIn[1].asDouble(),
				    extentsIn[2].asDouble());
      
      RB.shape = std::unique_ptr<btCollisionShape>(new btBoxShape(extents));
      
    } else { //add other shape types here
      std::cout << "unknown shape: " << bo["shape"].asString() << std::endl;
      exit(1);
    }
    btVector3 inertiaTensor(0,0,0);

    RB.shape->calculateLocalInertia(mass, inertiaTensor);
    RB.bulletBody = std::unique_ptr<btRigidBody> (new btRigidBody(mass, 
								  RB.motionState.get(), 
								  RB.shape.get(),
								  inertiaTensor));

      
    bulletWorld.addRigidBody(RB.bulletBody.get());
    
  }
  std::cout << "read in " << rigidBodies.size() << " rbs" << std::endl;
  computeConstraints();
  countConstraints();
  
  std::function<void(btBroadphasePair&,
					 btCollisionDispatcher&, 
					 const btDispatcherInfo&)> collisionCallback = 
	[this](btBroadphasePair& collisionPair,
		   btCollisionDispatcher& _dispatcher,
		   const btDispatcherInfo& dispatcherInfo){
	constraintCullingCallback(collisionPair, _dispatcher, dispatcherInfo);
  };
  dispatcher->setNearCallback(collisionCallback);

}

void World::timeStep(){
 //replace these two with a single CG solve

  //make sure to call countConstraints every time the number of constraints changes
  //make sure to call computeCrossProductMatrices before calling the JV and 
  //JTLambda multiply functions
	computeCrossProductMatrices();

  computeFemVelocities(); // this is really just forces
  //bulletWorld.stepSimulationVelocitiesOnly(dt); 

  

  for(auto& fem : femObjects){
	fem.copyStateToBulletParticles();
  }
  bulletWorld.preCoupledSolve(dt);

  //apply constant forces:
  for(auto& rb : rigidBodies){
	rb.bulletBody->applyCentralForce(rb.constantForce);
  }

  //std::cout << "pre solve" << std::endl;
  //for(auto& fem: femObjects) {fem.dump();}
  solve();

  //std::cout << "post solve" << std::endl;
  //for(auto& fem: femObjects) {fem.dump();}
	//bulletWorld.stepSimulation(dt);

  for(auto& fem : femObjects){
	fem.copyStateToBulletParticles();
	//std::cout <<"pre postcopule " << std::endl;
	//fem.dump();
  }
  bulletWorld.postCoupledSolve(dt);
  
  //updateFemPositions();
  for(auto& fem : femObjects){
	fem.copyVelocitiesFromBulletParticles();
	fem.updatePositions(dt);
	fem.fracture();
	//std::cout << "updating pos" << std::endl;
	//fem.dump();
  }

}

void World::timeStepRigidCollisions(){


  //replace these two with a single CG solve

  //make sure to call countConstraints every time the number of constraints changes
  //make sure to call computeCrossProductMatrices before calling the JV and 
  //JTLambda multiply functions
	computeCrossProductMatrices();

  computeFemVelocities(); // this is really just forces
  //bulletWorld.stepSimulationVelocitiesOnly(dt); 
  bulletWorld.preCoupledSolve(dt);

  solve();

	//bulletWorld.stepSimulation(dt);
  for(auto& fem : femObjects){
	fem.updateBulletShapes();
  }
  bulletWorld.postCoupledSolve(dt);

  //updateFemPositions();
  for(auto& fem : femObjects){
	fem.stitchTets();
	fem.fracture();
  }
}


void World::solve() {
	double delta_new, delta_old, alpha, beta, *dptr1, *dptr2, *dptr3;
	double tol=1e-40;
	unsigned int iter=0, max_iter = 10000;
	nfemdof = 0;
	nrbdof = 0;
	unsigned int nrbs = 0;

	unsigned int n = 0; 
	for (unsigned int i=0; i<femObjects.size(); i++) {
		femObjects[i].firstNodeIndex = n/3;
		n += 3*femObjects[i].nv;
		nfemdof += 3*femObjects[i].nv;
	}

	for (unsigned int i=0; i<rigidBodies.size(); i++) {
		if (rigidBodies[i].rbType == RigidBody::RBType::RB_PLANE) continue;
		n += 6;
		nrbs++;
		nrbdof += 6;
	}

	double *work = new double [4*n + nfemdof];
	double *x = work;     // solution
	double *r = work+n;   // residual
	double *q = work+2*n; // temp storage
	double *d = work+3*n; // search direction
	double *p = work+4*n; // preconditioner
	
	for (unsigned int i=0; i<4*n+nfemdof; i++) work[i] = 0.0;

	// set up fem preconditioner
	int offset = 0;
	for (unsigned int i=0; i<femObjects.size(); i++) {
		femObjects[i].gm->setPrecon(p+offset);
		offset += 3*femObjects[i].nv;
	}

	// zero out initial guess
	dptr1 = x;
	for (unsigned int i=0; i<n; i++, dptr1++) (*dptr1) = 0.0;

	dptr1 = r;
	// fill in rhs (momentum + time integrated forces), fem first
	for (unsigned int i=0; i<femObjects.size(); i++) {
		SlVector3 *vptr = femObjects[i].frc;
		for (unsigned int j=0; j<femObjects[i].nv; j++, vptr++) {
			(*(dptr1++)) = (*vptr)[0]; (*(dptr1++)) = (*vptr)[1]; (*(dptr1++)) = (*vptr)[2];
		}
	}

	// rigid bodies
	for (unsigned int i=0; i<rigidBodies.size(); i++) {
		if (rigidBodies[i].rbType == RigidBody::RBType::RB_PLANE) continue;
		double mass = 1.0 / rigidBodies[i].bulletBody->getInvMass();
    (*(dptr1++)) = mass*rigidBodies[i].bulletBody->getLinearVelocity()[0];
    (*(dptr1++)) = mass*rigidBodies[i].bulletBody->getLinearVelocity()[1];
    (*(dptr1++)) = mass*rigidBodies[i].bulletBody->getLinearVelocity()[2];

    //rotational part, there must be a less hideous way of doing this.
		double angularvelocity[3];
		angularvelocity[0] = rigidBodies[i].bulletBody->getAngularVelocity()[0];
		angularvelocity[1] = rigidBodies[i].bulletBody->getAngularVelocity()[1];
		angularvelocity[2] = rigidBodies[i].bulletBody->getAngularVelocity()[2];
    bulletOverwriteMultiply(angularvelocity, dptr1,
			    rigidBodies[i].bulletBody->getInvInertiaTensorWorld().inverse());
		dptr1 += 3;
	}

	// constraint springs, double check that these are the right rigid body positions
	size_t rbIndex = nfemdof;
  for(auto& rb : rigidBodies){
		if (rb.rbType == RigidBody::RBType::RB_PLANE) continue;
    for(auto& c : rb.constraints){
      size_t femIndex = 3*(femObjects[c.femIndex].firstNodeIndex + c.nodeIndex);
			btVector3 f, t;
			SlVector3 &fempos = femObjects[c.femIndex].pos[c.nodeIndex];
			btVector3 rbpos = rb.bulletBody->getCenterOfMassTransform()*c.localPosition;
			f[0] = dt*1e4*(fempos[0] - rbpos[0]);
			f[1] = dt*1e4*(fempos[1] - rbpos[1]);
			f[2] = dt*1e4*(fempos[2] - rbpos[2]);
			t = f.cross(rb.bulletBody->getCenterOfMassPosition()-rbpos);
			r[femIndex  ] -= f[0]; r[femIndex+1] -= f[1]; r[femIndex+2] -= f[2];
			r[rbIndex  ] += f[0];  r[rbIndex+1] += f[1];  r[rbIndex+2] += f[2];
			r[rbIndex+3] += t[0];  r[rbIndex+4] += t[1];  r[rbIndex+5] += t[2];
		}
		rbIndex+=6;
	}

	// apply preconditioner to the residual to get intial search direction, d = p*r
	dptr1 = d; dptr2 = p; dptr3 = r;
	for (unsigned int i=0; i<nfemdof; i++, dptr1++, dptr2++, dptr3++) {
		(*dptr1) = (*dptr2) * (*dptr3);
	}
	applyRigidPreconditioner(r+nfemdof, d+nfemdof);
	applySpringPreconditioner(r, d);

	delta_new = cblas_ddot(n, (double*)r, 1, (double*)d, 1);
	tol *= delta_new;

	if (delta_new > tol) {
		while (iter++ < max_iter) {
			// apply matrix
			offset = 0;
			for (unsigned int i=0; i<femObjects.size(); i++) {
				femObjects[i].gm->mvmul((SlVector3 *)(d+offset),(SlVector3*) (q+offset));
				offset += 3*femObjects[i].nv;																	
			}
			mulRigidMassMatrix(d+offset, q+offset);
			mulSpringMatrix(d,q);

			alpha = delta_new/cblas_ddot(n, d, 1, q, 1);
			cblas_daxpy(n, alpha, d, 1, x, 1);
			cblas_daxpy(n, -alpha, q, 1, r, 1);
			
			// apply preconditioner q = P*r
			dptr1 = q; dptr3 = p; dptr2 = r;
			for (unsigned int i=0; i<nfemdof; i++, dptr1++, dptr2++, dptr3++) {
				(*dptr1) = (*dptr2) * (*dptr3);
			}
			applyRigidPreconditioner(r+nfemdof, q+nfemdof);
			applySpringPreconditioner(r, q);

			delta_old = delta_new;
			delta_new = cblas_ddot(n, r, 1, q, 1);
			beta = delta_new / delta_old;
			
			cblas_daxpy(n, beta, d, 1, q, 1);
			cblas_dcopy(n, q, 1, d, 1);

			//std::cout<<iter<<" "<<delta_new<<" "<<tol<<std::endl;
			if (delta_new < tol) break;
		}
	}

	dptr1 = x;
	for (unsigned int i=0; i<femObjects.size(); i++) {
		SlVector3 *vptr = femObjects[i].vel;
		for (unsigned int j=0; j<femObjects[i].nv; j++, vptr++) {
			(*vptr)[0] = (*(dptr1++));  (*vptr)[1] = (*(dptr1++));  (*vptr)[2] = (*(dptr1++));
		}
	}
	for (unsigned int i=0; i<rigidBodies.size(); i++) {
		if (rigidBodies[i].rbType == RigidBody::RBType::RB_PLANE) continue;
		rigidBodies[i].bulletBody->setLinearVelocity(btVector3(dptr1[0], dptr1[1], dptr1[2]));
		rigidBodies[i].bulletBody->setAngularVelocity(btVector3(dptr1[3], dptr1[4], dptr1[5]));
		dptr1 += 6;
	}
}


void World::computeFemVelocities(){
  
  for(auto& femObject : femObjects){
    if(femObject.active){
      femObject.setForces(dt, SlVector3(gravity[0], gravity[1], gravity[2]));
    } else {
      femObject.setForces(dt, 0.0);
    }
		for (unsigned int v=0; v<femObject.nv; v++) {
			femObject.separationTensor[v] = 0.0;
			femObject.unbalancedTensLoad[v] = 0.0;
			femObject.unbalancedCompLoad[v] = 0.0;
		}
  }

  for(auto& femObject : femObjects){
    femObject.computeForces(dt);
    femObject.clearGlobalMatrix();
    femObject.computeGlobalMatrix(dt);
    //femObject.solveForVelocities(dt);
  }

}


void World::updateFemPositions(){
  for(auto& femObject : femObjects){
    femObject.updatePositions(dt);
		femObject.fracture();
  }

}


void World::dumpFrame(){

  char framestring[80];
  sprintf(framestring, "frames/foo-%%03i.%%04i.obj");
  std::cout << "writing frame: " << currentFrame << std::endl;

  int objectCount = 0;
  char fname[80];
  for(auto& femObject : femObjects){
    sprintf(fname, framestring, objectCount, currentFrame);
    femObject.dumpObj(fname);
    objectCount++;
  }
  for(auto& rigidBody : rigidBodies){
    sprintf(fname, framestring, objectCount, currentFrame);
    rigidBody.dump(fname);
	objectCount++;
  }
  currentFrame++;

}

void World::computeConstraints(){
  std::cout << "computing constraints" << std::endl;
  //find the FEM nodes that are inside of a rigid body and constrain them.

  for( auto rbInd : range(rigidBodies.size())){
    auto& rb = rigidBodies[rbInd];
    auto mass = 1.0/rb.bulletBody->getInvMass();
    for(auto femInd : range(femObjects.size())){
      auto& fem = femObjects[femInd];

      //cull based on BB overlap
      btVector3 rbBoundMin, rbBoundMax;
      rb.bulletBody->getAabb(rbBoundMin, rbBoundMax);
      
      if(!(rbBoundMin.x() > fem.bbMax.x() ||
	   rbBoundMax.x() < fem.bbMin.x() ||
	   rbBoundMin.y() > fem.bbMax.y() ||
	   rbBoundMax.y() < fem.bbMin.y() ||
	   rbBoundMin.z() > fem.bbMax.z() ||
	   rbBoundMax.z() < fem.bbMin.z())){
	std::cout << "aabb overlap" << std::endl;
	for(auto i : range(fem.nv)){

	  auto worldPos = btVector3(fem.pos[i].x(),
				    fem.pos[i].y(),
				    fem.pos[i].z());
	  btVector3 relPos = rb.bulletBody->getWorldTransform().inverse()*worldPos;

	  //check if the relPos is in the object
	  
	  if(rb.rbType == RigidBody::RB_BOX){
	    const btBoxShape *box = dynamic_cast<btBoxShape*>(rb.bulletBody->getCollisionShape());
	    auto extents = box->getHalfExtentsWithMargin();
	    if( relPos.x() >= -extents.x()  &&
		relPos.x() <= extents.x() &&
		relPos.y() >= -extents.y() &&
		relPos.y() <= extents.y() &&
		relPos.z() >= -extents.z() &&
		relPos.z() <= extents.z()){
	      
	      rb.constraints.push_back({relPos, femInd, i, rbInd, Eigen::Matrix3d::Zero(), 1.0});//std::min(fem.mass[i], mass)-1e-6});
				fem.grip(femInd);
	      std::cout << "added constraint: " << relPos << ' ' << femInd << ' ' << i << std::endl;
	    }
	    

	  } else if(rb.rbType == RigidBody::RB_PLANE){
	    const auto* plane = dynamic_cast<btStaticPlaneShape*>(rb.bulletBody->getCollisionShape());
	    if(relPos.dot(plane->getPlaneNormal()) <= plane->getPlaneConstant()){
	      rb.constraints.push_back({relPos, femInd, i, rbInd, Eigen::Matrix3d::Zero(), std::min(fem.mass[i], mass)-1e-6});
				fem.grip(femInd);
	      std::cout << "added constraint: " << relPos << ' ' << femInd << ' ' << i << std::endl;
	    }
	    

	  }else {
	    std::cout << "unknown RB type" << std::endl;
	    exit(1);
	  }
	}
      }
    }
  }
}

void World::constraintCullingCallback(btBroadphasePair& collisionPair, 
									  btCollisionDispatcher& _dispatcher, 
									  const btDispatcherInfo& dispatcherInfo){
  //seems pretty slow, but here's a fist shot
  //there are probably fewer rigidbodies than tets, so lets's search there first
  
  auto* ptr0 = reinterpret_cast<btRigidBody*>(collisionPair.m_pProxy0->m_clientObject);
  auto* ptr1 = reinterpret_cast<btRigidBody*>(collisionPair.m_pProxy1->m_clientObject);


  auto predicate = [&](const RigidBody& rb){
	return (rb.bulletBody.get() == ptr0 ||
			rb.bulletBody.get() == ptr1);
  };

  auto it1 = std::find_if(rigidBodies.begin(), rigidBodies.end(),predicate);
  if(it1 != rigidBodies.end()){
	//std::cout << "found a RB" << std::endl;
	//one of them is a rigid body
	auto it2 = std::find_if(it1 + 1, rigidBodies.end(), predicate);
	if(it2 != rigidBodies.end()){
	  //they're both rigid bodies, so process collisions normally
	  dispatcher->defaultNearCallback(collisionPair, _dispatcher, dispatcherInfo);
	} else {
	  //see if they're joined by a constraint
	  auto& rb = *it1;
	  auto it3 = std::find_if(rb.constraints.begin(), rb.constraints.end(),
							  [&](const CouplingConstraint& constraint) -> bool{
								
								//auto* ptr = 
								//femObjects[constraint.femIndex].bulletTetBodies[constraint.nodeIndex].get();
								auto* ptr = 
								  femObjects[constraint.femIndex].particleRigidBodies[constraint.nodeIndex].get();
							 return(ptr == ptr0 ||
									ptr == ptr1);
						   });
	  if(it3 != rb.constraints.end()){
		//they're joined, so don't process the collision
		//std::cout << "not colliding constrained tet with rb" << std::endl;
	  } else {
		//there is no constraint, handle them naturally
		dispatcher->defaultNearCallback(collisionPair, _dispatcher, dispatcherInfo);
	  }
	}
  } else {
	//neither are rigid bodies, so process the collision normally
	dispatcher->defaultNearCallback(collisionPair, _dispatcher, dispatcherInfo);
	
  }
}



void World::computeCrossProductMatrices(){

  for(auto& rb : rigidBodies){
    for(auto& c : rb.constraints){
      c.crossProductMatrix = 
	crossProductMatrix(quatRotate(rb.bulletBody->getOrientation(),c.localPosition));
    }
  }
}


void World::applySpringPreconditioner(double *in, double *out) {
	// do nothing for now...
}

void World::mulSpringMatrix(double* in, double*out){
	size_t rbIndex = nfemdof;
	double mult = dt*dt*1e4;

  for(auto& rb : rigidBodies){
		if (rb.rbType == RigidBody::RBType::RB_PLANE) continue;
    for(auto& c : rb.constraints){
      size_t femIndex = 3*(femObjects[c.femIndex].firstNodeIndex + c.nodeIndex);

			out[femIndex    ] += mult*in[femIndex    ];
			out[femIndex + 1] += mult*in[femIndex + 1];
			out[femIndex + 2] += mult*in[femIndex + 2];
			out[femIndex    ] -= mult*in[rbIndex    ];
			out[femIndex + 1] -= mult*in[rbIndex + 1];
			out[femIndex + 2] -= mult*in[rbIndex + 2];

			out[rbIndex    ] += mult*in[rbIndex    ];
			out[rbIndex + 1] += mult*in[rbIndex + 1];
			out[rbIndex + 2] += mult*in[rbIndex + 2];
			out[rbIndex    ] -= mult*in[femIndex    ];
			out[rbIndex + 1] -= mult*in[femIndex + 1];
			out[rbIndex + 2] -= mult*in[femIndex + 2];

			double v[3] = {mult*in[rbIndex+3],mult*in[rbIndex+4],mult*in[rbIndex+5]};

			inPlaceMult(v, out+femIndex, c.crossProductMatrix);
			inPlaceMultTranspose(v, out+rbIndex, c.crossProductMatrix);
    }
		rbIndex += 6;
  }
  
}

void World::mulJV(double* in, double*out){
	size_t rbIndex = nfemdof;
  size_t constraintIndex = nrbdof + nfemdof;
	double v[3];

  for(auto& rb : rigidBodies){
		if (rb.rbType == RigidBody::RBType::RB_PLANE) continue;
    for(auto& c : rb.constraints){
      size_t femIndex = 3*(femObjects[c.femIndex].firstNodeIndex + c.nodeIndex);
      //the constaint block is I_3x3 for the FEM,
      // -I_3x3 for the rigid linear
      // and a cross product matrix the rotational part
      out[constraintIndex    ] += c.weight*(in[femIndex   ] - in[rbIndex   ]);
      out[constraintIndex + 1] += c.weight*(in[femIndex +1] - in[rbIndex +1]);
      out[constraintIndex + 2] += c.weight*(in[femIndex +2] - in[rbIndex +2]);
      
			v[0] = c.weight*in[rbIndex + 3];
			v[1] = c.weight*in[rbIndex + 4];
			v[2] = c.weight*in[rbIndex + 5];

      inPlaceMult(v,  
								out + constraintIndex, 
								c.crossProductMatrix);


      constraintIndex+=3;
    }
		rbIndex += 6;
  }
  
}

void World::mulJTLambda(double* in, double* out){
	size_t rbIndex = nfemdof;
  size_t constraintIndex = nrbdof+nfemdof;
	double v[3];

  for(auto& rb : rigidBodies){
		if (rb.rbType == RigidBody::RBType::RB_PLANE) continue;
    for(auto& c : rb.constraints){
      size_t femIndex = 3*(femObjects[c.femIndex].firstNodeIndex + c.nodeIndex);
      
			v[0] = c.weight*in[constraintIndex    ];
			v[1] = c.weight*in[constraintIndex + 1];
			v[2] = c.weight*in[constraintIndex + 2];

      out[femIndex    ] += v[0];
      out[femIndex + 1] += v[1];
      out[femIndex + 2] += v[2];
      
      out[rbIndex    ] -= v[0];
      out[rbIndex + 1] -= v[1];
      out[rbIndex + 2] -= v[2];

      inPlaceMultTranspose(v,
			 out + rbIndex + 3,
			 c.crossProductMatrix);

      constraintIndex+=3;
    }
		rbIndex+=6;
  }
}


//make sure in and out are the right size
//in = first rigid body velocity in teh vector
//out = first rigid body velocity in the vector
void World::mulRigidMassMatrix(double* in, double* out){
	int index = 0;
  for(auto pr : enumerate(rigidBodies)){
		if (pr.element.rbType == RigidBody::RBType::RB_PLANE) continue;
    auto mass = 1.0/pr.element.bulletBody->getInvMass();
    out[index    ] = mass*in[index    ];
    out[index + 1] = mass*in[index + 1];
    out[index + 2] = mass*in[index + 2];
    //rotational part
    bulletOverwriteMultiply(in + index + 3, 
			    out + index + 3,
			    pr.element.bulletBody->getInvInertiaTensorWorld().inverse());

		index+=6;
  }
}

void World::applyRigidPreconditioner(double* in, double* out){
  
	int index = 0;
  for(auto pr : enumerate(rigidBodies)){
    auto& rb = pr.element;
		if (rb.rbType == RigidBody::RBType::RB_PLANE) continue;
    out[index    ] = rb.bulletBody->getInvMass()*in[index    ];
    out[index + 1] = rb.bulletBody->getInvMass()*in[index + 1];
    out[index + 2] = rb.bulletBody->getInvMass()*in[index + 2];
    
    bulletOverwriteMultiply(in + index + 3,
			    out + index + 3,
			    rb.bulletBody->getInvInertiaTensorWorld());
		index+=6;
  }
}

void World::countConstraints(){
  totalNumConstraints = std::accumulate(rigidBodies.begin(),
					rigidBodies.end(),
					0,
					[&](size_t count, const RigidBody& rb){
					  return count + rb.constraints.size();
					});
}

void World::applyRegularizer(double* in, double* out){
  size_t constraintIndex = 0;

  for(auto& rb : rigidBodies){
		if (rb.rbType == RigidBody::RBType::RB_PLANE) continue;
    for(auto& c : rb.constraints){
			out[constraintIndex    ] = 0.0;//c.weight*in[constraintIndex    ];
			out[constraintIndex + 1] = 0.0;//c.weight*in[constraintIndex + 1];
			out[constraintIndex + 2] = 0.0;//c.weight*in[constraintIndex + 2];
			constraintIndex += 3;
		}
	}  
}

void World::applyRegularizerPreconditioner(double* in, double* out){
  size_t constraintIndex = 0;

  for(auto& rb : rigidBodies){
		if (rb.rbType == RigidBody::RBType::RB_PLANE) continue;
    for(auto& c : rb.constraints){
			out[constraintIndex    ] = 0.0;//in[constraintIndex    ] / (c.weight);
			out[constraintIndex + 1] = 0.0;//in[constraintIndex + 1] / (c.weight);
			out[constraintIndex + 2] = 0.0;//in[constraintIndex + 2] / (c.weight);
			constraintIndex += 3;
		}
	}  
}


bool output = false;
void World::project(double *in, RigidBody &rb, int rbIndex, CouplingConstraint &c) {
	Eigen::Vector3d v, CTv;
	size_t femIndex = 3*(femObjects[c.femIndex].firstNodeIndex + c.nodeIndex);
	double rbInvMass = rb.bulletBody->getInvMass();
	double femInvMass = 1.0/femObjects[c.femIndex].mass[c.nodeIndex];
	double sumInvMass = rbInvMass + femInvMass;

	// project velocity onto the constraint (C^T v)

	//CTv[0] = femInvMass*in[femIndex   ] - rbInvMass*in[rbIndex   ];
	//CTv[1] = femInvMass*in[femIndex +1] - rbInvMass*in[rbIndex +1];
	//CTv[2] = femInvMass*in[femIndex +2] - rbInvMass*in[rbIndex +2];
	CTv[0] = in[femIndex   ] - in[rbIndex   ];
	CTv[1] = in[femIndex +1] - in[rbIndex +1];
	CTv[2] = in[femIndex +2] - in[rbIndex +2];
  
	v[0] = in[rbIndex + 3];
	v[1] = in[rbIndex + 4];
	v[2] = in[rbIndex + 5];
	
	//CTv += c.crossProductMatrix * v;
	
	// compute the mass weighting matrix (C^T M^-1 C)^-1
	const btMatrix3x3 &btInvI = rb.bulletBody->getInvInertiaTensorWorld();
	Eigen::Matrix3d Iinv;
	Iinv << 
		btInvI[0][0], btInvI[0][1], btInvI[0][2],
		btInvI[1][0], btInvI[1][1], btInvI[1][2],
		btInvI[2][0], btInvI[2][1], btInvI[2][2];
	Eigen::Matrix3d CTMinvC; //= c.crossProductMatrix * Iinv * c.crossProductMatrix.transpose();
	CTMinvC << 0,0,0, 0,0,0, 0,0,0;
	CTMinvC(0,0) += sumInvMass; 
	CTMinvC(1,1) += sumInvMass; 
	CTMinvC(2,2) += sumInvMass; 
	
	Eigen::Vector3d CTMinvCinvCTv = CTMinvC.inverse()*CTv;
	

	Eigen::Vector3d CCTMinvCinvCTv[3];
	CCTMinvCinvCTv[0] = CTMinvCinvCTv;
	CCTMinvCinvCTv[1] = -CTMinvCinvCTv;
	CCTMinvCinvCTv[2] = c.crossProductMatrix.transpose()*CTMinvCinvCTv;
	
	CCTMinvCinvCTv[0] *= femInvMass;
	CCTMinvCinvCTv[1] *= rbInvMass;
	CCTMinvCinvCTv[2] = Iinv * CCTMinvCinvCTv[2];
	

	in[femIndex    ] -= CCTMinvCinvCTv[0][0];
	in[femIndex  +1] -= CCTMinvCinvCTv[0][1];
	in[femIndex  +2] -= CCTMinvCinvCTv[0][2];
	in[rbIndex    ] -= CCTMinvCinvCTv[1][0];
	in[rbIndex  +1] -= CCTMinvCinvCTv[1][1];
	in[rbIndex  +2] -= CCTMinvCinvCTv[1][2];
}

void World::project(double *in) {
	for (unsigned int i=0; i<4; i++) {

		size_t rbIndex = nfemdof;
		for(auto& rb : rigidBodies){
			if (rb.rbType == RigidBody::RBType::RB_PLANE) continue;
			for(auto& c : rb.constraints){
				project(in, rb, rbIndex, c);
			}
			rbIndex += 6;
		}

		rbIndex = nfemdof+nrbdof-6;
		for (int j = rigidBodies.size()-1; j >= 0; j--) {
			auto &rb = rigidBodies[j];
			if (rb.rbType == RigidBody::RBType::RB_PLANE) continue;
			for (int k = rb.constraints.size()-1; k >= 0; k--) {
				if (k==3) output = true;
				auto &c = rb.constraints[k];
				project(in, rb, rbIndex, c);
				output=false;
			}
			rbIndex -= 6;
		}
	}
}

class Preconditioner {
public:
	Preconditioner(World *w) {
		this->w = w;
		p = new double [w->nfemdof];
		int offset = 0;
		for (unsigned int i=0; i<w->femObjects.size(); i++) {
			w->femObjects[i].gm->setPrecon(p+offset);
			offset += 3*w->femObjects[i].nv;
		}
	}

	~Preconditioner() {
		delete [] p;
	}

	void Apply(const SimpleVector &x, SimpleVector &y) const {
		double *dptr1 = y.vals; double *dptr2 = p; double *dptr3 = x.vals;
		for (unsigned int i=0; i<w->nfemdof; i++, dptr1++, dptr2++, dptr3++) {
			(*dptr1) = (*dptr2) * (*dptr3);
		}
		w->applyRigidPreconditioner(x.vals+w->nfemdof, y.vals+w->nfemdof);
		w->applyRegularizerPreconditioner(x.vals+w->nfemdof+w->nrbdof, y.vals+w->nfemdof+w->nrbdof);
	}

private:
	World *w;
	double *p;
};
	
class Operator {
public:
	Operator(World *w) {
		this->w = w;
	}

	~Operator() {
	}
	
	void Apply(const SimpleVector &x, SimpleVector &y) const {
		int offset = 0;
		for (unsigned int i=0; i<w->femObjects.size(); i++) {
			w->femObjects[i].gm->mvmul((SlVector3 *)(x.vals+offset),(SlVector3*) (y.vals+offset));
			offset += 3*w->femObjects[i].nv;																	
		}
		w->mulRigidMassMatrix(x.vals+offset, y.vals+offset);
		w->applyRegularizer(x.vals+offset+w->nrbdof, y.vals+offset+w->nrbdof);
		w->mulJV(x.vals, y.vals);
		w->mulJTLambda(x.vals, y.vals);
	}
private:
	World *w;
};

void World::solveMinres() {

	nfemdof = 0;
	nrbdof = 0;
	unsigned int nrbs = 0;
	unsigned int n = 0; 

	for (unsigned int i=0; i<femObjects.size(); i++) {
			femObjects[i].firstNodeIndex = n/3;
			n += 3*femObjects[i].nv;
			nfemdof += 3*femObjects[i].nv;
		}
	
	for (unsigned int i=0; i<rigidBodies.size(); i++) {
		if (rigidBodies[i].rbType == RigidBody::RBType::RB_PLANE) continue;
		n += 6;
		nrbs++;
		nrbdof += 6;
	}
	n += 3*totalNumConstraints;	

	SimpleVector rhs(n), x(n);
	x = 0;

	double *dptr1 = rhs.vals;
	for (unsigned int i=0; i<femObjects.size(); i++) {
		SlVector3 *vptr = femObjects[i].frc;
		for (unsigned int j=0; j<femObjects[i].nv; j++, vptr++) {
			(*(dptr1++)) = (*vptr)[0]; (*(dptr1++)) = (*vptr)[1]; (*(dptr1++)) = (*vptr)[2];
		}
	}

	for (unsigned int i=0; i<rigidBodies.size(); i++) {
		if (rigidBodies[i].rbType == RigidBody::RBType::RB_PLANE) continue;
		double mass = 1.0 / rigidBodies[i].bulletBody->getInvMass();
    (*(dptr1++)) = mass*rigidBodies[i].bulletBody->getLinearVelocity()[0];
    (*(dptr1++)) = mass*rigidBodies[i].bulletBody->getLinearVelocity()[1];
    (*(dptr1++)) = mass*rigidBodies[i].bulletBody->getLinearVelocity()[2];

    //rotational part, there must be a less hideous way of doing this.
		double angularvelocity[3];
		angularvelocity[0] = rigidBodies[i].bulletBody->getAngularVelocity()[0];
		angularvelocity[1] = rigidBodies[i].bulletBody->getAngularVelocity()[1];
		angularvelocity[2] = rigidBodies[i].bulletBody->getAngularVelocity()[2];
    bulletOverwriteMultiply(angularvelocity, dptr1,
			    rigidBodies[i].bulletBody->getInvInertiaTensorWorld().inverse());
		dptr1 += 3;
	}
	for (unsigned int i=0; i<3*totalNumConstraints; i++) {
		(*(dptr1++)) = 0.0;
	}

	Operator op(this);
	//Preconditioner prec(this);
	Preconditioner *prec = NULL;

	double shift(0);
	int max_iter(100);
	double tol(1e-12);
	bool show(true);

	MINRES(op, x, rhs, prec, shift, max_iter, tol, show);

	dptr1 = x.vals;
	for (unsigned int i=0; i<femObjects.size(); i++) {
		SlVector3 *vptr = femObjects[i].vel;
		for (unsigned int j=0; j<femObjects[i].nv; j++, vptr++) {
			(*vptr)[0] = (*(dptr1++));
			(*vptr)[1] = (*(dptr1++));
			(*vptr)[2] = (*(dptr1++));
		}
	}
	for (unsigned int i=0; i<rigidBodies.size(); i++) {
		if (rigidBodies[i].rbType == RigidBody::RBType::RB_PLANE) continue;
		rigidBodies[i].bulletBody->setLinearVelocity(btVector3(dptr1[0], dptr1[1], dptr1[2]));
		rigidBodies[i].bulletBody->setAngularVelocity(btVector3(dptr1[3], dptr1[4], dptr1[5]));
		dptr1 += 6;
	}
}

#if 0
void World::solve() {
	double delta_new, delta_old, alpha, beta, *dptr1, *dptr2, *dptr3;
	double tol=1e-40;
	unsigned int iter=0, max_iter = 10000;
	nfemdof = 0;
	nrbdof = 0;
	unsigned int nrbs = 0;

	unsigned int n = 0; 
	for (unsigned int i=0; i<femObjects.size(); i++) {
		femObjects[i].firstNodeIndex = n/3;
		n += 3*femObjects[i].nv;
		nfemdof += 3*femObjects[i].nv;
	}
	for (unsigned int i=0; i<rigidBodies.size(); i++) {
		if (rigidBodies[i].rbType == RigidBody::RBType::RB_PLANE) continue;
		n += 6;
		nrbs++;
		nrbdof += 6;
	}
	n += 3*totalNumConstraints;	

	double *work = new double [5*n + nfemdof];
	double *x = work;     // solution
	double *r = work+n;   // residual
	double *q = work+2*n; // temp storage
	double *d = work+3*n; // search direction
	double *q2 = work+3*n; // search direction
	double *p = work+4*n; // preconditioner
	
	for (unsigned int i=0; i<5*n+nfemdof; i++) work[i] = 0.0;

	int offset = 0;
	for (unsigned int i=0; i<femObjects.size(); i++) {
		femObjects[i].gm->setPrecon(p+offset);
		offset += 3*femObjects[i].nv;
	}

	dptr1 = r;
	dptr2 = x;
	for (unsigned int i=0; i<n; i++, dptr2++) (*dptr2) = 0.0;

	for (unsigned int i=0; i<femObjects.size(); i++) {
		SlVector3 *vptr = femObjects[i].frc;
		for (unsigned int j=0; j<femObjects[i].nv; j++, vptr++) {
			(*(dptr1++)) = (*vptr)[0]; (*(dptr1++)) = (*vptr)[1]; (*(dptr1++)) = (*vptr)[2];
			//(*(dptr2++)) = (*vptr)[0] / femObjects[i].mass[j];
			//(*(dptr2++)) = (*vptr)[1] / femObjects[i].mass[j];
			//(*(dptr2++)) = (*vptr)[2] / femObjects[i].mass[j];
			//std::cout<<(*vptr)[1]/femObjects[i].mass[j]<<" "<<(*(dptr2++))<<" "<<femObjects[i].mass[j]<<std::endl;
		}
	}

	for (unsigned int i=0; i<rigidBodies.size(); i++) {
		if (rigidBodies[i].rbType == RigidBody::RBType::RB_PLANE) continue;
		double mass = 1.0 / rigidBodies[i].bulletBody->getInvMass();
    (*(dptr1++)) = mass*rigidBodies[i].bulletBody->getLinearVelocity()[0];
    (*(dptr1++)) = mass*rigidBodies[i].bulletBody->getLinearVelocity()[1];
    (*(dptr1++)) = mass*rigidBodies[i].bulletBody->getLinearVelocity()[2];
		//(*(dptr2++)) = rigidBodies[i].bulletBody->getLinearVelocity()[0];
		//(*(dptr2++)) = rigidBodies[i].bulletBody->getLinearVelocity()[1];
		//(*(dptr2++)) = rigidBodies[i].bulletBody->getLinearVelocity()[2];

    //rotational part, there must be a less hideous way of doing this.
		double angularvelocity[3];
		angularvelocity[0] = rigidBodies[i].bulletBody->getAngularVelocity()[0];
		angularvelocity[1] = rigidBodies[i].bulletBody->getAngularVelocity()[1];
		angularvelocity[2] = rigidBodies[i].bulletBody->getAngularVelocity()[2];
    bulletOverwriteMultiply(angularvelocity, dptr1,
			    rigidBodies[i].bulletBody->getInvInertiaTensorWorld().inverse());
		dptr1 += 3;
		//(*(dptr2++)) = angularvelocity[0];
		//(*(dptr2++)) = angularvelocity[1];
		//(*(dptr2++)) = angularvelocity[2];
	}
	for (unsigned int i=0; i<3*totalNumConstraints; i++) {
		(*(dptr1++)) = 0.0;
		//(*(dptr2++)) = 0.0;
	}

	//offset = 0;
	//for (unsigned int i=0; i<femObjects.size(); i++) {
	//	femObjects[i].gm->mvmul((SlVector3 *)(x+offset),(SlVector3*) (q+offset));
	//	offset += 3*femObjects[i].nv;																	
	//}
	//mulRigidMassMatrix(x+offset, q+offset);
	//applyRegularizer(d+offset+nrbdof, q+offset+nrbdof);
	//mulJV(d, q);
	//mulJTLambda(d, q);
	
	//cblas_daxpy(n, -1, q, 1, r, 1); 
	
	//std::cout<<"project r"<<std::endl;
	//project(r);

	// d = p*r
	dptr1 = d; dptr2 = p; dptr3 = r;
	for (unsigned int i=0; i<nfemdof; i++, dptr1++, dptr2++, dptr3++) {
		(*dptr1) = (*dptr2) * (*dptr3);
	}
	applyRigidPreconditioner(r+nfemdof, d+nfemdof);
	//applyRegularizerPreconditioner(r+nfemdof+nrbdof, d+nfemdof+nrbdof);

	//project(r);
	//project(d);
	delta_new = cblas_ddot(n, (double*)r, 1, (double*)d, 1);
	tol *= delta_new;

	if (delta_new > tol) {
		while (iter++ < max_iter) {
			// apply matrix
			project(d);
			offset = 0;
			for (unsigned int i=0; i<femObjects.size(); i++) {
				femObjects[i].gm->mvmul((SlVector3 *)(d+offset),(SlVector3*) (q+offset));
				offset += 3*femObjects[i].nv;																	
			}
			mulRigidMassMatrix(d+offset, q+offset);
			//applyRegularizer(d+offset+nrbdof, q+offset+nrbdof);
			//mulJV(d, q);
			//mulJTLambda(d, q);
			project(q);

			alpha = delta_new/cblas_ddot(n, d, 1, q, 1);
			cblas_daxpy(n, alpha, d, 1, x, 1);
			cblas_daxpy(n, -alpha, q, 1, r, 1);
			
			// apply preconditioner q = P*r
			dptr1 = q; dptr3 = p; dptr2 = r;
			for (unsigned int i=0; i<nfemdof; i++, dptr1++, dptr2++, dptr3++) {
				(*dptr1) = (*dptr2) * (*dptr3);
			}
			applyRigidPreconditioner(r+nfemdof, q+nfemdof);
			//applyRegularizerPreconditioner(r+nfemdof+nrbdof, q+nfemdof+nrbdof);

			delta_old = delta_new;
			delta_new = cblas_ddot(n, r, 1, q, 1);
			//project(q);
			beta = delta_new / delta_old;
			
			cblas_daxpy(n, beta, d, 1, q, 1);
			cblas_dcopy(n, q, 1, d, 1);
			//project(d);
			std::cout<<iter<<" "<<delta_new<<" "<<tol<<std::endl;
			if (delta_new < tol) break;// || delta_new > delta_old) break;
			
		}
	}

	dptr1 = x;
	for (unsigned int i=0; i<femObjects.size(); i++) {
		SlVector3 *vptr = femObjects[i].vel;
		for (unsigned int j=0; j<femObjects[i].nv; j++, vptr++) {
			(*vptr)[0] = (*(dptr1++));// / femObjects[i].mass[j]; 
			(*vptr)[1] = (*(dptr1++));// / femObjects[i].mass[j]; 
			(*vptr)[2] = (*(dptr1++));// / femObjects[i].mass[j];
		}
	}
	for (unsigned int i=0; i<rigidBodies.size(); i++) {
		if (rigidBodies[i].rbType == RigidBody::RBType::RB_PLANE) continue;
		rigidBodies[i].bulletBody->setLinearVelocity(btVector3(dptr1[0], dptr1[1], dptr1[2]));
		rigidBodies[i].bulletBody->setAngularVelocity(btVector3(dptr1[3], dptr1[4], dptr1[5]));
		dptr1 += 6;
	}
}
#endif

