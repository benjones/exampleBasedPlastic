
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

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#elif linux
extern "C" {
#include <cblas.h>
}
#endif

using iter::range;
using iter::enumerate;

World::World(std::string filename)
  :broadphaseInterface
   (std::unique_ptr<btBroadphaseInterface>(new btDbvtBroadphase())),
   collisionConfiguration
   (std::unique_ptr<btCollisionConfiguration>(new btDefaultCollisionConfiguration())),
   dispatcher
   (std::unique_ptr<btDispatcher>(new btCollisionDispatcher(collisionConfiguration.get()))),
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
    RB.shape = 
      std::unique_ptr<btCollisionShape>(new btStaticPlaneShape(btVector3(0,1,0), groundHeight));
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
  for(auto i : range(femObjectsIn.size())){
    auto& femObjectIn = femObjectsIn[i];
    auto fnameIn = femObjectIn["filename"];
    femObjects.emplace_back();
    femObjects.back().load(fnameIn.asString().c_str());
    
    if(!femObjectIn["offset"].isNull()){
      assert(femObjectIn["offset"].isArray() &&
	     femObjectIn["offset"].size() == 3);
      SlVector3 offset(femObjectIn["offset"][0].asDouble(),
		       femObjectIn["offset"][1].asDouble(),
		       femObjectIn["offset"][2].asDouble());
      auto& femObj = femObjects.back();
      femObj.bbMin += offset;
      femObj.bbMax += offset;
      for(auto  i : range(femObj.nv)){
	femObj.pos[i] += offset;
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
}

void World::solve() {
	double delta_new, delta_old, alpha, beta, *dptr1, *dptr2, *dptr3;
	double tol=1e-10;
	unsigned int iter=0, max_iter = 10000;
	unsigned int nfemdof = 0;

	unsigned int n = 6*rigidBodies.size()+totalNumConstraints;
	for (unsigned int i=0; i<femObjects.size(); i++) {
		n += 3*femObjects[i].nv;
		nfemdof += 3*femObjects[i].nv;
	}
	
	double *work = new double [4*n + nfemdof];
	double *x = work;     // solution
	double *r = work+n;   // residual
	double *q = work+2*n; // temp storage
	double *d = work+3*n; // search direction
	double *p = work+4*n; // preconditioner
	
	int offset = 0;
	for (unsigned int i=0; i<femObjects.size(); i++) {
		femObjects[i].gm->setPrecon(p+offset);
		offset += 3*femObjects[i].nv;
	}

	dptr1 = x;
	for (unsigned int i=0; i<n; i++, dptr1++) (*dptr1) = 0.0;

	dptr1 = r;
	for (unsigned int i=0; i<femObjects.size(); i++) {
		SlVector3 *vptr = femObjects[i].frc;
		for (unsigned int j=0; j<femObjects[i].nv; j++, vptr++) {
			(*(dptr1++)) = (*vptr)[0]; (*(dptr1++)) = (*vptr)[1]; (*(dptr1++)) = (*vptr)[2];
		}
	}
	
	// d = p*r
	dptr1 = d; dptr2 = p; dptr3 = r;
	for (unsigned int i=0; i<nfemdof; i++, dptr1++, dptr2++, dptr3++) {
		(*dptr1) = (*dptr2) * (*dptr3);
	}
	applyRigidPreconditioner(r+nfemdof, d+nfemdof);
	applyRegularizerPreconditioner(r+nfemdof+6*rigidBodies.size(), d+nfemdof+6*rigidBodies.size());


	delta_new = cblas_ddot(n, (double*)r, 1, (double*)d, 1);
	tol *= delta_new;
	
	if (cblas_dnrm2(n, (double*)r, 1) >= tol) {
		while (iter++ < max_iter) {
			// apply matrix
			offset = 0;
			for (unsigned int i=0; i<femObjects.size(); i++) {
				femObjects[i].gm->mvmul((SlVector3 *)(d+offset),(SlVector3*) (q+offset));
				offset += 3*femObjects[i].nv;																	
			}
			mulRigidMassMatrix(d+offset, q+offset);
			applyRegularizer(d+offset+6*rigidBodies.size(), q+offset+6*rigidBodies.size());
			mulJV(d, q);
			mulJTLambda(d, q);


			alpha = delta_new/cblas_ddot(n, d, 1, q, 1);
			cblas_daxpy(n, alpha, d, 1, x, 1);
			cblas_daxpy(n, -alpha, q, 1, r, 1);
			
			if (cblas_dnrm2(n, r, 1) < tol) break;
			
			// q = p*r
			dptr1 = q; dptr3 = p; dptr2 = r;
			for (unsigned int i=0; i<nfemdof; i++, dptr1++, dptr2++, dptr3++) {
				(*dptr1) = (*dptr2) * (*dptr3);
			}
			applyRigidPreconditioner(r+nfemdof, q+nfemdof);
			applyRegularizerPreconditioner(r+nfemdof+6*rigidBodies.size(), q+nfemdof+6*rigidBodies.size());

			delta_old = delta_new;
			delta_new = cblas_ddot(n, r, 1, q, 1);
			beta = delta_new / delta_old;
			
			cblas_daxpy(n, beta, d, 1, q, 1);
			cblas_dcopy(n, q, 1, d, 1);
		}
	}

	dptr1 = x;
	for (unsigned int i=0; i<femObjects.size(); i++) {
		SlVector3 *vptr = femObjects[i].vel;
		for (unsigned int j=0; j<femObjects[i].nv; j++, vptr++) {
			(*vptr)[0] = (*(dptr1++)); (*vptr)[1] = (*(dptr1++)); (*vptr)[2] = (*(dptr1++));
		}
	}
}

void World::timeStep(){


  //replace these two with a single CG solve

  //make sure to call countConstraints every time the number of constraints changes
  //make sure to call computeCrossProductMatrices before calling the JV and 
  //JTLambda multiply functions
  bulletWorld.stepSimulationVelocitiesOnly(dt); 

  computeFemVelocities();


  if(ground){
    for(auto& femObject : femObjects){
      for(auto j : range(femObject.nv)){
	//bullet uses y as up, I think, lets keep that convention
	auto x = femObject.pos[j][1] + dt*femObject.vel[j][1];
	if( x < 0){
	  auto impulse = -x/dt;
	  double nm= 0.0, m = sqrt(sqr(femObject.vel[j][0]) +
				   sqr(femObject.vel[j][2]));
	  if(m > 0){
	    nm = std::max(m - impulse*friction, 0.0) / m;
	  }
	  femObject.vel[j].set(nm*femObject.vel[j][0], 
			       femObject.vel[j][1] + impulse,
			       nm*femObject.vel[j][2]);
	}
      }
    }
  }
  
  /*
  for(auto &rb : rigidBodies){
    rb.couplingSolver.solveVelocityConstraints(*this, rb);
  }
  */

  bulletWorld.integrateTransforms(dt);
  bulletWorld.updateActivationState(dt);
  bulletWorld.synchronizeMotionStates();
  updateFemPositions();
  
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
    femObject.solveForVelocities(dt);
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
  }
  currentFrame++;

}

void World::computeConstraints(){
  std::cout << "computing constraints" << std::endl;
  //find the FEM nodes that are inside of a rigid body and constrain them.

  for( auto rbInd : range(rigidBodies.size())){
    auto& rb = rigidBodies[rbInd];
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
	      
	      rb.constraints.push_back({relPos, femInd, i, rbInd, Eigen::Matrix3d::Zero()});
	      std::cout << "added constraint: " << relPos << ' ' << femInd << ' ' << i << std::endl;
	      
	    }
	    

	  } else if(rb.rbType == RigidBody::RB_PLANE){
	    const auto* plane = dynamic_cast<btStaticPlaneShape*>(rb.bulletBody->getCollisionShape());
	    if(relPos.dot(plane->getPlaneNormal()) <= plane->getPlaneConstant()){
	      rb.constraints.push_back({relPos, femInd, i, rbInd, Eigen::Matrix3d::Zero()});
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


void World::computeCrossProductMatrices(){

  for(auto& rb : rigidBodies){
    for(auto& c : rb.constraints){
      c.crossProductMatrix = 
	crossProductMatrix(quatRotate(rb.bulletBody->getOrientation(),c.localPosition));
    }
  }
}

void World::mulJV(double* in, double*out){
  size_t constraintIndex = 0;
  size_t rbStart = femObjects.empty() ? 0 : 3*(femObjects.back().firstNodeIndex +
					       femObjects.back().nv);
  for(auto& rb : rigidBodies){
    for(auto& c : rb.constraints){
      size_t femIndex = femObjects[c.femIndex].firstNodeIndex + c.nodeIndex;
      size_t rbIndex = rbStart + 6*c.rbIndex;
      //the constaint block is I_3x3 for the FEM,
      // -I_3x3 for the rigid linear
      // and a cross product matrix the rotational part
      out[3*constraintIndex   ] = in[3*femIndex   ] - in[rbIndex];
      out[3*constraintIndex +1] = in[3*femIndex +1] - in[rbIndex +1];
      out[3*constraintIndex +2] = in[3*femIndex +2] - in[rbIndex +2];
      
      inPlaceMult(in + rbIndex + 3,  
		  out + 3*constraintIndex, 
		  c.crossProductMatrix);


      ++constraintIndex;
    }
  }
  
}


//in = &(lambda[0])
//out = &(velocitiesToSolveFor[0])
void World::mulJTLambda(double* in, double* out){
  size_t constraintIndex = 0;
  size_t rbStart = femObjects.empty() ? 0 : 3*(femObjects.back().firstNodeIndex +
					       femObjects.back().nv);
  
  for(auto& rb : rigidBodies){
    for(auto& c : rb.constraints){
      size_t femIndex = femObjects[c.femIndex].firstNodeIndex + c.nodeIndex;
      size_t rbIndex = rbStart + 6*c.rbIndex;
      
      out[3*femIndex    ] += in[constraintIndex*3    ];
      out[3*femIndex + 1] += in[constraintIndex*3 + 1];
      out[3*femIndex + 2] += in[constraintIndex*3 + 2];
      
      out[rbIndex    ] -= in[constraintIndex*3    ];
      out[rbIndex + 1] -= in[constraintIndex*3 + 1];
      out[rbIndex + 2] -= in[constraintIndex*3 + 2];

      inPlaceMultTranspose(in + constraintIndex*3,
			   out + rbIndex + 3,
			   c.crossProductMatrix);

      ++constraintIndex;
    }
  }
}

//make sure in and out are the right size
//in = first rigid body velocity in teh vector
//out = first rigid body velocity in the vector
void World::mulRigidMassMatrix(double* in, double* out){
  for(auto pr : enumerate(rigidBodies)){
    auto mass = 1.0/pr.element.bulletBody->getInvMass();
    out[6*pr.index    ] = mass*in[6*pr.index    ];
    out[6*pr.index + 1] = mass*in[6*pr.index + 1];
    out[6*pr.index + 2] = mass*in[6*pr.index + 2];
    //rotational part
    bulletOverwriteMultiply(in + 6*pr.index + 3, 
			    out + 6*pr.index + 3,
			    pr.element.bulletBody->getInvInertiaTensorWorld().inverse());

  }
}

void World::applyRigidPreconditioner(double* in, double* out){
  
  for(auto pr : enumerate(rigidBodies)){
    auto& rb = pr.element;
    out[6*pr.index    ] = rb.bulletBody->getInvMass()*in[6*pr.index    ];
    out[6*pr.index + 1] = rb.bulletBody->getInvMass()*in[6*pr.index + 1];
    out[6*pr.index + 2] = rb.bulletBody->getInvMass()*in[6*pr.index + 2];
    
    bulletOverwriteMultiply(in + 6*pr.index + 3,
			    out + 6*pr.index + 3,
			    rb.bulletBody->getInvInertiaTensorWorld());
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
  
  for(auto i : range(3*totalNumConstraints)){
    out[i] = regularizerAlpha*in[i];
  }
}

void World::applyRegularizerPreconditioner(double* in, double* out){

  auto recip = 1.0/regularizerAlpha;
  for(auto i : range(3*totalNumConstraints)){
    out[i] = recip*in[i];
  }
}
