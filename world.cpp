
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

void World::project(double *in, RigidBody &rb, int rbIndex, CouplingConstraint &c) {
	Eigen::Vector3d v, CTv;
	size_t femIndex = 3*(femObjects[c.femIndex].firstNodeIndex + c.nodeIndex);
	// project velocity onto the constraint (C^T v)

	CTv[0] = in[femIndex   ] - in[rbIndex   ];
	CTv[1] = in[femIndex +1] - in[rbIndex +1];
	CTv[2] = in[femIndex +2] - in[rbIndex +2];
  
	v[0] = in[rbIndex + 3];
	v[1] = in[rbIndex + 4];
	v[2] = in[rbIndex + 5];
	
	CTv += c.crossProductMatrix * v;
	
	// compute the mass weighting matrix (C^T M^-1 C)^-1
	const btMatrix3x3 &btInvI = rb.bulletBody->getInvInertiaTensorWorld();
	Eigen::Matrix3d Iinv;
	Iinv << 
		btInvI[0][0], btInvI[0][1], btInvI[0][2],
		btInvI[1][0], btInvI[1][1], btInvI[1][2],
		btInvI[2][0], btInvI[2][1], btInvI[2][2];
	Eigen::Matrix3d CTMinvC = c.crossProductMatrix * Iinv * c.crossProductMatrix.transpose();
	double rbInvMass = rb.bulletBody->getInvMass();
	double femInvMass = 1.0/femObjects[c.femIndex].mass[c.nodeIndex];
	double sumInvMass = rbInvMass + femInvMass;
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
	
	//std::cout<<in[femIndex]<<" "<<CCTMinvCinvCTv[0][0]<<std::endl;
	std::cout<<in[femIndex + 1]<<" "<<CCTMinvCinvCTv[0][1]<<std::endl;
	//std::cout<<in[femIndex + 2]<<" "<<CCTMinvCinvCTv[0][2]<<std::endl;
	//std::cout<<in[rbIndex + 0]<<" "<<CCTMinvCinvCTv[1][0]<<std::endl;
	std::cout<<in[rbIndex + 1]<<" "<<CCTMinvCinvCTv[1][1]<<std::endl;
	//std::cout<<in[rbIndex + 2]<<" "<<CCTMinvCinvCTv[1][2]<<std::endl;
	std::cout<<in[rbIndex + 3]<<" "<<CCTMinvCinvCTv[2][0]<<std::endl;
	std::cout<<in[rbIndex + 4]<<" "<<CCTMinvCinvCTv[2][1]<<std::endl;
	std::cout<<in[rbIndex + 5]<<" "<<CCTMinvCinvCTv[2][2]<<std::endl;
	std::cout<<"------------------"<<std::endl;
	in[femIndex    ] -= CCTMinvCinvCTv[0][0];
	in[femIndex  +1] -= CCTMinvCinvCTv[0][1];
	in[femIndex  +2] -= CCTMinvCinvCTv[0][2];
	in[rbIndex    ] -= CCTMinvCinvCTv[1][0];
	in[rbIndex  +1] -= CCTMinvCinvCTv[1][1];
	in[rbIndex  +2] -= CCTMinvCinvCTv[1][2];
	in[rbIndex  +3] -= CCTMinvCinvCTv[2][0];
	in[rbIndex  +4] -= CCTMinvCinvCTv[2][1];
	in[rbIndex  +5] -= CCTMinvCinvCTv[2][2];
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
				auto &c = rb.constraints[k];
				project(in, rb, rbIndex, c);
			}
			rbIndex -= 6;
		}
	}
}


void World::solve() {
	double delta_new, delta_old, alpha, beta, *dptr1, *dptr2, *dptr3;
	double tol=1e-100;
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

	dptr1 = x;
	for (unsigned int i=0; i<n; i++, dptr1++) (*dptr1) = 0.0;

	dptr1 = r;
	for (unsigned int i=0; i<femObjects.size(); i++) {
		SlVector3 *vptr = femObjects[i].frc;
		for (unsigned int j=0; j<femObjects[i].nv; j++, vptr++) {
			(*(dptr1++)) = (*vptr)[0]; (*(dptr1++)) = (*vptr)[1]; (*(dptr1++)) = (*vptr)[2];
			//std::cout<<(*vptr)[1]/femObjects[i].mass[j]<<" "<<(*(dptr2++))<<" "<<femObjects[i].mass[j]<<std::endl;
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
	for (unsigned int i=0; i<3*totalNumConstraints; i++) (*(dptr1++)) = 0.0;
	
	//project(r);
	
	// d = p*r
	dptr1 = d; dptr2 = p; dptr3 = r;
	for (unsigned int i=0; i<nfemdof; i++, dptr1++, dptr2++, dptr3++) {
		(*dptr1) = (*dptr2) * (*dptr3);
		std::cout<<(*dptr1)<<std::endl;
	}
	applyRigidPreconditioner(r+nfemdof, d+nfemdof);
	//applyRegularizerPreconditioner(r+nfemdof+nrbdof, d+nfemdof+nrbdof);

	project(d);

	delta_new = cblas_ddot(n, (double*)r, 1, (double*)d, 1);
	tol *= delta_new;

	/*
	for (unsigned int i=0; i<0; i++) {
		for (unsigned int j=0; j<n; j++) {
			d[j] = 0.0;
		}
		d[i] = 1.0;
		offset = 0;
		for (unsigned int i=0; i<femObjects.size(); i++) {
			femObjects[i].gm->mvmul((SlVector3 *)(d+offset),(SlVector3*) (q+offset));
			offset += 3*femObjects[i].nv;																	
		}
		mulRigidMassMatrix(d+offset, q+offset);
		applyRegularizer(d+offset+nrbdof, q+offset+nrbdof);
		mulJV(d, q);
		mulJTLambda(d, q);
		double foo = q[i];
		std::cout<<i<<": "<<q[i]<<" ";

		for (unsigned int j=0; j<n; j++) {
			d[j] = 1.0;
		}
		d[i] = 0.0;
		offset = 0;
		for (unsigned int i=0; i<femObjects.size(); i++) {
			femObjects[i].gm->mvmul((SlVector3 *)(d+offset),(SlVector3*) (q+offset));
			offset += 3*femObjects[i].nv;																	
		}
		mulRigidMassMatrix(d+offset, q+offset);
		applyRegularizer(d+offset+nrbdof, q+offset+nrbdof);
		mulJV(d, q);
		mulJTLambda(d, q);
		std::cout<<q[i]<<std::endl;
		if(std::abs(q[i]) > foo) {
			std::cout<<"crap!!!"<<std::endl;
		}
		}
	//exit(0);
	*/

	if (cblas_dnrm2(n, (double*)r, 1) > tol) {
		while (iter++ < max_iter) {
			// apply matrix
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
			
			if (cblas_dnrm2(n, r, 1) < tol) break;
			
			// apply preconditioner q = P*r
			dptr1 = q; dptr3 = p; dptr2 = r;
			for (unsigned int i=0; i<nfemdof; i++, dptr1++, dptr2++, dptr3++) {
				(*dptr1) = (*dptr2) * (*dptr3);
			}
			applyRigidPreconditioner(r+nfemdof, q+nfemdof);
			//applyRegularizerPreconditioner(r+nfemdof+nrbdof, q+nfemdof+nrbdof);

			delta_old = delta_new;
			delta_new = cblas_ddot(n, r, 1, q, 1);
			beta = delta_new / delta_old;
			
			cblas_daxpy(n, beta, d, 1, q, 1);
			cblas_dcopy(n, q, 1, d, 1);
			project(d);
		}
	}

	/*	
		offset = 0;
		for (unsigned int i=0; i<femObjects.size(); i++) {
			femObjects[i].gm->mvmul((SlVector3 *)(x+offset),(SlVector3*) (q+offset));
			offset += 3*femObjects[i].nv;																	
		}
		mulRigidMassMatrix(x+offset, q+offset);
		applyRegularizer(x+offset+nrbdof, q+offset+nrbdof);
		mulJV(x, q);
		mulJTLambda(x, q);

		//std::cout<<iter<<" "<<cblas_dnrm2(n, r, 1)<<std::endl;
		//for (unsigned int i=nfemdof+nrbdof; i<n; i++) {
		//std::cout<<x[i]<<" "<<r[i]<<" "<<q[i]<<std::endl; 
		//}

	for (unsigned int k=0; k<0; k++) {
		double diag;
		double sum = 0;
		for (unsigned int i=0; i<n; i++) {
			for (unsigned int j=0; j<n; j++) {
				d[j] = 0.0;
			}
			d[i] = 1.0;
			offset = 0;
			for (unsigned int i=0; i<femObjects.size(); i++) {
				femObjects[i].gm->mvmul((SlVector3 *)(d+offset),(SlVector3*) (q+offset));
				offset += 3*femObjects[i].nv;																	
			}
			mulRigidMassMatrix(d+offset, q+offset);
			applyRegularizer(d+offset+nrbdof, q+offset+nrbdof);
			mulJV(d, q);
			mulJTLambda(d, q);
			if (i == k) diag = q[k];
			else sum += std::abs(q[k]);
			if (std::abs(q[k]) > 0) std::cout<<k<<", "<<i<<" = "<<q[k]<<"; ";
			//if (std::fabs(q[n-3]) > 0)
			//std::cout<<n-3<<", "<<i<<": "<<q[n-3]<<" ";
		}
		std::cout<<std::endl;
		std::cout<<diag<<" "<<sum<<std::endl;
		if (sum > diag) std::cout<<"crap!!"<<std::endl;
	}
	std::cout<<std::endl;
	for (unsigned int i=0; i<n; i++) {
		for (unsigned int j=0; j<n; j++) {
			d[j] = 0.0;
		}
		d[i] = 1.0;
		offset = 0;
		for (unsigned int i=0; i<femObjects.size(); i++) {
			femObjects[i].gm->mvmul((SlVector3 *)(d+offset),(SlVector3*) (q+offset));
			offset += 3*femObjects[i].nv;																	
		}
		mulRigidMassMatrix(d+offset, q+offset);
		applyRegularizer(d+offset+nrbdof, q+offset+nrbdof);
		mulJV(d, q);
		mulJTLambda(d, q);
		if (std::fabs(q[n-2]) > 0)
			std::cout<<n-2<<", "<<i<<": "<<q[n-2]<<" ";
	}
		std::cout<<std::endl;
	for (unsigned int i=0; i<n; i++) {
		for (unsigned int j=0; j<n; j++) {
			d[j] = 0.0;
		}
		d[i] = 1.0;
		offset = 0;
		for (unsigned int i=0; i<femObjects.size(); i++) {
			femObjects[i].gm->mvmul((SlVector3 *)(d+offset),(SlVector3*) (q+offset));
			offset += 3*femObjects[i].nv;																	
		}
		mulRigidMassMatrix(d+offset, q+offset);
		applyRegularizer(d+offset+nrbdof, q+offset+nrbdof);
		mulJV(d, q);
		mulJTLambda(d, q);
		if (std::fabs(q[n-1]) > 0)
			std::cout<<n-1<<", "<<i<<": "<<q[n-1]<<" ";
	}
	std::cout<<x[176]<<" "<<x[182]<<" "<<x[183]<<" "<<x[184]<<" "<<x[278]<<std::endl;
	//exit(-1);
	*/

	dptr1 = x;
	for (unsigned int i=0; i<femObjects.size(); i++) {
		SlVector3 *vptr = femObjects[i].vel;
		for (unsigned int j=0; j<femObjects[i].nv; j++, vptr++) {
			(*vptr)[0] = (*(dptr1++)); (*vptr)[1] = (*(dptr1++)); (*vptr)[2] = (*(dptr1++));
		}
	}
	for (unsigned int i=0; i<rigidBodies.size(); i++) {
		if (rigidBodies[i].rbType == RigidBody::RBType::RB_PLANE) continue;
		rigidBodies[i].bulletBody->setLinearVelocity(btVector3(dptr1[0], dptr1[1], dptr1[2]));
		rigidBodies[i].bulletBody->setAngularVelocity(btVector3(dptr1[3], dptr1[4], dptr1[5]));
		dptr1 += 6;
	}
}

void World::timeStep(){


  //replace these two with a single CG solve

  //make sure to call countConstraints every time the number of constraints changes
  //make sure to call computeCrossProductMatrices before calling the JV and 
  //JTLambda multiply functions
	computeCrossProductMatrices();

  computeFemVelocities(); // this is really just forces
  bulletWorld.stepSimulationVelocitiesOnly(dt); 

	solve();




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
	      
	      rb.constraints.push_back({relPos, femInd, i, rbInd, Eigen::Matrix3d::Zero(), std::min(fem.mass[i], mass)-1e-6});
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


void World::computeCrossProductMatrices(){

  for(auto& rb : rigidBodies){
    for(auto& c : rb.constraints){
      c.crossProductMatrix = 
	crossProductMatrix(quatRotate(rb.bulletBody->getOrientation(),c.localPosition));
    }
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
      //out[constraintIndex    ] += c.weight*(-in[femIndex   ] - in[rbIndex   ]);
      //out[constraintIndex + 1] += c.weight*(-in[femIndex +1] - in[rbIndex +1]);
      //out[constraintIndex + 2] += c.weight*(-in[femIndex +2] - in[rbIndex +2]);
      
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


//in = &(lambda[0])
//out = &(velocitiesToSolveFor[0])
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
			out[constraintIndex    ] = 4*c.weight*in[constraintIndex    ];
			out[constraintIndex + 1] = 4*c.weight*in[constraintIndex + 1];
			out[constraintIndex + 2] = 4*c.weight*in[constraintIndex + 2];
			constraintIndex += 3;
		}
	}  
  //for(auto i : range(3*totalNumConstraints)){
	//out[i] = regularizerAlpha*in[i];
  //}
}

void World::applyRegularizerPreconditioner(double* in, double* out){
  size_t constraintIndex = 0;

  for(auto& rb : rigidBodies){
		if (rb.rbType == RigidBody::RBType::RB_PLANE) continue;
    for(auto& c : rb.constraints){
			out[constraintIndex    ] = in[constraintIndex    ] / (4*c.weight);
			out[constraintIndex + 1] = in[constraintIndex + 1] / (4*c.weight);
			out[constraintIndex + 2] = in[constraintIndex + 2] / (4*c.weight);
			constraintIndex += 3;
		}
	}  

  //auto recip = 1.0/regularizerAlpha;
  //for(auto i : range(3*totalNumConstraints)){
	//out[i] = recip*in[i];
  //}
}
