
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

#include "tminres.hpp"
#include "SimpleVector.hpp"

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
	
	//if (output) std::cout<<CTMinvC<<std::endl<<"CTMinvCinvCTv"<<CTMinvCinvCTv<<std::endl<<"CTv "<<CTv<<std::endl;

	Eigen::Vector3d CCTMinvCinvCTv[3];
	CCTMinvCinvCTv[0] = CTMinvCinvCTv;
	CCTMinvCinvCTv[1] = -CTMinvCinvCTv;
	CCTMinvCinvCTv[2] = c.crossProductMatrix.transpose()*CTMinvCinvCTv;
	
	CCTMinvCinvCTv[0] *= femInvMass;
	CCTMinvCinvCTv[1] *= rbInvMass;
	CCTMinvCinvCTv[2] = Iinv * CCTMinvCinvCTv[2];
	

	//in[femIndex    ] = (femInvMass*in[femIndex   ] - CCTMinvCinvCTv[0][0]);
	//in[femIndex  +1] = (femInvMass*in[femIndex +1] - CCTMinvCinvCTv[0][1]);
	//in[femIndex  +2] = (femInvMass*in[femIndex +2] - CCTMinvCinvCTv[0][2]);
	//in[rbIndex    ] = (rbInvMass*in[rbIndex   ] - CCTMinvCinvCTv[1][0]);
	//in[rbIndex  +1] = (rbInvMass*in[rbIndex +1] - CCTMinvCinvCTv[1][1]);
	//in[rbIndex  +2] = (rbInvMass*in[rbIndex +2] - CCTMinvCinvCTv[1][2]);
	in[femIndex    ] -= CCTMinvCinvCTv[0][0];
	in[femIndex  +1] -= CCTMinvCinvCTv[0][1];
	in[femIndex  +2] -= CCTMinvCinvCTv[0][2];
	in[rbIndex    ] -= CCTMinvCinvCTv[1][0];
	in[rbIndex  +1] -= CCTMinvCinvCTv[1][1];
	in[rbIndex  +2] -= CCTMinvCinvCTv[1][2];
	//in[rbIndex  +3] -= CCTMinvCinvCTv[2][0];
	//in[rbIndex  +4] -= CCTMinvCinvCTv[2][1];
	//in[rbIndex  +5] -= CCTMinvCinvCTv[2][2];
	if (output) {
	//std::cout<<in[femIndex]<<" "<<CCTMinvCinvCTv[0][0]<<std::endl;
		//std::cout<<in[femIndex + 1]<<" "<<CCTMinvCinvCTv[0][1]<<std::endl;
	//std::cout<<in[femIndex + 2]<<" "<<CCTMinvCinvCTv[0][2]<<std::endl;
	//std::cout<<in[rbIndex + 0]<<" "<<CCTMinvCinvCTv[1][0]<<std::endl;
	//std::cout<<in[rbIndex + 1]<<" "<<CCTMinvCinvCTv[1][1]<<std::endl;
	//std::cout<<in[rbIndex + 2]<<" "<<CCTMinvCinvCTv[1][2]<<std::endl;
	//std::cout<<in[rbIndex + 3]<<" "<<CCTMinvCinvCTv[2][0]<<std::endl;
	//std::cout<<in[rbIndex + 4]<<" "<<CCTMinvCinvCTv[2][1]<<std::endl;
	//std::cout<<in[rbIndex + 5]<<" "<<CCTMinvCinvCTv[2][2]<<std::endl;
		std::cout<<in[rbIndex+0]-in[femIndex+0]<<" "<<in[rbIndex+1]-in[femIndex+1]<<" "<<in[rbIndex+2]-in[femIndex+2]<<std::endl;
		std::cout<<in[rbIndex+0]<<" "<<in[rbIndex+1]<<" "<<in[rbIndex+2]<<std::endl;
	std::cout<<"------------------"<<std::endl;
	}}

void World::project(double *in) {
	//return;
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
		//w->applyRegularizer(x.vals+offset+w->nrbdof, y.vals+offset+w->nrbdof);
		w->mulJV(x.vals, y.vals);
		w->mulJTLambda(x.vals, y.vals);
	}
private:
	World *w;
};

#if 1
void World::solve() {

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
			rhs[femIndex  ] -= f[0];
			rhs[femIndex+1] -= f[1];
			rhs[femIndex+2] -= f[2];
			rhs[rbIndex  ] += f[0];
			rhs[rbIndex+1] += f[1];
			rhs[rbIndex+2] += f[2];
			rhs[rbIndex+3] += t[0];
			rhs[rbIndex+4] += t[1];
			rhs[rbIndex+5] += t[2];
		}
		rbIndex+=6;
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

	//SimpleVector foo(n);
	//double *b = new double[n];
	//for (unsigned int i=0; i<n; i++)
	//b[i]=0.0;

	//op.Apply(x,foo);
	//mulJV(x.vals, b);	
	//for (unsigned int i=0; i<n; i++)
	//std::cout<<b[i]<<std::endl;
}

#else
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

	project(r);
	project(d);
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
			project(q);
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


void World::computeCrossProductMatrices(){

  for(auto& rb : rigidBodies){
    for(auto& c : rb.constraints){
      c.crossProductMatrix = 
	crossProductMatrix(quatRotate(rb.bulletBody->getOrientation(),c.localPosition));
    }
  }
}

#if 1
void World::mulJV(double* in, double*out){
	size_t rbIndex = nfemdof;
	double v[3];
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
#else
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
#endif

#if 1
//in = &(lambda[0])
//out = &(velocitiesToSolveFor[0])
void World::mulJTLambda(double* in, double* out){
}
#else
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
#endif

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
			out[constraintIndex    ] = 0.0;//4*c.weight*in[constraintIndex    ];
			out[constraintIndex + 1] = 0.0;//4*c.weight*in[constraintIndex + 1];
			out[constraintIndex + 2] = 0.0;//4*c.weight*in[constraintIndex + 2];
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
			out[constraintIndex    ] = 0.0;//in[constraintIndex    ] / (4*c.weight);
			out[constraintIndex + 1] = 0.0;//in[constraintIndex + 1] / (4*c.weight);
			out[constraintIndex + 2] = 0.0;//in[constraintIndex + 2] / (4*c.weight);
			constraintIndex += 3;
		}
	}  

  //auto recip = 1.0/regularizerAlpha;
  //for(auto i : range(3*totalNumConstraints)){
	//out[i] = recip*in[i];
  //}
}
