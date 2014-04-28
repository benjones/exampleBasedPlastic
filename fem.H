#ifndef FEM_H
#define FEM_H

//use the c++11 versions if available
//#ifdef __clang__
#include <unordered_map>
#include <unordered_set>
//#else
//#include <tr1/unordered_map>
//#include <tr1/unordered_set>
//#endif

//#define HASH_MAP std::tr1::unordered_map
//#define HASH_SET std::tr1::unordered_set

#include "grip.H"
#include "obstacle.H"
#include <slVector.H>
#include <slMatrix.H>
#include <slUtil.H>
#include <assert.h>

//#include <dynamicsurface.h> //ignore collisions for now

//bullet stuff
#include <BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h>
#include <BulletCollision/Gimpact/btGImpactShape.h>
//#include <BulletCollision/CollisionShapes/btConvexHullShape.h>
#include <BulletCollision/CollisionShapes/btSphereShape.h>
#include <BulletDynamics/Dynamics/btRigidBody.h>
#include <LinearMath/btDefaultMotionState.h>

#include "cppitertools/range.hpp"
#include "utils.h"
using iter::range;

class btDiscreteDynamicsWorld;
class MaterialProperties;

/*class World;

struct BulletProxyHolder{
  BulletProxyHolder(btBroadphaseInterface* _interface, World& _world, 
					size_t _femIndex, size_t _nodeIndex);
  ~BulletProxyHolder();
  

  World& world;
  size_t femIndex, nodeIndex;
  btBroadphaseInterface* interface;
  std::unique_ptr<btBroadphaseProxy> proxy;
};
*/

class Tet {
public: 
	unsigned int vertices[4];

	inline unsigned int &operator()(const unsigned int &i) { return vertices[i]; }
	inline unsigned int operator()(const unsigned int &i) const { return vertices[i]; }


  //also sets linear and angular momentum
  /*void setInertiaTerms(SlVector3* pos,
					   SlVector3* vel,
					   double density,
					   btConvexHullShape* shape,
					   btRigidBody* body);
					   
  */
};

#include "globalMatrix.H"

class TetMatrix {
public:
	inline SlMatrix3x3 &operator()(const unsigned int &i, const unsigned int &j) { return data[index[i][j]];};
	inline SlMatrix3x3 operator()(const unsigned int &i, const unsigned int &j) const { return data[index[i][j]];};
	
	inline void mvmul (const SlVector3 *x, SlVector3 *y);
private:
	static const unsigned char index[4][4];
	SlMatrix3x3 data[10];
};


class FrcOffsets {
public:
	inline SlVector3 &operator()(const unsigned int &i) {return offsets[i];};
	inline SlVector3 operator()(const unsigned int &i) const {return offsets[i];};
	inline SlVector3 &operator[](const unsigned int &i) {return offsets[i];};
	inline SlVector3 operator[](const unsigned int &i) const {return offsets[i];};
private:
	SlVector3 offsets[4]; 
};

class TetNormals {
public:
	inline SlVector3 &operator()(const unsigned int &i) {return normals[i];};
	inline SlVector3 operator()(const unsigned int &i) const {return normals[i];};
	inline SlVector3 &operator[](const unsigned int &i) {return normals[i];};
	inline SlVector3 operator[](const unsigned int &i) const {return normals[i];};
private:
	SlVector3 normals[4]; 
};

class FemObject {
public:
	FemObject();
  
  FemObject(const FemObject& other) = delete; //no copying
  FemObject& operator=(const FemObject& other) = delete;

  //move only
  FemObject(FemObject&& other) = default;
  FemObject& operator=(FemObject&& other) = default;
  


  void load(const char *fname, const MaterialProperties &mp, btDiscreteDynamicsWorld* world);
	void loadNodeEle(const char *fname, const MaterialProperties &mp);
	~FemObject();
	void clearGlobalMatrix() {gm->clear();};
	void computeGlobalMatrix(double dt);
	void solveForVelocities(double dt);
	void updatePositions(double dt);
	void exsolve(double dt); // explicit integration
	void applygrips(SlVector3 *dv, double dt); // enforce constraints
	void filtergrips(SlVector3 *v);  // filtering constraints during implicit solve
	void computeForces(double dt);
	void collide(double dt);  // collide with the ground
	void setForces(double dt, const SlVector3 &gravity); // initialize force accumulator
	void computeBoundingBox(SlVector3 &bbMin, SlVector3 &bbMax);
	void computeK(int t);
	void applyPlasticity(int t, double dt, SlMatrix3x3 &stress, SlMatrix3x3 &S);
	void computeSeparationTensor(int t, double dt, const SlMatrix3x3 &S, const SlVector3 f[4]);
	void fracture();

	std::vector<FemGrip*> grips;  
	unsigned int nv, av; // number of vertices, allocated space
	SlVector3 *pos;  // position
	SlVector3 *vel;  // velocity
	SlVector3 *npos; // new position
	SlVector3 *nvel; // new velocity
	SlVector3 *bpos; // backup position
	SlVector3 *bvel; // backup velocity
	double* mass;    // node mass
	double* massInv; // inverse mass
	SlVector3 *frc;  // force accumulator
	unsigned char *flags; // 4 = constrained;
	inline bool gripped(int i){return ((flags[i]) & 4);}
	inline void grip(int i) {flags[i] &= 4;}

	unsigned int ntets;  // # tets
	Tet* tets; // tets

	double* tetMass; // tet mass
	SlMatrix3x3* beta;	// tet beta
	TetNormals *tetNormals; // tet normals
	TetMatrix *k;  // this is the stiffness matrix for the tet
	SlMatrix3x3 *Qx, *Qv;  // rotations for position and velocity (co-rotated formulation)

	//FrcOffsets *frcOffsets; // force offsets
	GlobalMatrix *gm;  // global system matrix

	double *density, *lambda, *mu, *scale;  // material parameters, phi = scale*lambda, psi = scale*mu
	double *yieldStress, *plasticModulus, *alpha, *flowrate;  // plasticity parameters
	double *toughness, *vertexToughness; // fracture parameters
	double alphaCap;

	double friction;

	SlVector3 *unbalancedCompLoad, *unbalancedTensLoad; // per-node fracture variables
	SlMatrix3x3 *separationTensor; // per-node fracture variable
	int *vertexToTetMap;
	int *vertexToTetMapStart;
	int *vertexToTetMapEnd;

	SlInt3 *tris;  // surface triangles
	unsigned int ntris;
	void dumpObj(const char* fname);
	void dumpMesh(const char* fname);

	bool active;

  size_t firstNodeIndex; //index of the first node in the global matrix with all FEM objects

  SlVector3 bbMin, bbMax;


  //void setupBulletMesh();
  //  void setupRigidProperties();
  //void updateBulletShapes();

  //  void stitchTets(); //stitch tets together after doing bullet collision stuff


  void setupBulletParticles();
  void copyStateToBulletParticles();
  void copyVelocitiesFromBulletParticles();
  
  void dumpVel(){
	std::cout << "vel" << std::endl;
	for(auto i : range(nv)){
	  std::cout << vel[i] << std::endl;
	}
  }

  void dumpPos(){
	std::cout << "pos" << std::endl;
	for(auto i : range(nv)){
	  std::cout << pos[i] << std::endl;
	}
  }

  void dumpTrans(){
	std::cout << "btrans" << std::endl;
	for(auto& tb : particleRigidBodies){
	  std::cout << tb->getCenterOfMassTransform() << std::endl;
	}
  }
  void dumpBulletVel(){
	std::cout << "bVel" << std::endl;
	for(auto& tb : particleRigidBodies){
	  std::cout << tb->getLinearVelocity() << std::endl;
	}
  }

  void dump(){
	dumpVel();
	dumpPos();
	dumpTrans();
	dumpBulletVel();
  }

  /*std::unique_ptr<btTriangleIndexVertexArray> bulletVertexArray;
  std::unique_ptr<btGImpactMeshShape> bulletMeshShape;
  std::unique_ptr<btDefaultMotionState> bulletMotionState;
  */

  //  std::vector<std::unique_ptr<btConvexHullShape>> bulletTetShapes;
  std::vector<std::unique_ptr<btRigidBody>> bulletTetBodies;
  std::vector<std::unique_ptr<btDefaultMotionState>> bulletMotionStates;
  //std::vector<std::unique_ptr<BulletProxyHolder>> bulletCollisionProxies;

  //this will be a 0 radius, 0 margin circle, tentatively
  std::unique_ptr<btCollisionShape> particleCollisionShape;
  std::vector<std::unique_ptr<btRigidBody>> particleRigidBodies;
  std::vector<std::unique_ptr<btDefaultMotionState>> particleMotionStates;

  btDiscreteDynamicsWorld* bulletWorld;
  double totalMass;


private:
//Allocated space for computation
	bool *solver_gripped;
	SlVector3 *solver_p;
	SlVector3 *solver_r;
	SlVector3 *solver_z;
	SlVector3 *solver_s;
	SlVector3 *solver_precon;

  std::vector<SlVector3> stitchedPositions;
  

};

class FemSimulator {
public:
	FemSimulator();
	~FemSimulator();
  //void initCollisions();
  //void handleSelfCollisions(double dt); //punt on this
	void computeVelocities(double dt, const SlVector3 &gravity);
	void updatePositions(double dt);

	unsigned int nobjects;
	FemObject *objects;

	double computeForcesTime, computeMatrixTime, solveTime, collisionTime;

private:
  //DynamicSurface *collisionSurface;
	//SlInt2 *collisionSurfaceToVertices;
	//unsigned int nCollisionSurfaceVertices;
};

class MaterialProperties {
public:
	double density, lambda, mu, scale, yieldStress, plasticModulus, flowrate, toughness;
};

inline void TetMatrix::mvmul (const SlVector3 *x, SlVector3 *y) {
	y[0] = data[0]*x[0] + transmult(data[1],x[1]) + transmult(data[3],x[2]) + transmult(data[6],x[3]);
	y[1] = data[1]*x[0] + data[2]*x[1] + transmult(data[4],x[2]) + transmult(data[7],x[3]);
	y[2] = data[3]*x[0] + data[4]*x[1] + data[5]*x[2] + transmult(data[8],x[3]);
	y[3] = data[6]*x[0] + data[7]*x[1] + data[8]*x[2] + data[9]*x[3];
}



#endif

