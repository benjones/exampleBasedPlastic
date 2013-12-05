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

class Tet {
public: 
	unsigned int vertices[4];

	inline unsigned int &operator()(const unsigned int &i) { return vertices[i]; }
	inline unsigned int operator()(const unsigned int &i) const { return vertices[i]; }
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
	void load(const char *fname);
	void loadNodeEle(const char *fname);
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
	void applyPlasticity(int t, double dt, SlMatrix3x3 &S);

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
	double *yieldStress, *plasticModulus, *alpha, *flowrate;  
	double alphaCap;

	double friction;

	SlInt3 *tris;  // surface triangles
	unsigned int ntris;
	void dumpObj(const char* fname);
	void dumpMesh(const char* fname);

	bool active;

  SlVector3 bbMin, bbMax;

private:
//Allocated space for computation
	bool *solver_gripped;
	SlVector3 *solver_p;
	SlVector3 *solver_r;
	SlVector3 *solver_z;
	SlVector3 *solver_s;
	SlVector3 *solver_precon;

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
	SlInt2 *collisionSurfaceToVertices;
	unsigned int nCollisionSurfaceVertices;
};

inline void TetMatrix::mvmul (const SlVector3 *x, SlVector3 *y) {
	y[0] = data[0]*x[0] + transmult(data[1],x[1]) + transmult(data[3],x[2]) + transmult(data[6],x[3]);
	y[1] = data[1]*x[0] + data[2]*x[1] + transmult(data[4],x[2]) + transmult(data[7],x[3]);
	y[2] = data[3]*x[0] + data[4]*x[1] + data[5]*x[2] + transmult(data[8],x[3]);
	y[3] = data[6]*x[0] + data[7]*x[1] + data[8]*x[2] + data[9]*x[3];
}



#endif

