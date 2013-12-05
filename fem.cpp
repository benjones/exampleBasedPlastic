#include "fem.H"
#include <iostream>
#include <fstream>
#include <slIO.H>
#include <slUtil.H>
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#elif linux
extern "C" {
#include <cblas.h>
}
#endif
#include <float.h>
#include <sys/time.h>

using std::cout;
using std::endl;

//#ifdef __clang__
using std::unordered_map;
//#else
//using std::tr1::unordered_map
//#endif

//#define TIMING
//#define OUTPUT

bool extractSurface (const Tet* tets, unsigned int ntets, SlInt3 * &tris, unsigned int &ntris);

void FemObject::load(const char *fname) {
	std::ifstream in(fname, std::ios::in);
	nv = 0;
	ntets = 0;

	char ch;
	while (in>>ch) {
		if (ch == 'v') {
			nv++;
			double junk;
			in>>junk>>junk>>junk;
		} else if (ch == 't') {
			ntets++;
			int ijunk;
			double djunk;
			in>>ijunk>>ijunk>>ijunk>>ijunk>>djunk>>djunk>>djunk>>djunk;
		}
	}

	in.close();

	pos = new SlVector3[nv];
	vel = new SlVector3[nv];
	npos = new SlVector3[nv];
	nvel = new SlVector3[nv];
	bpos = new SlVector3[nv];
	bvel = new SlVector3[nv];
	flags = new unsigned char[nv];
	frc = new SlVector3[nv];
	mass = new double[nv];
	massInv = new double[nv];

	tets = new Tet[ntets];
	tetMass = new double[ntets];
	beta = new SlMatrix3x3[ntets];
	tetNormals = new TetNormals[ntets];
	k = new TetMatrix[ntets];
	Qx = new SlMatrix3x3[ntets];
	Qv = new SlMatrix3x3[ntets];

  density = new double[ntets];
  lambda = new double[ntets];
	mu = new double[ntets];
	scale = new double[ntets];

	frcOffsets = new FrcOffsets[ntets];

	bbMin = SlVector3(DBL_MAX);
	bbMax = SlVector3(-DBL_MAX);
	std::cout<<"number of vertices: "<<nv<<"; number of tets: "<<ntets<<std::endl;

	in.open(fname, std::ios::in);
	for (unsigned int i=0; i<nv; i++) {
		eatChar('v', in);
		in>>pos[i][0]>>pos[i][1]>>pos[i][2];
		bbMin.minSet(pos[i]);
		bbMax.maxSet(pos[i]);
	}

	std::cout<<"Bounding box: "<<bbMin<<" x "<<bbMax<<std::endl;

	for (unsigned int i=0; i<ntets; i++) {
		eatChar('t', in);
		in>>tets[i](0)>>tets[i](1)>>tets[i](2)>>tets[i](3)>>density[i]>>lambda[i]>>mu[i]>>scale[i];
	}
	in.close();

	for (unsigned int i=0; i<nv; i++) {
		vel[i] = 0.0;
		mass[i] = 0.0;
		flags[i] = 0;
	}

	SlVector3 normals[4];

	for (unsigned int t=0; t<ntets; t++) {
		double lambda = this->lambda[t];
		double mu = this->mu[t];
		//double scale = this->scale[t];
		double density = this->density[t];
		SlVector3 &x0  = pos[tets[t](0)];
		SlVector3 &x1  = pos[tets[t](1)];
		SlVector3 &x2  = pos[tets[t](2)];
		SlVector3 &x3  = pos[tets[t](3)];

		SlMatrix3x3 X(x1-x0, x2-x0, x3-x0, SlMatrix3x3::col);
		double det = determinant(X);
		if (det < 0) cout<<"inverted element"<<endl;
		SlMatrix3x3 &b = beta[t];
		tetMass[t] = density * det / 6.0;

		b = inverse(X);
		SlMatrix3x3 bt(b);
		bt.inplaceTranspose();
		normals[0] = cross(x2-x1, x3-x1);
		normals[1] = cross(x3-x0, x2-x0); 
		normals[2] = cross(x1-x0, x3-x0); 
		normals[3] = cross(x2-x0, x1-x0); 
		tetNormals[t][0] = normals[1] * 0.5;
		tetNormals[t][1] = normals[2] * 0.5;
		tetNormals[t][2] = normals[3] * 0.5;

		for (unsigned int i=0; i<4; i++) {
			SlVector3 &ni = normals[i];
			for (unsigned int j=0; j<=i; j++) {
				SlVector3 &nj = normals[j];
				SlMatrix3x3 nnt(ni[0]*nj[0], ni[0]*nj[1], ni[0]*nj[2], 
												ni[1]*nj[0], ni[1]*nj[1], ni[1]*nj[2], 
												ni[2]*nj[0], ni[2]*nj[1], ni[2]*nj[2]);
				SlMatrix3x3 &m = k[t](i,j);
				m = lambda*nnt + mu*transpose(nnt);
				double diag = mu*dot(ni,nj);
				m(0,0) += diag;
				m(1,1) += diag;
				m(2,2) += diag;
				m /= 6.0*det;
			}
		}
		frcOffsets[t][0] = -(k[t](0,0)*x0 + transmult(k[t](0,1), x1) + transmult(k[t](0,2), x2) + transmult(k[t](0,3), x3));
		frcOffsets[t][1] = -(k[t](1,0)*x0 + k[t](1,1)*x1 + transmult(k[t](1,2), x2) + transmult(k[t](1,3), x3));
		frcOffsets[t][2] = -(k[t](2,0)*x0 + k[t](2,1)*x1 + k[t](2,2)*x2 + transmult(k[t](2,3), x3));
		frcOffsets[t][3] = -(k[t](3,0)*x0 + k[t](3,1)*x1 + k[t](3,2)*x2 + k[t](3,3)*x3);

		for (unsigned int i=0; i<4; i++) {
			mass[tets[t](i)] += 0.25*tetMass[t];
		}
	}

	for (unsigned int i=0; i<nv; i++) massInv[i] = 1.0/mass[i];

	extractSurface(tets, ntets, tris, ntris);
	dumpObj("extract.obj");

	gm = new GlobalMatrix(nv, ntets, tets);

	solver_gripped = new bool[nv];
	solver_p = new SlVector3[nv];
	solver_r = new SlVector3[nv];
	solver_z = new SlVector3[nv];
	solver_s = new SlVector3[nv];
	solver_precon = new SlVector3[nv];
}

FemObject::FemObject() {
	friction = 0.5;
	active = true;
}

FemObject::~FemObject() {
	delete [] pos;
	delete [] vel;
	delete [] npos;
	delete [] nvel;
	delete [] bpos;
	delete [] bvel;
	delete [] frc;
	delete [] flags;
	delete [] mass;
	delete [] massInv;

	delete [] tets;
	delete [] tetMass;
	delete [] beta;
	delete [] tetNormals;
	delete [] k;
	delete [] lambda;
	delete [] mu;
	delete [] scale;

	delete [] solver_gripped;
	delete [] solver_p;
	delete [] solver_r;
	delete [] solver_z;
	delete [] solver_s;
	delete [] solver_precon;
}

inline double sqr(double x) {return x*x;};

inline SlMatrix3x3 localToGlobal(SlMatrix3x3 global, const SlMatrix3x3 &Qx, const SlMatrix3x3 &Qv, const SlMatrix3x3 &QxT, const SlMatrix3x3 &QvT, double dt, double scale){
	SlMatrix3x3 a(global);
	a.inplaceMultPost(QxT);
	a.inplaceMultPre(Qx);
	global.inplaceMultPost(QvT);
	global.inplaceMultPre(Qv);
	return (sqr(dt)*a)+((dt*scale)*global);
}

inline void addMass(SlMatrix3x3 &mat, const double &mass) {
	mat(0,0) += mass;  mat(1,1) += mass;  mat(2,2) += mass;
}

void FemObject::computeGlobalMatrix(double dt) {
	SlMatrix3x3 DQT;
	for (unsigned int i=0; i<nv; i++) {
		addMass(gm->find(i,i), mass[i]);
	}

	for(uint t=0; t < ntets; t++) {
		const Tet& tet = tets[t];
		SlMatrix3x3 QxT(Qx[t]);
		SlMatrix3x3 QvT(Qv[t]);
		QxT.inplaceTranspose();
		QvT.inplaceTranspose();

		double scale = this->scale[t];

		for(uint i=0; i < 4; ++i) {
			for(uint j=0; j <= i; ++j) {
				DQT = localToGlobal(k[t](i,j),Qx[t], Qv[t], QxT, QvT, dt, scale);
				if( tet(j) <= tet(i)) { 
					gm->find(tet(i), tet(j)) += DQT;
				} else {
					DQT.inplaceTranspose();
					gm->find(tet(j), tet(i)) += DQT;
				}
			}
		}
	}
}

void applyPrecon(SlVector3 *precon, SlVector3 *x, SlVector3 *y, unsigned int n) {
	SlVector3 *a=precon;
	SlVector3 *b=x;
	SlVector3 *c=y;
	for (unsigned int i=0; i<n; i++, a++, b++, c++) {
		(*c)[0] = (*a)[0]*(*b)[0];
		(*c)[1] = (*a)[1]*(*b)[1];
		(*c)[2] = (*a)[2]*(*b)[2];
	}
};

void FemObject::solveForVelocities(double dt) {
	double sigma_new, sigma, alpha, beta;
	//static const double tol=1e-12;
	double tol=1e-10;
	unsigned int iter=0, max_iter = 10000;
	unsigned int n = 3*(nv);

	SlVector3 *p = solver_p;
	SlVector3 *r = solver_r;
	SlVector3 *z = solver_z;
	SlVector3 *s = solver_s;
	SlVector3 *precon = solver_precon;

	gm->setPrecon(precon);

	SlVector3 *vptr = p;
	for (unsigned int i=0; i<nv; i++, vptr++) {
		(*vptr) = 0.0;
	}
	applygrips(p, dt);

	gm->mvmul(p,z);
	vptr = r;
	SlVector3 *vptr1 = z;
	SlVector3 *vptr2 = frc;
	for (unsigned int i=0; i<nv; i++, vptr++, vptr2++, vptr1++) {
		(*vptr) = (*vptr2) - (*vptr1);
	}

	filtergrips(r);

	applyPrecon(precon,r,s,nv);
	filtergrips(s);

	sigma = cblas_ddot(n, (double*)s, 1, (double*)r, 1);
	tol *= sigma;
	
	if (cblas_dnrm2(n, (double*)r, 1) >= tol) {
		while (iter++ < max_iter) {
			gm->mvmul(s,z);
			filtergrips(z);
			alpha = sigma/cblas_ddot(n, (double*)z, 1, (double*)s, 1);
			cblas_daxpy(n, alpha, (double*)s, 1, (double*)p, 1);
			cblas_daxpy(n, -alpha, (double*)z, 1, (double*)r, 1);
			
			if (cblas_dnrm2(n, (double*)r, 1) < tol) break;
			
			applyPrecon(precon,r,z,nv);
			sigma_new = cblas_ddot(n, (double*)z, 1, (double*)r, 1);
			beta = sigma_new / sigma;
			
			cblas_daxpy(n, beta, (double*)s, 1, (double*)z, 1);
			cblas_dcopy(n, (double*)z, 1, (double*)s, 1);
			filtergrips(s);

			sigma = sigma_new;
		}
	}

	vptr = p;
	vptr2 = vel;
	for (unsigned int i=0; i<nv; i++, vptr++, vptr2++) (*vptr2) = (*vptr);
}

void FemObject::updatePositions (double dt) {
	SlVector3 *vptr = pos;
	SlVector3 *vptr2 = vel;
	for (unsigned int i=0; i<nv; i++, vptr++, vptr2++) (*vptr) += dt*(*vptr2);
}

void FemObject::applygrips(SlVector3 *v, double dt) {
	bool *gptr = solver_gripped;
	for(unsigned int i = 0; i < nv; i++, gptr++){ *gptr = false; }

	std::vector<FemGrip*>::iterator gi;
	for(gi = grips.begin(); gi != grips.end(); gi++){
		(*gi)->flagGrips(solver_gripped);
	}

	if (v != NULL) {
		for(gi = grips.begin(); gi != grips.end(); gi++){
			(*gi)->applyGrips(v, dt);
		}	
	}
	
	gptr = solver_gripped;
	unsigned char *fptr = flags;
	for(unsigned int i = 0; i < nv; i++, gptr++, fptr++){
		*fptr = *fptr | (*gptr ? 4 : 0); // 4 = constrained
	}
};

void FemObject::filtergrips(SlVector3 *v){
	SlVector3 *vptr = v;
	const bool *gptr = solver_gripped;
	for(unsigned int i = 0; i < nv; i++, vptr++, gptr++){
		*vptr = ( (*gptr) ? 0.0 : *vptr );
	}
}

void FemObject::setForces(double dt, const SlVector3 &gravity) {
	SlVector3 *vptr = frc;
	double *dptr = mass;
	for (unsigned int i=0; i<nv; i++, vptr++, dptr++) {
		(*vptr) = (*dptr)*(dt*gravity + vel[i]);
	}
}

void FemObject::computeForces(double dt) {
	SlMatrix3x3 Fx, Fv, QTx, QTv, S;
	SlVector3 xtilde[4], fxtilde[4];
	for(uint t=0; t < ntets; ++t) {
		SlMatrix3x3 &QX = Qx[t];
		SlMatrix3x3 &QV = Qv[t];

		const Tet& tet = tets[t];
		SlVector3 & x0 = pos[tet(0)];
		SlVector3 & x1 = pos[tet(1)];
		SlVector3 & x2 = pos[tet(2)];
		SlVector3 & x3 = pos[tet(3)];
		SlMatrix3x3 X(x1 - x0, x2 - x0, x3 - x0, SlMatrix3x3::col);
		Fx = X * beta[t];
		
		SlVector3 & v0 = vel[tet(0)];
		SlVector3 & v1 = vel[tet(1)];
		SlVector3 & v2 = vel[tet(2)];
		SlVector3 & v3 = vel[tet(3)];
		SlMatrix3x3 XDot(v1 - v0, v2 - v0, v3 - v0, SlMatrix3x3::col);
		Fv = XDot * beta[t];
	
		polar(Fx, QX, S);
		polar(Fv, QV, S);

		QTx = QX;
		QTx.inplaceTranspose();
		QTv = QV;
		QTv.inplaceTranspose();

		xtilde[0] = QTx * x0;
		xtilde[1] = QTx * x1;
		xtilde[2] = QTx * x2;
		xtilde[3] = QTx * x3;

		k[t].mvmul(xtilde, fxtilde);

		frc[tet(0)] -= QX * (dt*(fxtilde[0] + frcOffsets[t][0]));
		frc[tet(1)] -= QX * (dt*(fxtilde[1] + frcOffsets[t][1]));
		frc[tet(2)] -= QX * (dt*(fxtilde[2] + frcOffsets[t][2]));
		frc[tet(3)] -= QX * (dt*(fxtilde[3] + frcOffsets[t][3]));
	}
}

FemSimulator::FemSimulator() {
	nobjects = 0;
	objects = NULL;
	computeForcesTime = 0.0;
	computeMatrixTime = 0.0;
	solveTime = 0.0;
	collisionTime = 0.0;
	nCollisionSurfaceVertices = 0;
}

FemSimulator::~FemSimulator() {
	if (objects) delete [] objects;
}

void FemSimulator::updatePositions(double dt) {
	for (unsigned int i=0; i<nobjects; i++) {
		objects[i].updatePositions(dt);
	}		
}

void FemSimulator::computeVelocities(double dt, const SlVector3 &gravity) {
	timeval startTime, endTime;

	for (unsigned int i=0; i<nobjects; i++) {
		if (objects[i].active) objects[i].setForces(dt, gravity);
		else objects[i].setForces(dt, 0.0);
	}
	for (unsigned int i=0; i<nobjects; i++) {
		gettimeofday(&startTime, NULL);
		objects[i].computeForces(dt);
		gettimeofday(&endTime, NULL);
		computeForcesTime += (endTime.tv_sec - startTime.tv_sec) +
			(endTime.tv_usec-startTime.tv_usec)*1.0e-6;
		
		objects[i].clearGlobalMatrix();
			
		gettimeofday(&startTime, NULL);
		objects[i].computeGlobalMatrix(dt);
		gettimeofday(&endTime, NULL);
		computeMatrixTime += (endTime.tv_sec - startTime.tv_sec) +
			(endTime.tv_usec-startTime.tv_usec)*1.0e-6;
		
		gettimeofday(&startTime, NULL);
		objects[i].solveForVelocities(dt);
		gettimeofday(&endTime, NULL);
		solveTime += (endTime.tv_sec - startTime.tv_sec) +
			(endTime.tv_usec-startTime.tv_usec)*1.0e-6;

	}
}

void FemObject::dumpObj(const char *fname) {
	std::ofstream out(fname, std::ios::out);

	for (unsigned int i=0; i<nv; i++) {
		out<<"v "<<pos[i][0]<<" "<<pos[i][1]<<" "<<pos[i][2]<<endl;
	}

	for (unsigned int i=0; i<ntris; i++) {
		out<<"f "<<tris[i][0]+1<<" "<<tris[i][1]+1<<" "<<tris[i][2]+1<<" "<<endl;
	}
}

void FemObject::dumpMesh(const char *fname) {
	std::ofstream out(fname, std::ios::out);

	for (unsigned int i=0; i<nv; i++) {
		out<<"v "<<pos[i][0]<<" "<<pos[i][1]<<" "<<pos[i][2]<<endl;
	}

	for (unsigned int i=0; i<ntets; i++) {
		out<<"t "<<tets[i](0)<<" "<<tets[i](1)<<" "<<tets[i](2)<<" "<<tets[i](3)<<endl;
	}
}

const unsigned char TetMatrix::index[4][4] = {{0,1,3,6},{1,2,4,7},{3,4,5,8},{6,7,8,9}};

struct eqSlInt3 {
  bool operator()(const SlInt3 &t1, const SlInt3 &t2) const {
    return ((t1[0] == t2[0] && t1[1] == t2[1] && t1[2] == t2[2]) ||
	    (t1[1] == t2[0] && t1[0] == t2[1] && t1[2] == t2[2]) ||
	    (t1[2] == t2[0] && t1[1] == t2[1] && t1[0] == t2[2]) ||
	    (t1[0] == t2[0] && t1[2] == t2[1] && t1[1] == t2[2]) ||
	    (t1[1] == t2[0] && t1[2] == t2[1] && t1[0] == t2[2]) ||
	    (t1[2] == t2[0] && t1[0] == t2[1] && t1[1] == t2[2]));
  }
};

struct hashSlInt3 {
  size_t operator()(const SlInt3& t) const {
    return (t[0] ^ t[1] ^ t[2]);
  }
};

typedef unordered_map<SlInt3, SlInt2, hashSlInt3, eqSlInt3> TriToTetMap;

inline void updateTriToTetMap(TriToTetMap &triToTetMap, int n, int x, int y, int z) {
  SlInt3 tri(x, y, z);
  if (triToTetMap.count(tri) == 0) {
    SlInt2 e(n, -1);
    triToTetMap.insert(TriToTetMap::value_type(tri,e));
  } else {
    triToTetMap.find(tri)->second[1] = n;
  }
}


bool extractSurface (const Tet* tets, unsigned int ntets, 
										 SlInt3 * &tris, unsigned int &ntris) {
  TriToTetMap triToTetMap;
  unsigned int i;
	std::vector<SlInt3> vtris;
	const Tet *t = tets;

  for (i=0; i<ntets; t++, i++) {
    updateTriToTetMap(triToTetMap, i, (*t)(0), (*t)(2), (*t)(1));
    updateTriToTetMap(triToTetMap, i, (*t)(3), (*t)(0), (*t)(1));
    updateTriToTetMap(triToTetMap, i, (*t)(3), (*t)(2), (*t)(0));
    updateTriToTetMap(triToTetMap, i, (*t)(2), (*t)(3), (*t)(1));
  }

  for (TriToTetMap::const_iterator m=triToTetMap.begin(); m!=triToTetMap.end(); m++) {
    if (m->second[1] == -1) {
      vtris.push_back(m->first);
    }
  }
	tris = new SlInt3[vtris.size()];
	for (unsigned int i=0; i<vtris.size(); i++) {
		tris[i] = vtris[i];
	}
	ntris = vtris.size();

  return true;
}

void FemObject::computeBoundingBox(SlVector3 &bbMin, SlVector3 &bbMax) {
	bbMax.set(-DBL_MAX, -DBL_MAX, -DBL_MAX);
	bbMin.set(DBL_MAX, DBL_MAX, DBL_MAX);

	for (unsigned int i=0; i<nv; i++) {
		bbMax.maxSet(pos[i]);
		bbMin.minSet(pos[i]);
	}
}


#if 0
void FemObject::collide(double dt) {
	double bound = 5.0;
	for (unsigned int i=0; i<nv; i++) {
		double x = pos[i][2]+dt*vel[i][2];
		if (x < 0.0) {
			double impulse = -x / dt;
			double nm=0.0,m = sqrt(sqr(vel[i][0]) + sqr(vel[i][1]));
			if (m > 0) nm = std::max<double>(m - impulse*friction, 0.0) / m;
			vel[i].set(nm*vel[i][0], nm*vel[i][1], vel[i][2]+impulse);
		}
	}
}


// code to dump frame data
timeSinceLastFrame+=dt;
	if (timeSinceLastFrame > framerate) {
		timeSinceLastFrame = 0;
		char fname[80];
#ifndef SUPPRESS_OBJ_OUTPUT
		for (unsigned int i=0; i<nobjects; i++) {
			sprintf(fname, framestring, i, frame);
			cout<<"writing frame "<<fname<<endl;
			objects[i].dumpObj(fname);
		}
#endif

		frame++;
#ifdef TIMING
		if (frame > 0)
			{
				cout<<totalTime/frame<<" "<<computeForcesTime/frame<<" "<<" "<<computeMatrixTime/frame<<" "<<solveTime/frame<<" "<<collisionTime/frame<<endl;
			}
#endif
	}

// code for dealing with obstacles
	for (unsigned int i=0; i<obstacles.size(); i++) {
		obstacles[i]->step(dt);
	}

		if (ground) objects[i].collide(dt);
		for (unsigned int j=0; j<obstacles.size(); j++) {
			obstacles[j]->collide(objects[i], dt);
		}


	time += dt;
#endif


