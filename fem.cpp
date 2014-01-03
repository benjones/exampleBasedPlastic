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

#include "cppitertools/range.hpp"
using iter::range;

#include <BulletDynamics/Dynamics/btDiscreteDynamicsWorld.h>

static const double FEPS = 1e-6;

bool extractSurface (const Tet* tets, unsigned int ntets, SlInt3 * &tris, unsigned int &ntris);

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
				DQT = localToGlobal(k[t](i,j), Qx[t], Qv[t], QxT, QvT, dt, scale);
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

void FemObject::computeK(const int t) { 
	for (unsigned int i=0; i<4; i++) {
		SlVector3 &ni = tetNormals[t][i];
		for (unsigned int j=0; j<=i; j++) {
			SlVector3 &nj = tetNormals[t][j];
			SlMatrix3x3 nnt(ni[0]*nj[0], ni[0]*nj[1], ni[0]*nj[2], 
											ni[1]*nj[0], ni[1]*nj[1], ni[1]*nj[2], 
											ni[2]*nj[0], ni[2]*nj[1], ni[2]*nj[2]);
			SlMatrix3x3 &m = k[t](i,j);
			m = lambda[t]*nnt + mu[t]*transpose(nnt);
			double diag = mu[t]*dot(ni,nj);
			m(0,0) += diag;
			m(1,1) += diag;
			m(2,2) += diag;
			m /= 36.0*tetMass[t]/density[t];
		}
	}
}

bool decomposeF(const SlMatrix3x3 &F, SlVector3 &Fhat, 
								SlMatrix3x3 &U, SlMatrix3x3 &V) {
  SlSVDecomp(F, U , Fhat,V);

  if (determinant(V) < 0) {
    V(0,0) = -V(0,0);
    V(1,0) = -V(1,0);
    V(2,0) = -V(2,0);
  }

  unsigned int decider = 0;
  if (Fhat[0] < FEPS) decider |= 1;
  if (Fhat[1] < FEPS) decider |= 2;
  if (Fhat[2] < FEPS) decider |= 4;
  SlVector3 U0, U1, U2;

  switch(decider) {
  case 0:
    U = F*V;
    U(0,0) /= Fhat[0];    U(0,1) /= Fhat[1];    U(0,2) /= Fhat[2];
    U(1,0) /= Fhat[0];    U(1,1) /= Fhat[1];    U(1,2) /= Fhat[2];
    U(2,0) /= Fhat[0];    U(2,1) /= Fhat[1];    U(2,2) /= Fhat[2];
    break;
  case 1:
    U1[0] = (F(0,0)*V(0,1) + F(0,1)*V(1,1) + F(0,2)*V(2,1))/Fhat[1]; 
    U1[1] = (F(1,0)*V(0,1) + F(1,1)*V(1,1) + F(1,2)*V(2,1))/Fhat[1]; 
    U1[2] = (F(2,0)*V(0,1) + F(2,1)*V(1,1) + F(2,2)*V(2,1))/Fhat[1]; 

    U2[0] = (F(0,0)*V(0,2) + F(0,1)*V(1,2) + F(0,2)*V(2,2))/Fhat[2]; 
    U2[1] = (F(1,0)*V(0,2) + F(1,1)*V(1,2) + F(1,2)*V(2,2))/Fhat[2]; 
    U2[2] = (F(2,0)*V(0,2) + F(2,1)*V(1,2) + F(2,2)*V(2,2))/Fhat[2]; 

    U0 = cross(U1, U2);
    U = SlMatrix3x3(U0,U1,U2,SlMatrix3x3::col);
    break;

  case 2:
    U0[0] = (F(0,0)*V(0,0) + F(0,1)*V(1,0) + F(0,2)*V(2,0))/Fhat[0]; 
    U0[1] = (F(1,0)*V(0,0) + F(1,1)*V(1,0) + F(1,2)*V(2,0))/Fhat[0]; 
    U0[2] = (F(2,0)*V(0,0) + F(2,1)*V(1,0) + F(2,2)*V(2,0))/Fhat[0]; 

    U2[0] = (F(0,0)*V(0,2) + F(0,1)*V(1,2) + F(0,2)*V(2,2))/Fhat[2]; 
    U2[1] = (F(1,0)*V(0,2) + F(1,1)*V(1,2) + F(1,2)*V(2,2))/Fhat[2]; 
    U2[2] = (F(2,0)*V(0,2) + F(2,1)*V(1,2) + F(2,2)*V(2,2))/Fhat[2]; 

    U1 = cross(U2, U0);
    U = SlMatrix3x3(U0,U1,U2,SlMatrix3x3::col);
    break;

  case 3:
    U2[0] = (F(0,0)*V(0,2) + F(0,1)*V(1,2) + F(0,2)*V(2,2))/Fhat[2]; 
    U2[1] = (F(1,0)*V(0,2) + F(1,1)*V(1,2) + F(1,2)*V(2,2))/Fhat[2]; 
    U2[2] = (F(2,0)*V(0,2) + F(2,1)*V(1,2) + F(2,2)*V(2,2))/Fhat[2]; 

    U0 = cross(U2, SlVector3(1,0,0));
    if (sqrMag(U0) < FEPS) U0 = cross(U2, SlVector3(0,1,0));
    normalize(U0);
    U1 = cross(U2, U0);
    U = SlMatrix3x3(U0,U1,U2,SlMatrix3x3::col);
    break;

  case 4:
    U0[0] = (F(0,0)*V(0,0) + F(0,1)*V(1,0) + F(0,2)*V(2,0))/Fhat[0]; 
    U0[1] = (F(1,0)*V(0,0) + F(1,1)*V(1,0) + F(1,2)*V(2,0))/Fhat[0]; 
    U0[2] = (F(2,0)*V(0,0) + F(2,1)*V(1,0) + F(2,2)*V(2,0))/Fhat[0]; 

    U1[0] = (F(0,0)*V(0,1) + F(0,1)*V(1,1) + F(0,2)*V(2,1))/Fhat[1]; 
    U1[1] = (F(1,0)*V(0,1) + F(1,1)*V(1,1) + F(1,2)*V(2,1))/Fhat[1]; 
    U1[2] = (F(2,0)*V(0,1) + F(2,1)*V(1,1) + F(2,2)*V(2,1))/Fhat[1]; 

    U2 = cross(U0, U1);
    U = SlMatrix3x3(U0,U1,U2,SlMatrix3x3::col);
    break;

  case 5:
    U1[0] = (F(0,0)*V(0,1) + F(0,1)*V(1,1) + F(0,2)*V(2,1))/Fhat[1]; 
    U1[1] = (F(1,0)*V(0,1) + F(1,1)*V(1,1) + F(1,2)*V(2,1))/Fhat[1]; 
    U1[2] = (F(2,0)*V(0,1) + F(2,1)*V(1,1) + F(2,2)*V(2,1))/Fhat[1]; 

    U2 = cross(U1, SlVector3(1,0,0));
    if (sqrMag(U2) < FEPS) U2 = cross(U1, SlVector3(0,1,0));
    normalize(U2);
    U0 = cross(U1, U2);
    U = SlMatrix3x3(U0,U1,U2,SlMatrix3x3::col);
    break;

  case 6:
    U0[0] = (F(0,0)*V(0,0) + F(0,1)*V(1,0) + F(0,2)*V(2,0))/Fhat[0]; 
    U0[1] = (F(1,0)*V(0,0) + F(1,1)*V(1,0) + F(1,2)*V(2,0))/Fhat[0]; 
    U0[2] = (F(2,0)*V(0,0) + F(2,1)*V(1,0) + F(2,2)*V(2,0))/Fhat[0]; 

    U1 = cross(U0, SlVector3(1,0,0));
    if (sqrMag(U1) < FEPS) U1 = cross(U0, SlVector3(0,1,0));
    normalize(U1);
    U2 = cross(U0, U1);
    U = SlMatrix3x3(U0,U1,U2,SlMatrix3x3::col);
    break;

  case 7:
    cout<<"case 7"<<endl;
    U = SlMatrix3x3(1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0);
    break;

  default:
    cout<<"oh no"<<endl;
    exit(-1);
  }
    

  if (determinant(F) < 0) {
    if (Fhat[0] < Fhat[1] && Fhat[0] < Fhat[2]) {
      Fhat[0] = -Fhat[0];
      U(0,0) = -U(0,0);
      U(1,0) = -U(1,0);
      U(2,0) = -U(2,0);
    } else if (Fhat[1] < Fhat[2]) {
      Fhat[1] = -Fhat[1];
      U(0,1) = -U(0,1);
      U(1,1) = -U(1,1);
      U(2,1) = -U(2,1);
    } else {
      Fhat[2] = -Fhat[2];
      U(0,2) = -U(0,2);
      U(1,2) = -U(1,2);
      U(2,2) = -U(2,2);
    }
  }
  return true;
}

// somewhere this should probably be S+I or something...
void FemObject::applyPlasticity(int t, double dt, SlMatrix3x3 &stress, SlMatrix3x3 &S) {
	SlVector3 Fhat;
	SlMatrix3x3 V;
	double m = norm(stress);
	double thresh = m - yieldStress[t] - plasticModulus[t]*std::min(alpha[t], alphaCap);
	//std::cout<<S;
	//std::cout<<m<<" "<<yieldStress[t]<<std::endl;
	if (thresh < 0.0) return;
	SlSymetricEigenDecomp(S, Fhat, V);

	if (Fhat[0] < FEPS || Fhat[1] < FEPS || Fhat[2] < FEPS) return;

	SlMatrix3x3 VT(V);
	VT.inplaceTranspose();
	SlVector3 Fstarhat = Fhat / cbrt(Fhat[0]*Fhat[1]*Fhat[2]);
	double gamma = std::min(std::min(flowrate[t] * dt, 1.0) * thresh / m, 1.0);
	SlVector3 Fphat(pow(Fstarhat[0], gamma), pow(Fstarhat[1], gamma), pow(Fstarhat[2], gamma));
	beta[t] = beta[t] * V * diagonal(SlVector3(1.0/Fphat[0], 1.0/Fphat[1], 1.0/Fphat[2])) * VT;
	alpha[t] += dt*norm(V*diagonal(Fphat)*VT);
	Fhat[0] /= Fphat[0];
	Fhat[1] /= Fphat[1];
	Fhat[2] /= Fphat[2];
	S = V*diagonal(Fhat)*VT;
	
	// recompute stress
	stress = S;
	stress(0,0) -= 1.0; stress(1,1) -= 1.0; stress(2,2) -= 1.0;
	double trace = lambda[t]*(stress(0,0)+stress(1,1)+stress(2,2));
	stress *= 2.0*mu[t];
	stress(0,0) += trace; stress(1,1) += trace; stress(2,2) += trace;

	// compute new tetnormals
	SlMatrix3x3 B (inverse(beta[t]));
	
	tetNormals[t][0] = 0.5*(cross(SlVector3(B(0,1)-B(0,0), B(1,1)-B(1,0), B(2,1)-B(2,0)),
																SlVector3(B(0,2)-B(0,0), B(1,2)-B(1,0), B(2,2)-B(2,0))));
	tetNormals[t][1] = 0.5*(cross(SlVector3(B(0,2),B(1,2),B(2,2)),
																SlVector3(B(0,1),B(1,1),B(2,1))));
	tetNormals[t][2] = 0.5*(cross(SlVector3(B(0,0),B(1,0),B(2,0)),
																SlVector3(B(0,2),B(1,2),B(2,2))));
	tetNormals[t][3] = 0.5*(cross(SlVector3(B(0,1),B(1,1),B(2,1)),
																SlVector3(B(0,0),B(1,0),B(2,0))));
	computeK(t);
}

void FemObject::computeForces(double dt) {
	SlMatrix3x3 Fx, Fv, QTx, QTv, S, U, V;
	SlVector3 Fhat;


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
	
		polar(Fv, QV, S);
		polar(Fx, QX, S);
		bool inverted = false;
		if (determinant(S) < FEPS) {
			decomposeF(Fx, Fhat, U, V);
			QX = U*transpose(V);
			S = V*diagonal(Fhat)*transpose(V); 
			inverted = true;
		}

		SlMatrix3x3 stress(S);
		stress(0,0) -= 1.0; stress(1,1) -= 1.0; stress(2,2) -= 1.0;
		double trace = lambda[t]*(stress(0,0)+stress(1,1)+stress(2,2));
		stress *= 2.0*mu[t];
		stress(0,0) += trace; stress(1,1) += trace; stress(2,2) += trace;

		if (!inverted) {
			applyPlasticity(t, dt, stress, S);
		}

		SlMatrix3x3 QS(QX * stress);
		SlVector3 f[4];

		f[0] = QS * (dt*(tetNormals[t][0]));
		f[1] = QS * (dt*(tetNormals[t][1]));
		f[2] = QS * (dt*(tetNormals[t][2]));
		f[3] = QS * (dt*(tetNormals[t][3]));
		frc[tet(0)] += f[0];
		frc[tet(1)] += f[1];
		frc[tet(2)] += f[2];
		frc[tet(3)] += f[3];
		computeSeparationTensor(t, dt, stress, f);
		//if (t == 0) std::cout<<frc[tet(0)]<<std::endl;
	}
}

inline void M(const SlVector3 &x, const double &val, SlMatrix3x3 &m) {
	SlVector3 y(val*x);
	m(0,0) = x[0]*y[0];
	m(0,1) = m(1,0) = x[0]*y[1];
	m(0,2) = m(2,0) = x[0]*y[2];
	m(1,1) = x[1]*y[1];
	m(1,2) = m(2,1) = x[2]*y[1];
	m(2,2) = x[2]*y[2];
}

void FemObject::computeSeparationTensor(int t, double dt, const SlMatrix3x3 &S, const SlVector3 f[4]) {
	static const double EPS = 1e-6;
	SlVector3 vals;
	SlMatrix3x3 vecs, compression(0.0), tension(0.0), m;
	SlVector3 tensionF, compressF;
	double scale;

	SlSymetricEigenDecomp(S, vals, vecs);

	for (unsigned int i=0; i<3; i++) {
		if (vals[i] > -EPS && vals[i] < EPS) continue;
		M(SlVector3(vecs(0,i), vecs(1,i), vecs(2,i)), vals[i], m);
		(vals[i] > 0) ? tension += m : compression += m;
	}

	SlMatrix3x3 QS(Qx[t] * tension);
		
	for (unsigned int i=0; i<4; i++) {
		tensionF = QS * tetNormals[t][i];
		compressF = (f[i]/dt) - tensionF;
		scale = mag(tensionF);
		if (scale > EPS) {
			scale = 1.0/scale;
			M(tensionF, scale, m);
			separationTensor[tets[t](i)] += m;
		}
		scale = mag(compressF);
		if (scale > EPS) {
			scale = 1.0/scale;
			M(compressF, scale, m);
			separationTensor[tets[t](i)] -= m;
		}
		unbalancedTensLoad[tets[t](i)] += tensionF;
		unbalancedCompLoad[tets[t](i)] += compressF;
	}
}

inline void swap(int *array, int i1, int i2) {
	int tmp = array[i1];
	array[i1] = array[i2];
	array[i2] = tmp;
}

void FemObject::fracture() {
	static const double EPS = 1e-6;
	if (av == nv) {
		std::cout<<"Warning: we've stopped fracture"<<std::endl;
		return;
	}

	SlVector3 vals;
	SlMatrix3x3 vecs;
	int maxi;
	double maxval;
	int oldnv = nv;
	SlMatrix3x3 m;
	double scale;

	for (unsigned int v=0; v<oldnv; v++) {
		if (gripped(v)) continue;
		scale = mag(unbalancedTensLoad[v]);
		if (scale > EPS) {
			M(unbalancedTensLoad[v], 1.0/scale, m);
			separationTensor[v] -= m;
		}
		scale = mag(unbalancedCompLoad[v]);
		if (scale > EPS) {
			M(unbalancedCompLoad[v], 1.0/scale, m);
			separationTensor[v] += m;
		}

		SlSymetricEigenDecomp(separationTensor[v], vals, vecs);
		maxval = vals[0];
		maxi = 0;
		if (vals[1] > maxval) {
			maxval = vals[1];
			maxi = 1;
		} 
		if (vals[2] > maxval) {
			maxval = vals[2];
			maxi = 2;
		}

		if (maxval > vertexToughness[v]) {
			SlVector3 n(vecs(0,maxi),vecs(1,maxi),vecs(2,maxi));
			int below = vertexToTetMapStart[v];
			int above = vertexToTetMapEnd[v]-1;
			int pointabove[4];
			double massa=0.0, massb=0.0;
			double toughnessa=0.0, toughnessb=0.0;
			int vertex=-1;
			while (above >= below) {
				int pointsaboveplane = 0;
				int t = vertexToTetMap[below];
				double offset = dot(n, pos[v]);
				for (unsigned int j=0; j<4; j++) {
					if (tets[t](j) == v) {
						pointabove[j] = 0;
						vertex = j;
					} else if (dot(n, pos[tets[t](j)]) > offset) {
						pointsaboveplane++;
						pointabove[j] = 1;
						massa += tetMass[t];
						toughnessa += tetMass[t]*toughness[t];
					} else {
						pointabove[j] = -1;
						massb += tetMass[t];
						toughnessb += tetMass[t]*toughness[t];
					}
				}

				if (pointsaboveplane == 3) {
					tets[t](vertex) = nv;
					swap(vertexToTetMap, below, above--);
				} else if (pointsaboveplane == 0) {
					below++;
				} else {
					// determine whether tet should be above or below.  geometric nightmare.  for now we vote.
					if (pointsaboveplane > 1) {
						tets[t](vertex) = nv;
						swap(vertexToTetMap, below, above--);
					} else {
						below++;
					}
				}
			}
			if (below == vertexToTetMapEnd[v] || below == vertexToTetMapStart[v]) {
				for (unsigned int i=vertexToTetMapStart[v]; i<vertexToTetMapEnd[v]; i++) {
					int t = vertexToTetMap[i];
					for (unsigned int j=0; j<4; j++) {
						if (tets[t](j) == nv) {
							tets[t](j) = v;
						}
					}
				}
				continue;
			}
			vertexToTetMapStart[nv] = below;
			vertexToTetMapEnd[nv] = vertexToTetMapEnd[v];
			vertexToughness[nv] = toughnessa / massa;
			mass[nv] = massa * 0.25;
			massInv[nv] = 1.0/mass[nv];
			vertexToTetMapEnd[v] = below;
			vertexToughness[v] = toughnessb / massb;
			mass[v] = massb * 0.25;
			massInv[v] = 1.0/mass[v];
			pos[nv] = pos[v];
			vel[nv] = vel[v];
			flags[nv] = flags[v];
			nv++;
			if (av == nv) {
				std::cout<<"Warning: we've stopped fracture"<<std::endl;
				return;
			}
		}
	}
	if (oldnv != nv) {
		std::cout<<"fracture! oldnv = "<<oldnv<<" nv = "<<nv<<std::endl;
		delete gm;
		gm = new GlobalMatrix(nv, ntets, tets);
		extractSurface(tets, ntets, tris, ntris);
	}
}

FemSimulator::FemSimulator() {
	nobjects = 0;
	objects = NULL;
	computeForcesTime = 0.0;
	computeMatrixTime = 0.0;
	solveTime = 0.0;
	collisionTime = 0.0;
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
		for (unsigned int v=0; v<objects[i].nv; v++) {
			objects[i].separationTensor[v] = 0.0;
			objects[i].unbalancedTensLoad[v] = 0.0;
			objects[i].unbalancedCompLoad[v] = 0.0;
		}
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
	if (tris) delete [] tris;
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


/*****************************************************************************/
/*                                                                           */
/*  inputtextline()   Read a nonempty line from a file.                      */
/*                                                                           */
/*  A line is considered "nonempty" if it contains something that looks like */
/*  a number.  Comments (prefaced by `#') are ignored.                       */
/*                                                                           */
/*****************************************************************************/


char *inputtextline(char *string,
                    FILE *infile,
                    char *infilename)
{
  char *result;

  /* Search for something that looks like a number. */
  do {
    result = fgets(string, 1024, infile);
    if (result == (char *) NULL) {
      printf("  Error:  Unexpected end of file in %s.\n", infilename);
      exit(1);
    }

    /* Skip anything that doesn't look like a number, a comment, */
    /*   or the end of a line.                                   */
    while ((*result != '\0') && (*result != '#') &&
           (*result != '.') && (*result != '+') && (*result != '-') &&
           ((*result < '0') || (*result > '9'))) {
      result++;
    }
  /* If it's a comment or end of line, read another line and try again. */
  } while ((*result == '#') || (*result == '\0'));

  return result;
}


/*****************************************************************************/
/*                                                                           */
/*  inputfindfield()   Find the next field of a string.                      */
/*                                                                           */
/*  Jumps past the current field by searching for whitespace, then jumps     */
/*  past the whitespace to find the next field.                              */
/*                                                                           */
/*****************************************************************************/


char *inputfindfield(char *string)
{
  char *result;

  result = string;
  /* Skip the current field.  Stop upon reaching whitespace. */
  while ((*result != '\0') && (*result != '#') &&
         (*result != ' ') && (*result != '\t')) {
    result++;
  }

  /* Now skip the whitespace and anything else that doesn't look like a */
  /*   number, a comment, or the end of a line.                         */
  while ((*result != '\0') && (*result != '#') &&
         (*result != '.') && (*result != '+') && (*result != '-') &&
         ((*result < '0') || (*result > '9'))) {
    result++;
  }

  /* Check for a comment (prefixed with `#'). */
  if (*result == '#') {
    *result = '\0';
  }

  return result;
}

void FemObject::loadNodeEle(const char *fname) {
  totalMass = 0;
  
  char inputline[1024];
  char *stringptr;
  FILE *infile = NULL;
  unsigned int dim = 0;
	unsigned int nnodeattributes = 0;
	unsigned int neleattributes = 0;
  int markflag = 0;
	char nodefilename[80];
	char elefilename[80];
	unsigned int firstindex = 0;

	sprintf(nodefilename, "%s.node", fname);
	sprintf(elefilename, "%s.ele", fname);

	/* Read the vertices from a .node file. */
	printf("Opening %s.\n", nodefilename);

	infile = fopen(nodefilename, "r");
	if (infile == (FILE *) NULL) {
		printf("  Error:  Cannot access file %s.\n", nodefilename);
		exit(1);
	}

	/* Read number of vertices, number of dimensions, number of vertex */
	/*   attributes, and number of boundary markers.                   */
	stringptr = inputtextline(inputline, infile, nodefilename);
	nv = (int) strtol(stringptr, &stringptr, 0);
	av = 2*nv;

	stringptr = inputfindfield(stringptr);
	if (*stringptr == '\0') {
		dim = 3;
	} else {
		dim = (int) strtol(stringptr, &stringptr, 0);
	}
	stringptr = inputfindfield(stringptr);
	if (*stringptr == '\0') {
		nnodeattributes = 0;
	} else {
		nnodeattributes = (unsigned int) strtol(stringptr, &stringptr, 0);
	}
	stringptr = inputfindfield(stringptr);
	if (*stringptr == '\0') {
		markflag = 0;
	} else {
		markflag = (int) strtol(stringptr, &stringptr, 0);
	}

  if (nv < 4) {
    printf("Error:  Input must have at least four input vertices.\n");
    exit(1);
  }
  if (dim != 3) {
    printf("Error:  FEM only works with three-dimensional meshes.\n");
    exit(1);
  }

	pos = new SlVector3[av];
	vel = new SlVector3[av];
	npos = new SlVector3[av];
	nvel = new SlVector3[av];
	bpos = new SlVector3[av];
	bvel = new SlVector3[av];
	flags = new unsigned char[av];
	frc = new SlVector3[av];
	mass = new double[av];
	massInv = new double[av];
	vertexToughness = new double[av];
	separationTensor = new SlMatrix3x3[av];
	unbalancedTensLoad = new SlVector3[av];
	unbalancedCompLoad = new SlVector3[av];

	SlVector3 bbMin(DBL_MAX), bbMax(-DBL_MAX);
	std::cout<<"number of vertices: "<<nv<<"; number of tets: "<<ntets<<std::endl;
  for (unsigned int i = 0; i < nv; i++) {
    stringptr = inputtextline(inputline, infile, nodefilename);
		if (i == 0) firstindex = strtol(stringptr, &stringptr, 0);
		
    /* Read the vertex coordinates and attributes. */
    for (unsigned int j = 0; j < 3 + nnodeattributes; j++) {
      stringptr = inputfindfield(stringptr);
      if (*stringptr == '\0') {
				if (j < 3) {
          printf("Error:  Vertex %lu has no %c coordinate.\n",
                 (unsigned long) (i),
                 (j == 0) ? 'x' : (j == 1) ? 'y' : 'z');
          exit(1);
        }
      } else {
        pos[i][j] = strtod(stringptr, &stringptr);
      }
    }
		bbMin.minSet(pos[i]);
		bbMax.maxSet(pos[i]);

    if (markflag) {
      /* Read a vertex marker. */
      stringptr = inputfindfield(stringptr);
    }
  }

	fclose(infile);

	std::cout<<"Bounding box: "<<bbMin<<" x "<<bbMax<<std::endl;

  /* Read the tetrahedra from an .ele file. */
	printf("Opening %s.\n", elefilename);

  infile = fopen(elefilename, "r");
  if (infile == (FILE *) NULL) {
    printf("  Error:  Cannot access file %s.\n", elefilename);
    exit(1);
  }

  /* Read number of tetrahedra, number of vertices per tetrahedron, and */
  /*   number of tetrahedron attributes from .ele file.                 */
  stringptr = inputtextline(inputline, infile, elefilename);
  ntets = (int) strtol(stringptr, &stringptr, 0);

  stringptr = inputfindfield(stringptr);
  if (*stringptr == '\0') {
  } else {
    int elemnodes = (int) strtol(stringptr, &stringptr, 0);
    if (elemnodes < 4) {
      printf("Error:  Tetrahedra in %s must have at least 4 vertices.\n",
             elefilename);
      exit(1);
    }
  }
  stringptr = inputfindfield(stringptr);
  if (*stringptr == '\0') {
    neleattributes = 0;
  } else {
    neleattributes = (unsigned int) strtol(stringptr, &stringptr, 0);
  }

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
	yieldStress = new double[ntets];
	plasticModulus = new double[ntets];
	alpha = new double[ntets];
	flowrate = new double[ntets];
	toughness = new double[ntets];

  for (unsigned int i = 0; i < ntets; i++) {
    /* Read the tetrahedron's four vertices. */
    stringptr = inputtextline(inputline, infile, elefilename);
    for (unsigned int j = 0; j < 4; j++) {
      stringptr = inputfindfield(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Tetrahedron %lu is missing vertex %d in %s.\n",
               (unsigned long) i, j + 1, elefilename);
        exit(1);
      } else {
        tets[i](j) = (int) strtol(stringptr, &stringptr, 0) - firstindex;
      }
    }
		if (neleattributes == 0) {
			density[i] = 1000;
			lambda[i] = 5e4;
			mu[i] = 1e5;
			scale[i] = 1;
			yieldStress[i] = DBL_MAX;//10.0;
			plasticModulus[i] = 0.0;
			flowrate[i] = 100;
			toughness[i] = DBL_MAX;//100.0;
		} else if (neleattributes == 4) {
			stringptr = inputfindfield(stringptr);
			density[i] = strtod(stringptr, &stringptr);
			stringptr = inputfindfield(stringptr);
			lambda[i] = strtod(stringptr, &stringptr);
			stringptr = inputfindfield(stringptr);
			mu[i] = strtod(stringptr, &stringptr);
			stringptr = inputfindfield(stringptr);
			scale[i] = strtod(stringptr, &stringptr);
			yieldStress[i] = DBL_MAX;
			plasticModulus[i] = 0.0;
			flowrate[i] = 0.0;
			toughness[i] = 100.0;
		} else if (neleattributes == 7) {
			stringptr = inputfindfield(stringptr);
			density[i] = strtod(stringptr, &stringptr);
			stringptr = inputfindfield(stringptr);
			lambda[i] = strtod(stringptr, &stringptr);
			stringptr = inputfindfield(stringptr);
			mu[i] = strtod(stringptr, &stringptr);
			stringptr = inputfindfield(stringptr);
			scale[i] = strtod(stringptr, &stringptr);
			stringptr = inputfindfield(stringptr);
			yieldStress[i] = strtod(stringptr, &stringptr);
			stringptr = inputfindfield(stringptr);
			plasticModulus[i] = strtod(stringptr, &stringptr);
			stringptr = inputfindfield(stringptr);
			flowrate[i] = strtod(stringptr, &stringptr);
		} else if (neleattributes == 8) {
			stringptr = inputfindfield(stringptr);
			density[i] = strtod(stringptr, &stringptr);
			stringptr = inputfindfield(stringptr);
			lambda[i] = strtod(stringptr, &stringptr);
			stringptr = inputfindfield(stringptr);
			mu[i] = strtod(stringptr, &stringptr);
			stringptr = inputfindfield(stringptr);
			scale[i] = strtod(stringptr, &stringptr);
			stringptr = inputfindfield(stringptr);
			yieldStress[i] = strtod(stringptr, &stringptr);
			stringptr = inputfindfield(stringptr);
			plasticModulus[i] = strtod(stringptr, &stringptr);
			stringptr = inputfindfield(stringptr);
			flowrate[i] = strtod(stringptr, &stringptr);
			stringptr = inputfindfield(stringptr);
			toughness[i] = strtod(stringptr, &stringptr);
		}
  }
  fclose(infile);

	for (unsigned int i=0; i<nv; i++) {
		vel[i] = 0.0;
		mass[i] = 0.0;
		flags[i] = 0;
	}

	for (unsigned int t=0; t<ntets; t++) {
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

		totalMass += tetMass[t];

		b = inverse(X);
		SlMatrix3x3 bt(b);
		bt.inplaceTranspose();
		tetNormals[t][0] = 0.5 * cross(x2-x1, x3-x1);;
		tetNormals[t][1] = 0.5 * cross(x3-x0, x2-x0); 
		tetNormals[t][2] = 0.5 * cross(x1-x0, x3-x0); 
		tetNormals[t][3] = 0.5 * cross(x2-x0, x1-x0); 

		for (unsigned int i=0; i<4; i++) {
			mass[tets[t](i)] += 0.25*tetMass[t];
			vertexToughness[tets[t](i)] += 0.25*tetMass[t]*toughness[t];
		}
		computeK(t);
		alpha[t] = 0.0;
	}

	for (unsigned int i=0; i<nv; i++) {
		massInv[i] = 1.0/mass[i];
		vertexToughness[i] *= massInv[i];
	}

	// build vertex to tet map
	vertexToTetMap = new int[4*ntets];
	vertexToTetMapStart = new int[2*nv];
	vertexToTetMapEnd = new int[2*nv];
	std::vector<int> *tmpVertexToTetMap = new std::vector<int>[nv];
	for (unsigned int i=0; i<ntets; i++) {
		tmpVertexToTetMap[tets[i](0)].push_back(i);
		tmpVertexToTetMap[tets[i](1)].push_back(i);
		tmpVertexToTetMap[tets[i](2)].push_back(i);
		tmpVertexToTetMap[tets[i](3)].push_back(i);
	}
	int index = 0;
	for (unsigned int i=0; i<nv; i++) {
		vertexToTetMapStart[i] = index;
		for (unsigned int j=0; j<tmpVertexToTetMap[i].size(); j++) {
			vertexToTetMap[index++] = tmpVertexToTetMap[i][j];
		}
		vertexToTetMapEnd[i] = index;
	}
	delete [] tmpVertexToTetMap;

	extractSurface(tets, ntets, tris, ntris);
	dumpObj("extract.obj");
	setupBulletMesh();

	gm = new GlobalMatrix(nv, ntets, tets);

	solver_gripped = new bool[av];
	solver_p = new SlVector3[av];
	solver_r = new SlVector3[av];
	solver_z = new SlVector3[av];
	solver_s = new SlVector3[av];
	solver_precon = new SlVector3[av];
}

void FemObject::load(const char *fname, btDiscreteDynamicsWorld* world) {
  bulletWorld = world;

	if (!strstr(fname, ".mesh")) return loadNodeEle(fname);

	totalMass = 0;
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
	yieldStress = new double[ntets];
	plasticModulus = new double[ntets];
	alpha = new double[ntets];
	flowrate = new double[ntets];

	//frcOffsets = new FrcOffsets[ntets];

	SlVector3 bbMin(DBL_MAX), bbMax(-DBL_MAX);
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

	for (unsigned int t=0; t<ntets; t++) {
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
		
		totalMass += tetMass[t];


		b = inverse(X);
		SlMatrix3x3 bt(b);
		bt.inplaceTranspose();

		for (unsigned int i=0; i<4; i++) {
			mass[tets[t](i)] += 0.25*tetMass[t];
		}
		computeK(t);
		alpha[t] = 0.0;
	}

	for (unsigned int i=0; i<nv; i++) massInv[i] = 1.0/mass[i];

	extractSurface(tets, ntets, tris, ntris);
	dumpObj("extract.obj");

	setupBulletMesh();


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
	tris = NULL;
}

FemObject::~FemObject() {
	delete gm;
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
	delete [] yieldStress;
	delete [] plasticModulus;
	delete [] alpha;
	delete [] flowrate;
	delete [] toughness;
	delete [] vertexToughness;
	delete [] separationTensor;
	delete [] unbalancedTensLoad;
	delete [] unbalancedCompLoad;
	delete [] vertexToTetMap;
	delete [] vertexToTetMapStart;
	delete [] vertexToTetMapEnd;

	delete [] solver_gripped;
	delete [] solver_p;
	delete [] solver_r;
	delete [] solver_z;
	delete [] solver_s;
	delete [] solver_precon;
}




SlMatrix3x3 computeInertiaTensorPoints(double mass,
									   SlVector3 localVertices[4]){
  
  SlMatrix3x3 ret(0.0);
  SlMatrix3x3 identity;
  identity.setIdentity();
  for(auto i : range(4)){
	ret += 
	  (mass/4)*(dot(localVertices[i], localVertices[i])*identity -
				SlMatrixProduct(localVertices[i], localVertices[i]));
  }
  
  return ret;
}


void Tet::setInertiaTerms(SlVector3* pos,
						  SlVector3* vel,
						  double density,
						  btConvexHullShape* shape,
						  btRigidBody* body){
  
  SlVector3 centerOfMass(0.0);
  auto comVelocity = btVector3{0.0, 0.0, 0.0};
  for(auto vInd : vertices){
	centerOfMass += pos[vInd];
	comVelocity += btVector3{vel[vInd](0),
		vel[vInd](1),
		vel[vInd](2)};
  }
  centerOfMass /=4; //constant density over the tet
  comVelocity /= 4;
  
  SlVector3 localVertices[4]; //to align the tet with the principal axes and put it at the COM


  auto mass = 1.0/body->getInvMass();
  auto angularMomentum = SlVector3{0.0};
  for(auto i : range(4)){
	localVertices[i] = pos[vertices[i]] - centerOfMass;
	angularMomentum += (mass/4)*cross(localVertices[i], vel[vertices[i]]);
  }
  
  auto inertiaTensor = computeInertiaTensorPoints(mass,
												  localVertices);
  
  
  auto angularVelocity = inverse(inertiaTensor)*angularMomentum;
  


  SlMatrix3x3 eVecs;
  SlVector3 eVals;
  
  SlSymetricEigenDecomp(inertiaTensor,
						 eVals,
						 eVecs);

  //rotate the local vertices so that its aligned with the principal axes
  
  auto* shapePoints = shape->getUnscaledPoints();

  for(auto i  : range(4)){
	auto& v = localVertices[i];
	v = transmult(eVecs, v);
	
	//add the points to the shape
	shapePoints[i] = btVector3{v(0), v(1), v(2)}; 
  }
  assert(shape->getNumVertices() == 4);
  shape->recalcLocalAabb();

  //eVals are the principal inertia tensor components
  body->setMassProps(mass, btVector3{eVals(0), eVals(1), eVals(2)}); 
  
  auto rotMatrix = btMatrix3x3{eVecs(0,0), eVecs(0,1), eVecs(0,2),
							   eVecs(1,0), eVecs(1,1), eVecs(1,2),
							   eVecs(2,0), eVecs(2,1), eVecs(2,2)};

  body->setCenterOfMassTransform(btTransform{rotMatrix,
		btVector3{centerOfMass(0),
		  centerOfMass(1),
		  centerOfMass(2)}});
  
 
  body->setLinearVelocity(comVelocity);
  body->setAngularVelocity({angularVelocity(0),
		angularVelocity(1),
		angularVelocity(2)});
  


#define EXTRA_INERTIA_CHECKS 1
#if EXTRA_INERTIA_CHECKS
 
  //should double check that this transform is correct and gives us the world positions of vertices
  auto trans = body->getCenterOfMassTransform();
  for(auto i : range(4)){
	auto worldPos = trans(shape->getUnscaledPoints()[i]);
	assert(fabs(worldPos.x() - pos[vertices[i]](0)) < 1e-7);
	assert(fabs(worldPos.y() - pos[vertices[i]](1)) < 1e-7);
	assert(fabs(worldPos.z() - pos[vertices[i]](2)) < 1e-7);

  }

  //double check that the new inertia tensor is diagonal
  auto newInertiaTensor = computeInertiaTensorPoints(mass,
													 localVertices);

  std::cout << "new inertia tensor: \n" << newInertiaTensor << std::endl;

  //std::cout << "new det: " << determinant(newInertiaTensor) << std::endl;

  //newInertiaTensor = scale*newInertiaTensor;
  //std::cout << "scaled new inertia tensor: \n" << newInertiaTensor << std::endl;

  std::cout << "evecs: " << eVecs << std::endl;
  std::cout << "evals (scaled): " << eVals << std::endl;
  std::cout << "mass: " << mass << std::endl;

  std::cout << "RT before R " <<
	transpose(eVecs)*(inertiaTensor)*eVecs << std::endl;

  //double check that the new center of mass is close to 0
  SlVector3 localCOM = localVertices[0] + localVertices[1] + localVertices[2] + localVertices[3];
  assert(fabs(localCOM.x()) < 1e-7);
  assert(fabs(localCOM.y()) < 1e-7);
  assert(fabs(localCOM.z()) < 1e-7);


  assert(fabs(newInertiaTensor(0,1)) < 1e-7);
  assert(fabs(newInertiaTensor(0,2)) < 1e-7);
  assert(fabs(newInertiaTensor(1,0)) < 1e-7);
  assert(fabs(newInertiaTensor(1,2)) < 1e-7);
  assert(fabs(newInertiaTensor(2,0)) < 1e-7);
  assert(fabs(newInertiaTensor(2,1)) < 1e-7);
  assert(fabs(newInertiaTensor(0,0) - eVals(0)) < 1e-7);
  assert(fabs(newInertiaTensor(1,1) - eVals(1)) < 1e-7);
  assert(fabs(newInertiaTensor(2,2) - eVals(2)) < 1e-7);
		 
  std::cout << "inertia checks pass" << std::endl;
#endif
  
}



void FemObject::setupBulletMesh(){
  /*
  bulletVertexArray = 
    std::unique_ptr<btTriangleIndexVertexArray>{new btTriangleIndexVertexArray{ntris,
									       reinterpret_cast<int*>(tris),
									       sizeof(SlInt3),
									       nv,
									       reinterpret_cast<double*>(pos),
									       sizeof(SlVector3)}};
  bulletMeshShape = 
    std::unique_ptr<btGImpactMeshShape>{new btGImpactMeshShape{bulletVertexArray.get()}};

  bulletMeshShape->updateBound(); //update the AABBs
  */

  //create a rigid body per tet

  bulletTetShapes.reserve(ntets);
  bulletTetShapes.resize(0);
  for(auto i : range(ntets)){
	bulletTetShapes.push_back(std::unique_ptr<btConvexHullShape>{new btConvexHullShape{}});
	//add 4 points
	for(auto j : range(4)){
	  bulletTetShapes.back()->addPoint(btVector3{0,0,0}, false); //don't recalc AABBs
	}
	
	bulletMotionStates.
	  push_back(std::unique_ptr<btDefaultMotionState>{
		  new btDefaultMotionState{btTransform{btQuaternion::getIdentity(),
				btVector3(0,0,0)}}});
	
	bulletTetBodies.push_back(std::unique_ptr<btRigidBody>{
		new btRigidBody{tetMass[i], 
			bulletMotionStates.back().get(),
			bulletTetShapes.back().get()}});

	tets[i].setInertiaTerms(pos, vel, density[i], 
							bulletTetShapes.back().get(),
							bulletTetBodies.back().get());
	
	bulletWorld->addRigidBody(bulletTetBodies.back().get());
	

  }

}

void FemObject::updateBulletShapes(){
  for(auto i : range(ntets)){
	auto& t = tets[i];
	//compute moment of inertia, and set linear/angular velocities
	t.setInertiaTerms(pos, vel, density[i],
					  bulletTetShapes[i].get(),
					  bulletTetBodies[i].get());

  }

}


void FemObject::setupRigidProperties(){

  /*  btVector3 localInertia(0,0,0);
  bulletMeshShape->updateBound();
  bulletMeshShape->calculateLocalInertia(totalMass, localInertia);
  
  bulletMotionStates.
	push_back(std::unique_ptr<btDefaultMotionState>{
		new btDefaultMotionState{btTransform{btQuaternion::getIdentity(),
			  btVector3{0,0,0}}}};
	  
  */
}

//do simple mass weighted averaging of the vertex positions
void FemObject::stitchTets(){
  
  std::fill(pos, pos + nv, SlVector3{0.0});
  
  for(auto tInd : range(ntets)){
	auto& t = tets[tInd];
	auto& trans = bulletTetBodies[tInd]->getCenterOfMassTransform();
	for(auto i : range(4)){
	  auto worldPos = trans(bulletTetShapes[tInd]->getUnscaledPoints()[i]);
	  pos[t.vertices[i]] += (tetMass[tInd]/4)*SlVector3{worldPos.x(), worldPos.y(), worldPos.z()};
	}
  }
  for(auto pInd : range(nv)){
	pos[pInd] /= mass[pInd];
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


//this is incorrect, but may be worth pursuing later on

#if 0

double computeDiagonalTerm(SlVector3 p[4],
						   int i){
  
  return p[0](i)*(p[0](i) +
				  p[1](i) +
				  p[2](i) +
				  p[3](i)) +
	p[1](i)*(p[1](i) +
			 p[2](i) +
			 p[3](i)) +
	p[2](i)*(p[2](i) + p[3](i)) +
	p[3](i)*p[3](i);
}

double computeOffdiagTerm(SlVector3 p[4],
						  int i, int j){

  double ret = 0;
  for(auto k : range(4)){
	for(auto l : range(4)){	  
	  ret += ((k == l) ? 2*p[k](i)*p[l](j) :
			  p[k](i)*p[l](j));
	}
  }
  /*ret += p[0](i)*p[0](j) +
	p[1](i)*p[1](j) +
	p[2](i)*p[2](j) +
	p[3](i)*p[3](j);*/
  return ret;
}




SlMatrix3x3 computeInertiaTensor(double mass,
								 SlVector3 localVertices[4]){

  //use the formulas from here: http://www.thescipub.com/abstract/?doi=jmssp.2005.8.11
  //to compute the inertia tensor
  auto mass_60 = mass/10;
  auto mass_120 = mass/20;
  auto diagX = computeDiagonalTerm(localVertices, 0);
  auto diagY = computeDiagonalTerm(localVertices, 1);
  auto diagZ = computeDiagonalTerm(localVertices, 2);
  double a = mass_60*(diagY + diagZ);
  double b = mass_60*(diagX + diagZ);
  double c = mass_60*(diagX + diagY);
  double aPrime = mass_120*computeOffdiagTerm(localVertices, 1, 2);
  double bPrime = mass_120*computeOffdiagTerm(localVertices, 0, 2);
  double cPrime = mass_120*computeOffdiagTerm(localVertices, 0, 1);
  //there's probably an easier way to compute this
  
  return SlMatrix3x3 (a, -bPrime, -cPrime,
					  -bPrime, b, -aPrime,
					  -cPrime, -aPrime, c);
}
#endif
