#include "grip.H"

#include "fem.H"

#include <cstring>

PinGrip::PinGrip(const unsigned int *in_vi, unsigned int in_nv)
:nv(in_nv)
{
	vi = new unsigned int[nv];
	memcpy(vi, in_vi, nv*sizeof(unsigned int));
}

PinGrip::~PinGrip(){
	delete [] vi;
}

void PinGrip::flagGrips(bool *flags){
	for(unsigned int i = 0; i < nv; i++){
		flags[vi[i]] = true;
	}
}

void PinGrip::applyGrips(SlVector3 *dv, double dt){
	for(unsigned int i = 0; i < nv; i++){
		dv[vi[i]] = SlVector3(0.0);
	}
}

AbsPosGrip::AbsPosGrip(const FemObject *in_obj, unsigned int in_ngrips)
:obj(in_obj),ngrips(in_ngrips) {
	vi = new unsigned int[ngrips];
	vpos = new SlVector3[ngrips];	
}

AbsPosGrip::AbsPosGrip(const AbsPosGrip &g)
	:FemGrip(),obj(g.obj),ngrips(g.ngrips){
	vi = new unsigned int [ngrips];
	vpos = new SlVector3[ngrips];

	memcpy(vi, g.vi, sizeof(unsigned int)*ngrips);
	memcpy(vpos, g.vpos, sizeof(SlVector3)*ngrips);	
}

AbsPosGrip::~AbsPosGrip(){
	delete [] vi;
	delete [] vpos;
}

void AbsPosGrip::flagGrips(bool *flags){
	for(unsigned int i = 0; i < ngrips; i++){
		flags[vi[i]] = true;
	}
}

void AbsPosGrip::applyGrips(SlVector3 *dv, double dt){
	double dtInv = 1.0 / dt;
	for(unsigned int i = 0; i < ngrips; i++){
		unsigned int vindex = vi[i];
		SlVector3 oldPos = obj->pos[vindex];
		SlVector3 oldVel = obj->vel[vindex];
		dv[vindex] = dtInv*(vpos[i] - oldPos) - oldVel;
	}
}

KineGrip::KineGrip(const FemObject *in_obj, const unsigned int *in_vi, unsigned int in_nv, SlVector3 **in_vpos, unsigned int in_nframes) 
	:FemGrip(),obj(in_obj),nv(in_nv),nframes(in_nframes),theTime(0.0) {
	vi = new unsigned int[nv];
	memcpy(vi, in_vi, nv*sizeof(unsigned int));
	memory = new SlVector3[nframes*nv];
	vpos = new SlVector3*[nv];
	for (unsigned int i=0; i<nv; i++) {
		vpos[i] = memory + i*nframes;
		for (unsigned int j=0; j<nframes; j++) {
			vpos[i][j] = in_vpos[i][j];
		}
	}
}

KineGrip::KineGrip(const KineGrip &g)
	:FemGrip(),obj(g.obj),nv(g.nv),nframes(g.nframes),theTime(g.theTime){
	vi = new unsigned int [nv];
	memory = new SlVector3[nv*nframes];
	vpos = new SlVector3*[nv];

	memcpy(vi, g.vi, sizeof(unsigned int)*nv);
	for (unsigned int i=0; i<nv; i++) {
		vpos[i] = memory + i*nframes;
		for (unsigned int j=0; j<nframes; j++) {
			vpos[i][j] = g.vpos[i][j];
		}
	}
}

KineGrip::~KineGrip(){
	delete [] vi;
	delete [] memory;
	delete [] vpos;
}

void KineGrip::flagGrips(bool *flags){
	for(unsigned int i = 0; i < nv; i++){
		flags[vi[i]] = true;
	}
}

void KineGrip::applyGrips(SlVector3 *dv, double dt){
	theTime += dt;
	int frame = (int)floor(theTime / 30.0);
	double w = theTime - frame*30.0;

	double dtInv = 1.0 / dt;
	for(unsigned int i = 0; i < nv; i++){
		unsigned int vindex = vi[i];
		SlVector3 pos = (1.0-w)*vpos[i][frame] + w*vpos[i][frame+1];
		SlVector3 oldPos = obj->pos[vindex];
		SlVector3 oldVel = obj->vel[vindex];
		dv[vindex] = dtInv*(pos - oldPos) - oldVel;
	}
}

RotateGrip::RotateGrip(const FemObject *in_obj, const unsigned int *in_vi, unsigned int in_nv, double in_speed)
	:obj(in_obj),nv(in_nv),speed(in_speed),theTime(0.0){
	vi = new unsigned int[nv];
	opos = new SlVector3[nv];	
	memcpy(vi, in_vi, nv*sizeof(unsigned int));
	for (unsigned int i=0; i<nv; i++) {
		unsigned int j = vi[i];
		opos[i] = obj->pos[j];
	}
}

RotateGrip::~RotateGrip(){
	delete [] vi;
	delete [] opos;
}

void RotateGrip::flagGrips(bool *flags){
	for(unsigned int i = 0; i < nv; i++){
		flags[vi[i]] = true;
	}
}

void RotateGrip::applyGrips(SlVector3 *dv, double dt){
	double dtInv = 1.0 / dt;
	theTime += dt;
	double angVel = speed*2*M_PI;
	double theta = theTime*angVel;
	for(unsigned int i = 0; i < nv; i++){
		unsigned int j = vi[i];
		SlVector3 newPos;
		newPos[0] = opos[i][0];
		newPos[1] = cos(theta)*opos[i][1] - sin(theta)*opos[i][2];
		newPos[2] = sin(theta)*opos[i][1] + cos(theta)*opos[i][2];
		SlVector3 &oldPos = obj->pos[j];
		SlVector3 &oldVel = obj->vel[j];
		dv[j] = dtInv*(newPos - oldPos) - oldVel;
	}
}

SpiralGrip::SpiralGrip(const FemObject *in_obj, const unsigned int *in_vi, unsigned int in_nv, double in_speed)
	:obj(in_obj),nv(in_nv),speed(in_speed),theTime(0.0){
	vi = new unsigned int[nv];
	opos = new SlVector3[nv];	
	memcpy(vi, in_vi, nv*sizeof(unsigned int));
	for (unsigned int i=0; i<nv; i++) {
		unsigned int j = vi[i];
		opos[i] = obj->pos[j];
	}
}

SpiralGrip::~SpiralGrip(){
	delete [] vi;
	delete [] opos;
}

void SpiralGrip::flagGrips(bool *flags){
	if (theTime > 45) return;  // NEW
		for(unsigned int i = 0; i < nv; i++){
			flags[vi[i]] = true;
		}
}

void SpiralGrip::applyGrips(SlVector3 *dv, double dt){
	double l = 10.0;
	double dtInv = 1.0 / dt;
	theTime += dt;
	double theta = theTime*speed;
	//if (theta > 1) theta = 1;
	if (theTime > 45) return;
	for(unsigned int i = 0; i < nv; i++){
		unsigned int j = vi[i];
		SlVector3 newPos;
		if (theta <= 0) {
			newPos = opos[i];
		}	else if (theta <= 1) {
			if (opos[i][0] < l/2.0) {
				newPos[0] = (l/2.0)*theta;
			} else {
				newPos[0] = l - (l/2.0)*theta;
			}
			if (opos[i][0] < l/2.0) {
				newPos[1] = opos[i][1] + 0.2*theta;
			} else {
				newPos[1] = opos[i][1] - 0.2*theta;
			}				
			newPos[2] = theta * (l/(2*M_PI)-2*opos[i][2]) + opos[i][2];
		} else {
			if (opos[i][0] < l/2.0) {
				newPos[1] = opos[i][1] + 0.2*theta;
			} else {
				newPos[1] = opos[i][1] - 0.2*theta;
			}				
			theta -= 1;	
			if (opos[i][0] > l/2.0) {
				newPos[0] = l/2.0 + (1-theta/30.0)*sin(-theta)*(l/(2*M_PI))-sin(-theta)*opos[i][2];// - opos[i][2]); // was / 30
				newPos[2] = (1-theta/30.0)*cos(-theta)*(l/(2*M_PI)) - cos(-theta)*opos[i][2];
			} else {
				newPos[0] = l/2.0 + (1-theta/30.0)*sin(theta)*(l/(2*M_PI))-sin(theta)*opos[i][2];// - opos[i][2]);
				newPos[2] = (1-theta/30.0)*cos(theta)*(l/(2*M_PI)) - cos(theta)*opos[i][2];
			}
			theta += 1;
		}
		//std::cout<<opos[i]<<"->"<<newPos<<std::endl;
		
		SlVector3 &oldPos = obj->pos[j];
		SlVector3 &oldVel = obj->vel[j];
		dv[j] = dtInv*(newPos - oldPos) - oldVel;
		//if (theta > 1) std::cout<<oldPos<<"->"<<newPos<<std::endl;
	}
}

