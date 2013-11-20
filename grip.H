#ifndef GRIP_H_
#define GRIP_H_

#include <slVector.H>

class FemObject;

class FemGrip {
public:
 	virtual void flagGrips(bool *flags) = 0;
	virtual void applyGrips(SlVector3 *dv, double dt) = 0;
};

class PinGrip : public FemGrip {
private:
	unsigned int *vi;
	unsigned int nv;
public:
	PinGrip(const unsigned int *vi, unsigned int nv);
	virtual ~PinGrip();
	virtual void flagGrips(bool *flags);
	virtual void applyGrips(SlVector3 *dv, double dt);
};

class AbsPosGrip : public FemGrip {
public:
	unsigned int *vi;
	SlVector3 *vpos;

	const FemObject *obj;
	const unsigned int ngrips;

	AbsPosGrip(const FemObject *obj, unsigned int ngrips);
	AbsPosGrip(const AbsPosGrip &g);
	virtual ~AbsPosGrip();

	virtual void flagGrips(bool *flags);
	virtual void applyGrips(SlVector3 *dv, double dt);
};

class KineGrip : public FemGrip {
private:
	unsigned int *vi;
	SlVector3 **vpos;
	SlVector3 *memory;
	const FemObject *obj;
	const	unsigned int nv;
	const unsigned int nframes;
	double theTime;

public:

	KineGrip(const FemObject *obj, const unsigned int *vi, unsigned int nv, SlVector3 **vpos, unsigned int nframes);
	KineGrip(const KineGrip &g);
	virtual ~KineGrip();

	virtual void flagGrips(bool *flags);
	virtual void applyGrips(SlVector3 *dv, double dt);
};

class RotateGrip :public FemGrip {
private:
	unsigned int *vi;
	const FemObject *obj;
	SlVector3 *opos;
	unsigned int nv;
	double speed; // in rotations per second
	double theTime;

public:
	RotateGrip(const FemObject *obj, const unsigned int *vi, unsigned int nv, double speed); 
	virtual ~RotateGrip();

	virtual void flagGrips(bool *flags);
	virtual void applyGrips(SlVector3 *dv, double dt);
};

class SpiralGrip :public FemGrip {
private:
	unsigned int *vi;
	const FemObject *obj;
	SlVector3 *opos;
	unsigned int nv;
	double speed; // in rotations per second
	double theTime;

public:
	SpiralGrip(const FemObject *obj, const unsigned int *vi, unsigned int nv, double speed); 
	virtual ~SpiralGrip();

	virtual void flagGrips(bool *flags);
	virtual void applyGrips(SlVector3 *dv, double dt);
};
	

#endif
