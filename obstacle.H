#ifndef OBSTACLE_H_
#define OBSTACLE_H_

#include <slVector.H>

class FemObject;

class FemObstacle {
public:
	virtual void step(double dt) = 0;
	virtual void collide(FemObject &o, double dt) = 0;
};

class SphereObstacle : public FemObstacle {
private:
	SlVector3 center, vel;
	double radius;
public:
	SphereObstacle(const SlVector3 &center, const SlVector3 &vel, double radius):center(center), vel(vel), radius(radius){};
	void step(double dt);
	void collide(FemObject &o, double dt);
};

#endif
