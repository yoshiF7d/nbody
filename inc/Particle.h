#ifndef PARTICLE_H
#define PARTICLE_H
#include <iostream>

class Particle{
  public:
	int id;
	double m;
	double q;
	double x;
	double y;
	double px;
	double py;
	double fx;
	double fy;
	
	Particle(int id=0,double m=0, double q=0, double x=0, double y=0, double px=0, double py=0);
	void print(std::ostream &out=std::cout);
};
#endif
