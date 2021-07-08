#ifndef PARTICLE_H
#define PARTICLE_H
#include <iostream>
#include <array>
#include <vector>
#include <cstring>
#include <cmath>
const double epsilon=0.1;
const double rmin = 0.5;

struct Particle{
	int id;
	int pair;
	int step;
	double mass;
	double charge;
	position.()
	std::array<double,6> position;
	std::array<double,6> momentum;
	std::array<double,6> force;
	std::array<double,6> yank;
	Particle(
		int id = 0,
		double mass = 0, 
		double charge = 0,
		std::array<double,6> position = {0,0,0,0,0,0},
		std::array<double,6> momentum = {0,0,0,0,0,0}
	);
	void set(
		int id,
		double mass,
		double charge,
		std::array<double,6> position = {0,0,0,0,0,0},
		std::array<double,6> momentum = {0,0,0,0,0,0}
	);
	double error();
	double interact(Particle&);
	void evaluate(Particle&);
	void predict(double);
	void correct(double);
	void push();
	void print(std::ostream &out=std::cout);
	static void serialize(std::vector<Particle>&,int,char**,size_t*);
	static std::vector<Particle> deserialize(char*,size_t);
	static void deserialize_fill(std::vector<Particle>&,int*,char*,size_t);
};
#endif
