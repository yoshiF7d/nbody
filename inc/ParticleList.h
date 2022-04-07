#ifndef PARTICLE_LIST_H
#define PARTICLE_LIST_H
#include <Particle.h>
#include <iostream>
#include <vector>

struct ParticleList{
	std::vector<Particle> list;
	char *buffer = nullptr;
	size_t bufferSize = 0;
	
	ParticleList(int size);
	~ParticleList();
	Particle &operator[] (int n);
	
	void evaluate();
	void evaluate(int index);
	void evaluate(int size, int index);
	void evaluate(ParticleList &particles);
	void evaluate(ParticleList &particles,int index);
	void evaluate(ParticleList &particles,int size,int index);
	
	void evaluate2();
	void evaluate2(ParticleList &particles,int size,int index);
	
	void boundary(double range[][2]);
	void init(double range[][2]);
		
	void setBuffer();
	void setBuffer(int size);
	
	void pack();
	void pack(char *data);
	void pack(int size);
	void pack(char *data, int size);
	void pack(char *data, int size, int offset);
	
	void unpack();
	void unpack(char *data);
	void unpack(int size);
	void unpack(char *data, int size);
	void unpack(char *data, int size, int offset);
};

#endif