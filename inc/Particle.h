#ifndef PARTICLE_H
#define PARTICLE_H
#include <iostream>
#include <fstream>
#include <Constants.h>

#define X position[0]
#define Y position[1]
#define Z position[2]

#define PX momentum[0]
#define PY momentum[1]
#define PZ momentum[2]

#define FX force[0]
#define FY force[1]
#define FZ force[2]

#define YX yank[0]
#define YY yank[1]
#define YZ yank[2]

struct Particle{
	int id;
	int step;
	double mass;
	double charge;

	double position[DIMENSION][ORDER];
	double momentum[DIMENSION][ORDER];
	double force[DIMENSION][ORDER];
	double yank[DIMENSION][ORDER];
	
	/*position[0][0] : x0*/
	/*position[1][0] : x1*/
	/*position[2][0] : x2*/

	/*position[0][1] : y0*/
	/*position[1][1] : y1*/
	/*position[2][1] : y2*/

	/*x0 is for output x1,x2,... are for calculation*/
	
	double error(int i, int j);
	void evaluate(Particle& particle,int k=0);
	void evaluate2(Particle& particle,int k=0);
	void predict(double);
	void correct(double);
	void move(int i, int j);

	void reflect(double range[][2]);
	void periodic(double range[][2]);
	
	static void header(std::ostream &out=std::cout);
	void print(std::ostream &out=std::cout);
	void write(std::ofstream &ofs);
	
	static void init(std::vector<Particle> &array, int rmax);
	static void pack(std::vector<Particle> &array, int length, int offset, char *data);
	static void unpack(std::vector<Particle> &array, int length, int offset, char *data);
};
#endif
