#include <Particle.h>

Particle::Particle(int id, double m, double q, double x, double y, double px, double py){
	this->id = id;
	this->m = m;
	this->q = q;
	this->x = x;
	this->y = y;
	this->px = px;
	this->py = py;
	this->fx = 0;
	this->fy = 0;
}

void Particle::print(std::ostream &out){
	out << id << "\t"
	<< m << "\t"
	<< q << "\t"
	<< x << "\t"
	<< y << "\t"
	<< px << "\t"
	<< py << "\t"
	<< fx << "\t"
	<< fy << "\t\n";
}
