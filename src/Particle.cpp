#include <Particle.h>
Particle::Particle(
	int id,
	double mass, 
	double charge, 
	std::array<double,6> position,
	std::array<double,6> momentum
){
	this->step = 1;
	this->pair = -1;
	this->force = {0,0,0,0,0,0};
	this->yank = {0,0,0,0,0,0};
	this->set(id,mass,charge,position,momentum);
}

void Particle::set(
	int id,
	double mass,
	double charge,
	std::array<double,6> position,
	std::array<double,6> momentum
){
	this->id = id;
	this->mass = mass;
	this->charge = charge;
	this->position = position;
	this->momentum = momentum;
}

double Particle::interact(Particle& p){
	double x = position[0] - p.position[0];
	double y = position[1] - p.position[1];
	double px = momentum[0] - p.momentum[0];
	double py = momentum[1] - p.momentum[1];
	double r = sqrt(x*x+y*y+epsilon);
		double r3 = r*r*r;
		double r5 = r3*r*r;
		force[0] += x*((charge)*(p.charge))/r3;
		force[1] += y*((charge)*(p.charge))/r3;
		yank[0] += px*((charge)*(p.charge))/r3;
		yank[0] -= 3*(px*x+py*y)*x*((charge)*(p.charge))/r5;
		yank[1] += py*((charge)*(p.charge))/r3;
		yank[1] -= 3*(px*x+py*y)*y*((charge)*(p.charge))/r5;
		return (charge)*(p.charge)/abs(r);
}

void Particle::evaluate(Particle& p){
	double x = position[3] - p.position[3];
	double y = position[4] - p.position[4];
	double px = momentum[3] - p.momentum[3];
	double py = momentum[4] - p.momentum[4];
	double r = sqrt(x*x+y*y+epsilon);
	double r3 = r*r*r;
	double r5 = r3*r*r;
	force[3] += x*((charge)*(p.charge))/r3;
	force[4] += y*((charge)*(p.charge))/r3;
	yank[3] += px*((charge)*(p.charge))/r3;
	yank[3] -= 3*(px*x+py*y)*x*((charge)*(p.charge))/r5;
	yank[4] += py*((charge)*(p.charge))/r3;
	yank[4] -= 3*(px*x+py*y)*y*((charge)*(p.charge))/r5;
}

void Particle::predict(double dt){
	momentum[3] = 0.5*dt*dt*yank[0] + dt*force[0] + momentum[0];
	momentum[4] = 0.5*dt*dt*yank[1] + dt*force[1] + momentum[1];
	position[3] = (0.16666667*dt*dt*dt*yank[0] + 0.5*dt*dt*force[0] 
					+ dt*momentum[0]) / mass + position[0];
	position[4] = (0.16666667*dt*dt*dt*yank[1] + 0.5*dt*dt*force[1] 
					+ dt*momentum[1]) / mass + position[1];
}

void Particle::correct(double dt){
	momentum[3] = momentum[0] + 0.5*dt*(force[3]+force[0])
						- 0.083333333*dt*dt*(yank[3]-yank[0]);
	momentum[4] = momentum[1] + 0.5*dt*(force[4]+force[1])
						- 0.083333333*dt*dt*(yank[4]-yank[1]);
	position[3] = position[0] + (
							0.5*dt*(momentum[3]+momentum[0]) -
							0.083333333*dt*dt*(force[3]-force[0])
						)/mass;
	position[4] = position[1] + (
							0.5*dt*(momentum[4]+momentum[1]) -
							0.083333333*dt*dt*(force[4]-force[1])
						)/mass;
}

void Particle::push(){
	momentum[0] = momentum[3];
	momentum[1] = momentum[4];
	position[0] = position[3];
	position[1] = position[4];
}

double Particle::error(){
	return (force[3] - force[0])*(force[3] - force[0]) + 
	(force[4] - force[1])*(force[4] - force[1]);
}

void Particle::print(std::ostream &out){
	out << id << "\t"
	<< step << "\t"
	<< mass << "\t"
	<< charge << "\t"
	<< position[0] << "\t"
	<< position[1] << "\t"
	<< momentum[0] << "\t"
	<< momentum[1] << "\t"
	<< force[0] << "\t"
	<< force[1] << "\t\n";
}

void Particle::serialize(std::vector<Particle> &array, int length, char **data, size_t *size){
	*size = 0;
	*size += sizeof(int); // id
	*size += sizeof(double); // mass
	*size += sizeof(double); // charge
	*size += 2*sizeof(double); // position
	*size += 2*sizeof(double); // momentum
	*size += 2*sizeof(double); // force
	*size *= length;
	char *head = *data = (char*)malloc(*size);
	for(int i=0;i<length;i++){
		memcpy(head,&(array[i].id),sizeof(int)); head += sizeof(int);
		memcpy(head,&(array[i].mass),sizeof(double)); head += sizeof(double);
		memcpy(head,&(array[i].charge),sizeof(double)); head += sizeof(double);
		memcpy(head,&(array[i].position[0]),sizeof(double)); head += sizeof(double);
		memcpy(head,&(array[i].position[1]),sizeof(double)); head += sizeof(double);
		memcpy(head,&(array[i].momentum[0]),sizeof(double)); head += sizeof(double);
		memcpy(head,&(array[i].momentum[1]),sizeof(double)); head += sizeof(double);
		memcpy(head,&(array[i].force[0]),sizeof(double)); head += sizeof(double);
		memcpy(head,&(array[i].force[1]),sizeof(double)); head += sizeof(double);
	}
}
std::vector<Particle> Particle::deserialize(char *data, size_t size){
	size_t s = 0;
	s += sizeof(int); // id
	s += sizeof(double); // mass
	s += sizeof(double); // charge
	s += 2*sizeof(double); // position
	s += 2*sizeof(double); // momentum
	s += 2*sizeof(double); // force
	std::vector<Particle> array(size/s);
	char *head = data;
	for(int i=0;i<array.size();i++){
		std::memcpy(&(array[i].id),head,sizeof(int)); head += sizeof(int);
		std::memcpy(&(array[i].mass),head,sizeof(double)); head += sizeof(double);
		std::memcpy(&(array[i].charge),head,sizeof(double)); head += sizeof(double);
		std::memcpy(&(array[i].position[0]),head,sizeof(double)); head += sizeof(double);
		std::memcpy(&(array[i].position[1]),head,sizeof(double)); head += sizeof(double);
		std::memcpy(&(array[i].momentum[0]),head,sizeof(double)); head += sizeof(double);
		std::memcpy(&(array[i].momentum[1]),head,sizeof(double)); head += sizeof(double);
		std::memcpy(&(array[i].force[0]),head,sizeof(double)); head += sizeof(double);
		std::memcpy(&(array[i].force[1]),head,sizeof(double)); head += sizeof(double);
	}
	free(data);
	return array;
}

void Particle::deserialize_fill(std::vector<Particle> &array, int *index, char *data, size_t size){
	size_t s = 0;
	s += sizeof(int); // id
	s += sizeof(double); // mass
	s += sizeof(double); // charge
	s += 2*sizeof(double); // position
	s += 2*sizeof(double); // momentum
	s += 2*sizeof(double); // force
	char *head = data;
	//std::cout << "deserialize_fill : start : " << *index << std::endl;
	//std::cout << "deserialize_fill : end : " << *index + size/s << std::endl;
	for(int i=*index;i<*index+size/s;i++){
		//array[i].print();
		std::memcpy(&(array[i].id),head,sizeof(int)); head += sizeof(int);
		std::memcpy(&(array[i].mass),head,sizeof(double)); head += sizeof(double);
		std::memcpy(&(array[i].charge),head,sizeof(double)); head += sizeof(double);
		std::memcpy(&(array[i].position[0]),head,sizeof(double)); head += sizeof(double);
		std::memcpy(&(array[i].position[1]),head,sizeof(double)); head += sizeof(double);
		std::memcpy(&(array[i].momentum[0]),head,sizeof(double)); head += sizeof(double);
		std::memcpy(&(array[i].momentum[1]),head,sizeof(double)); head += sizeof(double);
		std::memcpy(&(array[i].force[0]),head,sizeof(double)); head += sizeof(double);
		std::memcpy(&(array[i].force[1]),head,sizeof(double)); head += sizeof(double);
		//array[i].print();
		//std::cout << "===\n";
	}
	*index += size/s;
	free(data);
}