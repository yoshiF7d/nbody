#include <Particle.h>
#include <Constants.h>
#include <iomanip>
#include <random>

void Particle::evaluate(Particle& particle, int k){
	double x = X[k] - particle.X[k];
	double y = Y[k] - particle.Y[k];
	double px = PX[k] - particle.PX[k];
	double py = PY[k] - particle.PY[k];
	double r = sqrt(x*x+y*y+EPSILON);
	double r3 = r*r*r;
	double r5 = r3*r*r;
	double cc = charge * particle.charge;
	FX[k] += x*cc/r3;
	FY[k] += y*cc/r3;
	/*
	YX[k] += px*cc/r3;
	YX[k] -= 3*(px*x+py*y)*x*cc/r5;
	YY[k] += py*cc/r3;
	YY[k] -= 3*(px*x+py*y)*y*cc/r5;
	*/
}

void Particle::evaluate2(Particle& particle, int k){
	double x = X[k] - particle.X[k];
	double y = Y[k] - particle.Y[k];
	double px = PX[k] - particle.PX[k];
	double py = PY[k] - particle.PY[k];
	double r = sqrt(x*x+y*y+EPSILON);
	double r3 = r*r*r;
	double r5 = r3*r*r;
	double cc = charge * particle.charge;
	FX[k] += x*cc/r3;
	FY[k] += y*cc/r3;
	
	YX[k] += px*cc/r3;
	YX[k] -= 3*(px*x+py*y)*x*cc/r5;
	YY[k] += py*cc/r3;
	YY[k] -= 3*(px*x+py*y)*y*cc/r5;
}

/*taylor expansion*/
/*1/2,1*/
/*1/6,1/2,1*/
void Particle::predict(double dt){
	PX[1] = 0.5*dt*dt*YX[0] + dt*FX[0] + PX[0];
	PY[1] = 0.5*dt*dt*YY[0] + dt*FY[0] + PY[0];
	X[1] = (0.16666667*dt*dt*dt*YX[0] + 0.5*dt*dt*FX[0] + dt*PX[0]) / mass + X[0];
	Y[1] = (0.16666667*dt*dt*dt*YY[0] + 0.5*dt*dt*FY[0] + dt*PY[0]) / mass + Y[0];
}

/*1/2,-1/12*/
/*1/2,-1/12*/
void Particle::correct(double dt){
	PX[1] = PX[0] + 0.5*dt*(FX[1]+FX[0])- 0.083333333*dt*dt*(YX[1]-YX[0]);
	PY[1] = PY[0] + 0.5*dt*(FY[1]+FY[0])- 0.083333333*dt*dt*(YY[1]-YY[0]);
	X[1] = X[0] + (0.5*dt*(PX[1]+PX[0]) -0.083333333*dt*dt*(FX[1]-FX[0]))/mass;
	Y[1] = Y[0] + (0.5*dt*(PY[1]+PY[0]) -0.083333333*dt*dt*(FY[1]-FY[0]))/mass;
}

void Particle::move(int i, int j){
	PX[j] = PX[i];
	PY[j] = PY[i];
	X[j] = X[i];
	Y[j] = Y[i];
}

double Particle::error(int i, int j){
	return (PX[i] - PX[j])*(PX[i] - PX[j]) + (PY[i] - PY[j])*(PY[i] - PY[j]);
}

void Particle::reflect(double range[][2]){
	if(X[0] > range[0][1]){X[0] -= 2*(X[0]-range[0][1]); PX[0] *= -1;}
	if(X[0] < range[0][0]){X[0] += 2*(range[0][0]-X[0]); PX[0] *= -1;}
	if(Y[0] > range[1][1]){Y[0] -= 2*(Y[0]-range[1][1]); PY[0] *= -1;}
	if(Y[0] < range[1][0]){Y[0] += 2*(range[1][0]-Y[0]); PY[0] *= -1;}
}
void Particle::periodic(double range[][2]){
	if(X[0] > range[0][1]){X[0] -= range[0][1] - range[0][0];}
	if(X[0] < range[0][0]){X[0] += range[0][1] - range[0][0];}
	if(Y[0] > range[1][1]){Y[0] -= range[1][1] - range[1][0];}
	if(Y[0] < range[1][0]){Y[0] += range[1][1] - range[1][0];}
}

void Particle::header(std::ostream &out){
	out << std::setw(6) << "id"
	<< std::setw(6) << "step"
	<< std::setw(10) << "mass" 
	<< std::setw(10) << "charge"
	<< std::setw(10) << "x"
	<< std::setw(10) << "y"
	<< std::setw(10) << "px"
	<< std::setw(10) << "py"
	<< std::setw(10) << "fx"
	<< std::setw(10) << "fy" << "\n";
}

void Particle::print(std::ostream &out){
	out << std::setw(6) << id
	<< std::setw(6)  << step
	<< std::setw(10) << std::scientific << std::setprecision(2) << mass
	<< std::setw(10) << std::scientific << std::setprecision(2) << charge 
	<< std::setw(10) << std::scientific << std::setprecision(2) << X[0]
	<< std::setw(10) << std::scientific << std::setprecision(2) << Y[0]
	<< std::setw(10) << std::scientific << std::setprecision(2) << PX[0]
	<< std::setw(10) << std::scientific << std::setprecision(2) << PY[0]
	<< std::setw(10) << std::scientific << std::setprecision(2) << FX[0]
	<< std::setw(10) << std::scientific << std::setprecision(2) << FY[0] << "\n";
}

void Particle::write(std::ofstream &ofs){
	ofs.write((char*)&id,sizeof(int));
	ofs.write((char*)&mass,sizeof(double));
	ofs.write((char*)&charge,sizeof(double));
	ofs.write((char*)&X[0],sizeof(double));
	ofs.write((char*)&Y[0],sizeof(double));
	ofs.write((char*)&PX[0],sizeof(double));
	ofs.write((char*)&PY[0],sizeof(double));
	ofs.write((char*)&FX[0],sizeof(double));
	ofs.write((char*)&FY[0],sizeof(double));
}

void Particle::pack(std::vector<Particle> &array, int length, int offset, char *data){
	char *head = data;
	for(int i=offset;i<offset+length;i++){
		memcpy(head,&(array[i].id),sizeof(int)); head += sizeof(int);
		memcpy(head,&(array[i].mass),sizeof(double)); head += sizeof(double);
		memcpy(head,&(array[i].charge),sizeof(double)); head += sizeof(double);
		memcpy(head,&(array[i].X[0]),sizeof(double)); head += sizeof(double);
		memcpy(head,&(array[i].Y[0]),sizeof(double)); head += sizeof(double);
		memcpy(head,&(array[i].PX[0]),sizeof(double)); head += sizeof(double);
		memcpy(head,&(array[i].PY[0]),sizeof(double)); head += sizeof(double);
		memcpy(head,&(array[i].FX[0]),sizeof(double)); head += sizeof(double);
		memcpy(head,&(array[i].FY[0]),sizeof(double)); head += sizeof(double);
	}
}

void Particle::unpack(std::vector<Particle> &array, int length, int offset, char *data){
	char *head = data;
	for(int i=offset;i<offset+length;i++){
		memcpy(&(array[i].id),head,sizeof(int)); head += sizeof(int);
		memcpy(&(array[i].mass),head,sizeof(double)); head += sizeof(double);
		memcpy(&(array[i].charge),head,sizeof(double)); head += sizeof(double);
		memcpy(&(array[i].X[0]),head,sizeof(double)); head += sizeof(double);
		memcpy(&(array[i].Y[0]),head,sizeof(double)); head += sizeof(double);
		memcpy(&(array[i].PX[0]),head,sizeof(double)); head += sizeof(double);
		memcpy(&(array[i].PY[0]),head,sizeof(double)); head += sizeof(double);
		memcpy(&(array[i].FX[0]),head,sizeof(double)); head += sizeof(double);
		memcpy(&(array[i].FY[0]),head,sizeof(double)); head += sizeof(double);
	}
}