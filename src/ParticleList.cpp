#include <ParticleList.h>
#include <Constants.h>
#include <iomanip>
#include <random>

ParticleList::ParticleList(int size){
	this->list.resize(size);
}

ParticleList::~ParticleList(){
	free(this->buffer);
	this->buffer = nullptr;
}

Particle &ParticleList::operator[] (int n){
	return this->list[n];
}

void ParticleList::evaluate(){
	this->evaluate(*this,this->list.size(),0);
}
void ParticleList::evaluate(int index){
	this->evaluate(*this,this->list.size(),index);
}
void ParticleList::evaluate(int size, int index){
	this->evaluate(*this,size,index);
}
void ParticleList::evaluate(ParticleList &particles){
	this->evaluate(particles,particles.list.size(),0);
}
void ParticleList::evaluate(ParticleList &particles,int index){
	this->evaluate(particles,particles.list.size(),index);
}
void ParticleList::evaluate(ParticleList &particles,int size,int index){
	for(int i=0;i<this->list.size();i++){
		for(int j=0;j<size;j++){
			if(i==j){continue;}
			this->list[i].evaluate(particles.list[j],index);
		}
	}
}

void ParticleList::evaluate2(){
	this->evaluate2(*this,this->list.size(),0);
}
void ParticleList::evaluate2(ParticleList &particles,int size,int index){
	for(int i=0;i<this->list.size();i++){
		for(int j=0;j<size;j++){
			if(i==j){continue;}
			this->list[i].evaluate2(particles.list[j],index);
		}
	}
}


void ParticleList::init(double range[][2]){
	std::mt19937 mt(SEED);
	std::uniform_real_distribution<> radius_random(0,range[0][1]);
	std::uniform_real_distribution<> theta_random(0,2*M_PI);
	std::uniform_real_distribution<> momentum_set(MOMENTUM_MIN,MOMENTUM_MAX);
	int n = this->list.size();
	for(int i=0;i<n/2;i++){
		this->list[i].id = i;
		this->list[i].step = 1;
		this->list[i].mass = MASS_PROTON;
		this->list[i].charge = CHARGE_PROTON;
		double r = radius_random(mt);
		double th = theta_random(mt);
		this->list[i].X[0] = r*cos(th);
		this->list[i].Y[0] = r*sin(th);
		this->list[i].PX[0] = momentum_set(mt);
		this->list[i].PY[0] = momentum_set(mt);
	}
	for(int i=n/2;i<n;i++){
		this->list[i].id = i;
		this->list[i].step = 1;
		this->list[i].mass = MASS_ELECTRON;
		this->list[i].charge = CHARGE_ELECTRON;
		double r = radius_random(mt);
		double th = theta_random(mt);
		this->list[i].X[0] = r*cos(th);
		this->list[i].Y[0] = r*sin(th);
		this->list[i].PX[0] = momentum_set(mt);
		this->list[i].PY[0] = momentum_set(mt);
	};
}
void ParticleList::boundary(double range[][2]){
	#ifdef REFLECT
	for(int i;i<this->list.size();i++) {
		this->list[i].reflect(range);
	}
	#endif
	#ifdef PERIODIC
	for(int i;i<this->list.size();i++) {
		this->list[i].periodic(range);
	}
	#endif
}

void ParticleList::setBuffer(){
	this->setBuffer(this->list.size());
}
void ParticleList::setBuffer(int size){
	this->bufferSize = (sizeof(int) + sizeof(double)*8)*size;
	this->buffer = (char*)malloc(this->bufferSize);
}

void ParticleList::pack(){
	this->pack(this->buffer,this->list.size(),0);
}
void ParticleList::pack(char *data){
	this->pack(this->buffer,this->list.size(),0);
}
void ParticleList::pack(int size){
	this->pack(this->buffer,size,0);
}
void ParticleList::pack(char *data, int size){
	this->pack(this->buffer,size,0);
}
void ParticleList::pack(char *data, int size, int offset){
	char *head = data;
	for(int i=offset;i<offset+size;i++){
		memcpy(head,&(this->list[i].id),sizeof(int)); head += sizeof(int);
		memcpy(head,&(this->list[i].mass),sizeof(double)); head += sizeof(double);
		memcpy(head,&(this->list[i].charge),sizeof(double)); head += sizeof(double);
		memcpy(head,&(this->list[i].X[0]),sizeof(double)); head += sizeof(double);
		memcpy(head,&(this->list[i].Y[0]),sizeof(double)); head += sizeof(double);
		memcpy(head,&(this->list[i].PX[0]),sizeof(double)); head += sizeof(double);
		memcpy(head,&(this->list[i].PY[0]),sizeof(double)); head += sizeof(double);
		memcpy(head,&(this->list[i].FX[0]),sizeof(double)); head += sizeof(double);
		memcpy(head,&(this->list[i].FY[0]),sizeof(double)); head += sizeof(double);
	}
}

void ParticleList::unpack(){
	this->unpack(this->buffer,this->list.size(),0);
}
void ParticleList::unpack(char *data){
	this->unpack(data,this->list.size(),0);
}
void ParticleList::unpack(int size){
	this->unpack(this->buffer,size,0);
}
void ParticleList::unpack(char *data, int size){
	this->unpack(data,size,0);
}
void ParticleList::unpack(char *data, int size, int offset){
	char *head = data;
	for(int i=offset;i<offset+size;i++){
		memcpy(&(this->list[i].id),head,sizeof(int)); head += sizeof(int);
		memcpy(&(this->list[i].mass),head,sizeof(double)); head += sizeof(double);
		memcpy(&(this->list[i].charge),head,sizeof(double)); head += sizeof(double);
		memcpy(&(this->list[i].X[0]),head,sizeof(double)); head += sizeof(double);
		memcpy(&(this->list[i].Y[0]),head,sizeof(double)); head += sizeof(double);
		memcpy(&(this->list[i].PX[0]),head,sizeof(double)); head += sizeof(double);
		memcpy(&(this->list[i].PY[0]),head,sizeof(double)); head += sizeof(double);
		memcpy(&(this->list[i].FX[0]),head,sizeof(double)); head += sizeof(double);
		memcpy(&(this->list[i].FY[0]),head,sizeof(double)); head += sizeof(double);
	}
}
