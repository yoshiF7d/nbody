#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sstream>
#include <iomanip>
#include <cmdline.h>
#include <mpi.h>
#include <cmath>

#include <Constants.h>
#include <ParticleList.h>
#include <Particle.h>

int main(int argc, char **argv){
	int nproc=1,rank=0;
#ifdef MPI
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
	int size,n0,n,N;
	double te,dt,radius;
	double tout,out_interval;
	bool cout_on=true;
	std::string dir;
	std::ofstream ofs;
	std::stringstream ss;
	ParticleList particles;
	ParticleList particles2;

	if(rank == 0){
		cmdline::parser parser;
		parser.add<int>("number",'n',"particle number",false,10);
		parser.add<int>("interval",'i',"output interval",false,1);
		parser.add<double>("time",'t',"simulation time",false,2);
		parser.add<double>("radius",'r',"wall radius",false,10);
		parser.add<double>("step",'s',"time step",false,1./512);
		parser.footer("outdir");
		parser.parse_check(argc,argv);
		size = n0 = n = N = parser.get<int>("number");
		te = parser.get<double>("time");
		tout = out_interval = parser.get<int>("interval");
		radius = parser.get<double>("radius");
		dt = parser.get<double>("step");
		if(parser.rest().size() > 0){
			dir=parser.rest()[0];
			mkdir(dir.c_str(),0777);
			cout_on=false;
		}
		if(cout_on){
			std::cout << std::setw(8) << "time";
			Particle::header();
		}
	}
#ifdef MPI
	/*
		array size:
			rank0 : n0
			others : n
		total size:
			N : n0 + nproc*n (n0 particles at front)
	*/
	n = N/nproc;
	n0 = n + N%nproc;
	MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&n0,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&te,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&radius,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&out_interval,1,MPI_INT,0,MPI_COMM_WORLD);
#endif

	double range[2][2] = {{-radius,radius},{-radius,radius}};
	if(rank==0){
		//rank0 particles stores every particles info
		particles.list.resize(N);	
		particles.init(range);
		size = n0;
	}else{
		particles.list.resize(n);
		size = n;
	}
	//particles2 is for the communication.
	particles2.list.resize(n0);
	particles2.setBuffer();

#ifdef MPI
	if(rank==0){
		/*distribute*/
		int offset = n0;
		for(int i=1;i<nproc;i++,offset+=n){
			particles.pack(particles2.buffer,n,offset);
			MPI_Ssend(particles2.buffer,particles2.bufferSize,MPI_CHAR,i,0,MPI_COMM_WORLD);
		}
	}else{
		MPI_Recv(particles2.buffer,particles2.bufferSize,MPI_CHAR,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		particles.unpack(particles2.buffer,n);
	}
#endif

	int counter=0;
	//main loop
	for(double t=0;t<=te;t+=dt,tout+=dt){
		//initialize force
		for(int i=0;i<size;i++){
			particles[i].FX[0] = particles[i].FY[0] = 0;
		}
		//calculate force
		#ifdef EULER
		//explicit euler
		#ifdef MPI
		for(int i=0;i<nproc;i++){
			if(i!=rank){
				MPI_Recv(particles2.buffer,particles2.bufferSize,MPI_CHAR,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				particles2.unpack();
				if(i==0){
					particles.evaluate(particles2,n0,0);
				}else{
					particles.evaluate(particles2,n,0);
				}
			}else{
				particles.pack(particles2.buffer,size);
				for(int j=0;j<nproc;j++){
					if(j!=i){
						MPI_Ssend(particles2.buffer,particles2.bufferSize,MPI_CHAR,j,i,MPI_COMM_WORLD);
					}
				}
			}
		}
		#endif
		
		particles.evaluate();
		
		for(int i=0;i<n0;i++){
			particles[i].PX[0] += dt*particles[i].FX[0];
			particles[i].PY[0] += dt*particles[i].FY[0];
			particles[i].X[0] += dt*particles[i].PX[0]/particles[i].mass;
			particles[i].Y[0] += dt*particles[i].PY[0]/particles[i].mass;
		}
		#endif
		
		//boundary condition
		particles.boundary(range);
		//output
		if(tout>=out_interval){
			tout=0;
				
			if(rank==0){
				int offset = n0;
				#ifdef MPI
				for(int i=1;i<nproc;i++,offset+=n){
					MPI_Recv(particles2.buffer,particles.bufferSize,MPI_CHAR,i,i,MPI_COMM_WORLD,&status);
					particles.unpack(particles2.buffer,n,offset);
				}
				#endif
				if(cout_on){
					for(int i=0;i<N;i++){
						std::cout << std::setw(5) << std::scientific << std::setprecision(2) <<  t; 
						particles[i].print();
					}
					std::cout << "\n";
				}else{
					ss.str("");
					ss << std::setfill('0') << std::setw(4) << counter;
					ofs.open(dir + std::string("/") + ss.str() + std::string(".data"),std::ios::binary);
					for(int i=0;i<N;i++){
						ofs.write((char*)&t,sizeof(double));
						particles[i].write(ofs);
					}
					ofs.close();
				}
				std::cout << "[";
				for(int i=0;i<100;i++){
					if(i<100*t/te){std::cout << "=";}
					else if(i == std::ceil(100*t/te)){std::cout << "🐌";}
					else{std::cout << " ";}
				}
				std::cout << "]\n"; 
				std::cout << "t : " << t << "/" << te << "\n";
				std::cout << "\033[F\033[J";
				std::cout << "\033[F\033[J";
			}else{
				#ifdef MPI
				particles.pack(particles2.buffer,n);
				MPI_Ssend(particles2.buffer,size,MPI_CHAR,0,rank,MPI_COMM_WORLD);
				#endif
			}
			counter++;
		}
	}
#ifdef MPI
	MPI_Finalize();
#endif
	return 0;
}
