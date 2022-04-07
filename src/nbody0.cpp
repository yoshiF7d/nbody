#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sstream>
#include <iomanip>
#include <cmdline.h>
#include <cmath>

#include <Constants.h>
#include <ParticleList.h>
#include <Particle.h>

int main(int argc, char **argv){
	int n,N;
	double te,dt,radius;
	double tout,out_interval;
	bool cout_on=true;
	std::string dir;
	std::ofstream ofs;
	std::stringstream ss;

	cmdline::parser parser;
	parser.add<int>("number",'n',"particle number",false,10);
	parser.add<double>("interval",'i',"output interval",false,1);
	parser.add<double>("time",'t',"simulation time",false,2);
	parser.add<double>("radius",'r',"wall radius",false,10);
	parser.add<double>("step",'s',"time step",false,1./512);
	parser.footer("outdir");
	parser.parse_check(argc,argv);
	n = N = parser.get<int>("number");
	te = parser.get<double>("time");
	tout = out_interval = parser.get<double>("interval");
	radius = parser.get<double>("radius");
	dt = parser.get<double>("step");
	if(parser.rest().size() > 0){
		dir=parser.rest()[0];
		mkdir(dir.c_str(),0777);
		cout_on=false;
	}

	double range[2][2] = {{-radius,radius},{-radius,radius}};
	//initialize particle parameters
	ParticleList particles(n);
	particles.init(range);

	if(cout_on){
		std::cout << std::setw(8) << "time"; 
		Particle::header();
	}

	//main loop
	int counter=0;
	for(double t=0;t<=te;t+=dt,tout+=dt){
		//initialize force
		
		#if defined(SYMP2) || defined(SYMP4) 
		//std::cout << "rank : " << rank << " energy : " << energy << std::endl;
		//symplectic integration
		//http://cms.phys.s.u-tokyo.ac.jp/~naoki/CIPINTRO/CIP/symplectic.html
		for(int k=0;k<symp_r;k++){
			for(int i=0;i<n;i++){
				particles[i].FX[0] = particles[i].FY[0] = 0;
			}
			particles.evaluate();
			for(int i=0;i<n;i++){
				particles[i].PX[0] += dt*symp_u[k]*particles[i].FX[0];
				particles[i].PY[0] += dt*symp_u[k]*particles[i].FY[0];
			}
			for(int i=0;i<n;i++){
				particles[i].X[0] += dt*symp_k[k]*particles[i].PX[0]/particles[i].mass;
				particles[i].Y[0] += dt*symp_k[k]*particles[i].PY[0]/particles[i].mass;
			}
		}
		#endif
		
		#ifdef EULER
		//explicit euler
		for(int i=0;i<n;i++){
			particles[i].FX[0] = particles[i].FY[0] = 0;
		}
		particles.evaluate();
		for(int i=0;i<n;i++){
			particles[i].PX[0] += dt*particles[i].FX[0];
			particles[i].PY[0] += dt*particles[i].FY[0];
			particles[i].X[0] += dt*particles[i].PX[0]/particles[i].mass;
			particles[i].Y[0] += dt*particles[i].PY[0]/particles[i].mass;
		}
		#endif
		
		#ifdef ERK4
		for(int k=0;k<erk_n-1;k++){
			for(int i=0;i<n;i++){
				particles[i].FX[k] = particles[i].FY[k] = 0;
			}
			particles.evaluate(k);
			for(int i=0;i<n;i++){
				particles[i].X[k+1] = particles[i].X[0];
				particles[i].Y[k+1] = particles[i].Y[0];
				particles[i].PX[k+1] = particles[i].PX[0];
				particles[i].PY[k+1] = particles[i].PY[0];
				
				for(int l=0;l<k+1;l++){
					particles[i].X[k+1] += erk_a[k][l]*dt*particles[i].PX[l]/particles[i].mass;
					particles[i].Y[k+1] += erk_a[k][l]*dt*particles[i].PY[l]/particles[i].mass;
					particles[i].PX[k+1] += erk_a[k][l]*dt*particles[i].FX[l];
					particles[i].PY[k+1] += erk_a[k][l]*dt*particles[i].FY[l];
				}
			}
		}
		for(int i=0;i<n;i++){
			for(int k=0;k<erk_n;k++){
				particles[i].X[0] += erk_b[k]*dt*particles[i].PX[k]/particles[i].mass;
				particles[i].Y[0] += erk_b[k]*dt*particles[i].PY[k]/particles[i].mass;
				particles[i].PX[0] += erk_b[k]*dt*particles[i].FX[k];
				particles[i].PY[0] += erk_b[k]*dt*particles[i].FY[k];
			}
		}
		#endif

		#ifdef PEC
		//predictor corrector
		//http://jun.artcompsci.org/papers/bussei-nbody/node7.html
		for(int i=0;i<n;i++){
			particles[i].FX[0] = particles[i].FY[0] = 0;
			particles[i].YX[0] = particles[i].YY[0] = 0;
		}
		particles.evaluate2();
		for(int i=0;i<n;i++){
			particles[i].move(0,2);
		}
			
		for(int i=0;i<n;i++){
			retry:
			for(int k=0;k<particles[i].step;k++){
				for(int j=0;j<n;j++){
					particles[j].predict(dt/particles[i].step);
				}
				particles[i].FX[1] = particles[i].FY[1] = 0;
				particles[i].YX[1] = particles[i].YY[1] = 0;
				
				for(int j=0;j<n;j++){
					if(i==j){continue;}
					particles[i].evaluate2(particles[j],1);
				}
				particles[i].correct(dt/particles[i].step);
				if(particles[i].error(0,1) > ERROR_MAX){
					if(particles[i].step >= STEP_MAX){break;}
					particles[i].step <<= 1;
					particles[i].move(2,0);
					goto retry;
				}
			}
			if(particles[i].error(0,1) < ERROR_MIN && particles[i].step > 1){
				particles[i].step >>= 1;
			}
			particles[i].move(1,0);
		}
		#endif

		//boundary condition
		particles.boundary(range);
		//output
		if(tout>=out_interval){
			tout=0;
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
			//std::cout <<  "energy : " << std::setprecision(6) << std::scientific << energy << "\n";
			//std::cout << "\033[F\033[J";
			std::cout << "\033[F\033[J";
			std::cout << "\033[F\033[J";
			counter++;
		}
	}
	return 0;
}
