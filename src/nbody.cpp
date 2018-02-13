#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <list>
#include <Particle.h>
#include <sys/stat.h>
#include <sstream>
#include <iomanip>
#include <cmdline.h>

const double dt = 10.0/256;
const double te = 1000;
const double wall_radius = 10;
const double epsilon=0.001;
const int out_interval = 100;

const int symp_r = 2;
const double symp_k[2]={
	0.5,
	0.5
};
const double symp_u[2]={
	0,
	1
};

int main(int argc, char **argv){
	std::random_device seed_gen;
	std::mt19937 mt(seed_gen());
	std::uniform_real_distribution<> mass_set(0.1,1.0);
	std::uniform_real_distribution<> charge_set(-0.1,0.1);
	std::uniform_real_distribution<> position_set(-1,1);
	std::uniform_real_distribution<> momentum_set(-0.1,0.1);
	
	std::ofstream ofs;
	std::stringstream ss;
	cmdline::parser parser;
	
	parser.add<int>("number",'n',"particle number",false,10);
	parser.parse_check(argc,argv);
	
	int N = 10;
	N=parser.get<int>("number");
	
	bool cout_on=true;

	if(parser.rest().size() > 0){
		mkdir(parser.rest()[0].c_str(),0777);
		cout_on=false;
	}
	
	if(cout_on){std::cout << "time\tid\tmass\tcharge\tx\ty\tpx\tpy\tfx\tfy\n";}
	
	std::vector<Particle*> array(N);
	//initialize particle parameters	
	for(int i=0;i<N;i++){
		array[i] = new Particle(
			i,
			mass_set(mt),
			charge_set(mt),
			position_set(mt),
			position_set(mt),
			momentum_set(mt),
			momentum_set(mt)
		);
	}
	
	double r,x,y,th,energy;
	int counter=0;
	//main loop
	for(double t=0;t<=te;t+=dt){
		//calculate inter perticle force
		energy = 0;
		for(int i=0;i<N;i++){
			array[i]->px = 0;
			array[i]->py = 0;
			energy += (array[i]->px)*(array[i]->px)/(2*array[i]->m);
			energy += (array[i]->py)*(array[i]->py)/(2*array[i]->m);
			//coulomb force
			for(int j=i+1;j<N;j++){
				x = array[i]->x - array[j]->x;
				y = array[i]->y - array[j]->y;
				r = sqrt(x*x+y*y+epsilon);
				array[i]->fx += x*((array[i]->q)*(array[j]->q))/(r*r*r);
				array[i]->fy += y*((array[i]->q)*(array[j]->q))/(r*r*r);
				array[j]->fx -= x*((array[i]->q)*(array[j]->q))/(r*r*r);
				array[j]->fy -= y*((array[i]->q)*(array[j]->q))/(r*r*r);
				energy += (array[i]->q)*(array[j]->q)/abs(r);
			}
			//force from potential well
			x = array[i]->x;
			y = array[i]->y;
			r = wall_radius;
			th = tanh((x*x+y*y-r*r));
			array[i]->fx -= 10*x*(1-th*th);
			array[i]->fy -= 10*y*(1-th*th);
			energy += 10*th;
		}
		//symplectic integration
		//http://www-cms.phys.s.u-tokyo.ac.jp/~naoki/CIPINTRO/CIP/symplectic.html
		for(int r=0;r<symp_r;r++){
			for(int i=0;i<N;i++){
				array[i]->px += dt*symp_u[r]*array[i]->fx;
				array[i]->py += dt*symp_u[r]*array[i]->fy;
			}
			for(int i=0;i<N;i++){
				array[i]->x += dt*symp_k[r]*array[i]->px/array[i]->m;
				array[i]->y += dt*symp_k[r]*array[i]->py/array[i]->m;
			}
		}
		//explicit euler
		/*
		for(int i=0;i<N;i++){
			array[i]->px += dt*array[i]->fx;
			array[i]->py += dt*array[i]->fy;
			array[i]->x += dt*array[i]->px/array[i]->m;
			array[i]->y += dt*array[i]->py/array[i]->m;
		}
		*/
		//display current progress to stdout
		counter++;
		if(counter%10==0){
			std::cout << "[";
			for(int i=0;i<100;i++){
				if(i<100*t/te){std::cout << "=";}
				else if(i == ceil(100*t/te)){std::cout << "🐌";}
				else{std::cout << " ";}
			}
				
			std::cout << "]\n"; 
			std::cout << "t : " << t << "/" << te << "\n";
			std::cout <<  "energy : " << std::setprecision(6) << std::scientific << energy << "\n";
			std::cout << "\033[F\033[J";
			std::cout << "\033[F\033[J";
			std::cout << "\033[F\033[J";
		}
		//output
		if(counter%out_interval==0){
			if(cout_on){
				for(int i=0;i<N;i++){
					std::cout << t << "\t"; 
					array[i]->print();
				}
				std::cout << "\n";
			}else{
				ss.str("");
				ss << std::setfill('0') << std::setw(5) << counter;
				ofs.open(std::string(argv[1]) + std::string("/") + ss.str() + std::string(".dat"),std::ios::binary);
				for(int i=0;i<N;i++){
					ofs.write((char*)&t,sizeof(double));
					ofs.write((char*)&(array[i]->id),sizeof(int));
					ofs.write((char*)&(array[i]->m),sizeof(double));
					ofs.write((char*)&(array[i]->q),sizeof(double));
					ofs.write((char*)&(array[i]->x),sizeof(double));
					ofs.write((char*)&(array[i]->y),sizeof(double));
					ofs.write((char*)&(array[i]->px),sizeof(double));
					ofs.write((char*)&(array[i]->py),sizeof(double));
					ofs.write((char*)&(array[i]->fx),sizeof(double));
					ofs.write((char*)&(array[i]->fy),sizeof(double));
				}
				ofs.close();	
			}
		}
	}
	return 0;
}
