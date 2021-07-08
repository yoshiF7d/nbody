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
#include <mpi.h>
#define SYMP

const double mass_proton=100;
const double mass_electron=0.1;
const double charge_proton=10;
const double charge_electron=-10;
const double momentum_max = 0.001;
const double momentum_min = -0.001;
const double error_max = 1;
const double error_min = 1e-5;
const double step_max = 1<<8;

const int symp_r = 2;
const double symp_k[2] = {
	0.5,
	0.5
};	
	const double symp_u[2] = {
	0,
	1
};

int main(int argc, char **argv){

#ifdef MPI
	int nproc,rank,offset;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
	std::random_device seed_gen;
	std::mt19937 mt(seed_gen());
	int n,N,out_interval;
	double te,dt,wall_radius;
	bool cout_on=true;
	std::string dir;
	std::ofstream ofs;
	std::stringstream ss;
#ifdef MPI
	if(rank == 0){
#endif
	cmdline::parser parser;
	parser.add<int>("number",'n',"particle number",false,10);
	parser.add<int>("interval",'i',"output interval",false,1);
	parser.add<double>("time",'t',"simulation time",false,100);
	parser.add<double>("radius",'r',"wall radius",false,10);
	parser.add<double>("step",'s',"time step",false,1./256);
	parser.parse_check(argc,argv);
    n = N = parser.get<int>("number");
	te = parser.get<double>("time");
	out_interval = parser.get<int>("interval");
	wall_radius = parser.get<double>("radius");
	dt = parser.get<double>("step");
	if(parser.rest().size() > 0){
		dir=parser.rest()[0];
		mkdir(dir.c_str(),0777);
		cout_on=false;
	}
	if(cout_on){std::cout << "time\tid\tstep\tmass\tcharge\tx\ty\tpx\tpy\tfx\tfy\n";}
#ifdef MPI
	}
#endif

#ifdef MPI
	std::vector<Particle> array;
	if(rank==0){
		n = N/nproc;
		MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&te,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&wall_radius,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&out_interval,1,MPI_INT,0,MPI_COMM_WORLD);
		n += N%nproc;
		/*
		for(int i=1;i<nproc;i++){
			offset = n + i*(N/nproc);
			MPI_Send(&offset,1,MPI_INT,i,i,MPI_COMM_WORLD);
		}
		*/
		offset = 0;
		//std::cout << "N : " << N << std::endl;
		array.resize(N);	//rank0 processor's array stores every particles for output purpose
	}else{
		MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&te,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&wall_radius,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&out_interval,1,MPI_INT,0,MPI_COMM_WORLD);
		//MPI_Recv(&offset,1,MPI_INT,0,rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		array.resize(n);
	}
	offset = 0;
#else
	std::vector<Particle> array(n);
#endif
	std::uniform_real_distribution<> radius_random(0,wall_radius);
	std::uniform_real_distribution<> theta_random(0,2*M_PI);
	std::uniform_real_distribution<> momentum_set(momentum_min,momentum_max);
	double xmax = wall_radius;
	double xmin = -wall_radius;
	double ymax = wall_radius;
	double ymin = -wall_radius;
	//initialize particle parameters
#ifdef MPI
	if(rank%2==0){
		for(int i=0;i<n;i++){
			array[i].set(offset + i,mass_proton,charge_proton,
				{radius_random(mt)*cos(theta_random(mt)),radius_random(mt)*sin(theta_random(mt)),0,0,0,0},
				{0,0,0,0,0,0}
			);
		}
		if(rank==0 && nproc%2!=0){
			for(int i=n/2;i<n;i++){
				array[i].set(offset + i,mass_electron,charge_electron,
					{radius_random(mt)*cos(theta_random(mt)),radius_random(mt)*sin(theta_random(mt)),0,0,0,0},
					{0,0,0,0,0,0}
				);
			}
		}
	}else{
		for(int i=0;i<n;i++){
			array[i].set(offset + i,mass_electron,charge_electron,
				{radius_random(mt)*cos(theta_random(mt)),radius_random(mt)*sin(theta_random(mt)),0,0,0,0},
				{0,0,0,0,0,0}
			);
		}
	}
#else
	for(int i=0;i<n/2;i++){
		array[i].set(
			i,mass_proton,charge_proton,
			{radius_random(mt)*cos(theta_random(mt)),radius_random(mt)*sin(theta_random(mt)),0,0,0,0},
			{0,0,0,0,0,0}
		);
	}
	for(int i=n/2;i<n;i++){
		array[i].set(
			i,mass_electron,charge_electron,
			{radius_random(mt)*cos(theta_random(mt)),radius_random(mt)*sin(theta_random(mt)),0,0,0,0},
			{0,0,0,0,0,0}
		);
	}
#endif
	double energy;
	int counter=0;
	std::vector<std::array<double, 2> > enarray;
	//main loop
	for(double t=0;t<=te;t+=dt){
		//calculate inter particle force
		energy = 0;
		for(int i=0;i<n;i++){
			array[i].force[0] = array[i].force[1] = 0;
			array[i].yank[0] = array[i].yank[1] = 0;
			energy += (array[i].momentum[0])*(array[i].momentum[0])/(2*array[i].mass);
			energy += (array[i].momentum[1])*(array[i].momentum[1])/(2*array[i].mass);
		}
		for(int i=0;i<n;i++){
			//coulomb force
			for(int j=0;j<n;j++){
				if(i==j){continue;}
				energy += array[i].interact(array[j]);
			}
		}
#ifdef MPI
		char *data;
		int isize;
		size_t size;
		std::vector<Particle> array2;
		for(int i=0;i<nproc;i++){
			if(i!=rank){
				MPI_Status status;
				MPI_Probe(i,i,MPI_COMM_WORLD,&status);
				MPI_Get_count(&status,MPI_CHAR,&isize);
				size = isize;
				data = (char*)malloc(size);
				MPI_Recv(data,size,MPI_CHAR,i,i,MPI_COMM_WORLD,&status);
				array2 = Particle::deserialize(data,size);
				for(int i=0;i<n;i++){
					//coulomb force from particles of other process.
					for(int j=0;j<array2.size();j++){
						energy += array[i].interact(array2[j]);
					}
				}
			}else{
				Particle::serialize(array,n,&data,&size);
				for(int j=0;j<nproc;j++){
					if(j!=i){
						MPI_Ssend(data,size,MPI_CHAR,j,i,MPI_COMM_WORLD);
					}
				}
			}
		}
#endif

#ifdef PC
		//predictor corrector
		//http://jun.artcompsci.org/papers/bussei-nbody/node7.html#SECTION00070000000000000000
		for(int i=0;i<n;i++){
			std::array<double,6> position = array[i].position;
			std::array<double,6> momentum = array[i].momentum;
			retry:
			for(int k=0;k<array[i].step;k++){
				for(int j=0;j<n;j++){
					array[j].predict(dt/array[i].step);
				}
				array[i].force[3] = array[i].force[4] = 0;
				array[i].yank[3] = array[i].yank[4] = 0;
				for(int j=0;j<n;j++){
					if(i==j){continue;}
					array[i].evaluate(array[j]);
				}
				array[i].correct(dt/array[i].step);
				array[i].push();
				//std::cout << "error : " << array[i].error() << "\n";
				if(array[i].error() > error_max){
					//std::cout << array[i].step << "\n";
					if(array[i].step >= step_max){break;}
					array[i].step <<= 1;
					array[i].position = position;
					array[i].momentum = momentum;
					goto retry;
				}
			}
			if(array[i].error() < error_min && array[i].step > 1){
				//std::cout << array[i].step << "\n";
				array[i].step >>= 1;
			}
		}
#ifdef MPI
		for(int i=0;i<nproc;i++){
			if(i!=rank){
				MPI_Status status;
				MPI_Probe(i,i,MPI_COMM_WORLD,&status);
				MPI_Get_count(&status,MPI_CHAR,&isize);
				size = isize;
				data = (char*)malloc(size);
				MPI_Recv(data,size,MPI_CHAR,i,i,MPI_COMM_WORLD,&status);
				array2 = Particle::deserialize(data,size);
				for(int i=0;i<n;i++){
					for(int j=0;j<array2.size();j++){
						array[i].evaluate(array2[j]);
					}
				}
			}else{
				Particle::serialize(array,n,&data,&size);
				for(int j=0;j<nproc;j++){
					if(j!=i){
						MPI_Ssend(data,size,MPI_CHAR,j,i,MPI_COMM_WORLD);
					}
				}
			}
		}
#endif
		for(int i=0;i<n;i++){
			array[i].correct(dt);
			array[i].push();
		}
#endif
#ifdef SYMP
		//std::cout << "rank : " << rank << " energy : " << energy << std::endl;
		//symplectic integration
		//http://cms.phys.s.u-tokyo.ac.jp/~naoki/CIPINTRO/CIP/symplectic.html
		for(int k=0;k<symp_r;k++){
			for(int i=0;i<n;i++){
				array[i].momentum[0] += dt*symp_u[k]*array[i].force[0];
				array[i].momentum[1] += dt*symp_u[k]*array[i].force[1];
			}
			for(int i=0;i<n;i++){
				array[i].position[0] += dt*symp_k[k]*array[i].momentum[0]/array[i].mass;
				array[i].position[1] += dt*symp_k[k]*array[i].momentum[1]/array[i].mass;
			}
		}
#endif
		//explicit euler
		/*
		for(int i=0;i<N;i++){
			array[i]->px += dt*array[i]->fx;
			array[i]->py += dt*array[i]->fy;
			array[i]->x += dt*array[i]->px/array[i]->m;
			array[i]->y += dt*array[i]->py/array[i]->m;
		}
		*/
		//boundary condition
		for(int i=0;i<n;i++){
			if(array[i].position[0] > xmax){array[i].position[0] -= xmax - xmin;}
			if(array[i].position[0] < xmin){array[i].position[0] += xmax - xmin;}
			if(array[i].position[1] > ymax){array[i].position[1] -= ymax - ymin;}
			if(array[i].position[1] < ymin){array[i].position[1] += ymax - ymin;}
		}
		counter++;
		if(counter%out_interval==0){
#ifdef MPI
		if(rank==0){
#endif
			enarray.push_back(std::array<double,2>{t,energy});
#ifdef MPI
		}
#endif

#ifdef MPI
		if(rank==0){
			//std::cout << "n : " << n << std::endl;
			int nn = n;
			int isize;
			for(int i=1;i<nproc;i++){
				MPI_Status status;
				MPI_Probe(i,i,MPI_COMM_WORLD,&status);
				MPI_Get_count(&status,MPI_CHAR,&isize);
				size = isize;
				data = (char*)malloc(size);
				MPI_Recv(data,size,MPI_CHAR,i,i,MPI_COMM_WORLD,&status);
				Particle::deserialize_fill(array,&nn,data,size);
				double energy2;
				MPI_Recv(&energy2,1,MPI_DOUBLE,i,i+nproc,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				energy += energy2;
			}
#endif
			if(cout_on){
				for(int i=0;i<N;i++){
					std::cout << t << "\t"; 
					array[i].print();
				}
				std::cout << "\n";
			}else{
				ss.str("");
				ss << std::setfill('0') << std::setw(5) << counter/out_interval;
				ofs.open(dir + std::string("/") + ss.str() + std::string(".data"),std::ios::binary);
				for(int i=0;i<N;i++){
					ofs.write((char*)&t,sizeof(double));
					ofs.write((char*)&(array[i].id),sizeof(int));
					ofs.write((char*)&(array[i].mass),sizeof(double));
					ofs.write((char*)&(array[i].charge),sizeof(double));
					ofs.write((char*)&(array[i].position[0]),sizeof(double));
					ofs.write((char*)&(array[i].position[1]),sizeof(double));
					ofs.write((char*)&(array[i].momentum[0]),sizeof(double));
					ofs.write((char*)&(array[i].momentum[1]),sizeof(double));
					ofs.write((char*)&(array[i].force[0]),sizeof(double));
					ofs.write((char*)&(array[i].force[1]),sizeof(double));
				}
				ofs.close();
			}
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
#ifdef MPI
		}else{
			Particle::serialize(array,n,&data,&size);
			MPI_Ssend(data,size,MPI_CHAR,0,rank,MPI_COMM_WORLD);
			MPI_Ssend(&energy,1,MPI_DOUBLE,0,rank+nproc,MPI_COMM_WORLD);
		}
#endif
		}
	}
#ifdef MPI
	if(rank==0){
#endif
	ofs.open(std::string("energy.txt"));
	for(int i=0;i<enarray.size();i++){
		ofs << enarray[i][0] << "\t" << enarray[i][1] << "\n";
	}
	ofs.close();
#ifdef MPI
	}
#endif

#ifdef MPI
	MPI_Finalize();
#endif
	return 0;
}
