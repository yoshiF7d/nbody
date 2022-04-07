#ifndef CONSTANTS_H
#define CONSTANTS_H

#define DIMENSION 2
#define ORDER 4

//#define EULER
#define PEC
#define REFLECT

const double POTENTIAL_EPSILON=0.1;
const double RADIUS_MIN = 0.5;
const unsigned int SEED = 2284556121;
const double EPSILON=0.1;

const double MASS_PROTON=100;
const double MASS_ELECTRON=0.1;
const double CHARGE_PROTON=10;
const double CHARGE_ELECTRON=-10;
const double MOMENTUM_MAX=0.1;
const double MOMENTUM_MIN=-0.1;

const double ERROR_MAX = 1;
const double ERROR_MIN = 1e-5;
const double STEP_MAX = 1<<8;

#ifdef SYMP2
const int symp_r = 2;
const double symp_k[2] = {
	0.5,
	0.5
};	
const double symp_u[2] = {
	0,
	1
};
#endif
	
#ifdef SYMP4
const int symp_m = 4;
const int symp_r = 4;
	
const double symp_k[4] = {
	+0.67560359597982888591,
	-0.17560359597982888591,
	-0.17560359597982888591,
	+0.67560359597982888591
};
	
const double symp_u[4] = {
	+0.00000000000000000000,
	+1.35120719195965777182,
	-1.70241438391931554364,
	+1.35120719195965777182
};
#endif
	
#ifdef ERK4
#define ORDER 4
const int erk_n = 4;
const double erk_a[3][3] = {
	{0.5,0  ,0},
	{0  ,0.5,0},
	{0  ,0  ,1}
};
const double erk_b[4] = {
	0.166667,0.333333,0.333333,0.166667
};
#endif

#endif