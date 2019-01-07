#ifndef GLOBAL
#define GLOBAL //Global.hpp


//Define global fixed parameters
#define N (1000)                  //number of cells
#define M (500)                   //number of cells that mate in first round
#define numSave (100)
#define numt (20000000)            // maximum generations for the simulation
#define Pi (3.141592653589793238462643383279502884197169399375)
#define num_m (100);

//define global parameters that can be defined and modified in main
extern double gR;         // degradation rate for receptors
extern double gL;         // degradaton rate for ligand
extern double gLR;         // degradation rate for receptor/ligand pair

extern double epsilon;

extern double kon;
extern double koff;
extern double kb;
extern double f;
extern double K;
extern double n;

extern double maxAlpha;

extern double mu;
extern double sigma;

// mutation size for recombinant locues
extern double mu_r;
extern double sigma_r;

// mutation rates
extern double nu_bR;
extern double nu_bL;
extern double nu_br;
extern double nu_bl;
extern double nu_recomb;
extern double r_fixed; //fixed recombination rate

extern double normal_base;

extern double dv0;


#endif
