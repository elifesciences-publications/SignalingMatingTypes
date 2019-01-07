#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <sys/stat.h>

#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"
#include "WriteFiles.hpp"
#include "Signal.hpp"


using namespace std;

// declare parameters
// system parameters
double gR = 0.5;         // degradation rate for receptors
double gL = 0.5;         // degradaton rate for ligand
double gLR = 0.5;         // degradation rate for receptor/ligand pair

double maxAlpha = 0.75;


double kon = 1.0;
double koff = 1.0;
double f=1.0;
double kb = f*kon/koff;
double K=20.0;
double n=1.0;

double epsilon = 0.0000000000001;//0.00000001;
double epsilonR = 0.001;
double mu = 0.0;

double mu_r = 0.0;
double sigma_r = 0.1;

// mutation rates
double nu_bR = 0.005;
double nu_bL = 0.005;
double nu_br = 0.00;
double nu_bl = 0.00;

double sigma = 0.1;

double nu_recomb = 0.00;
double r_fixed = 0.0;


double s_before = 10.0;
double s_now = 1.0;
double dr = 10.0;
double r_before = 100.0; // the mean recombination rate before

int stabCheck = 0;
int countStab = 0;

double normal_base = 2.0;

// initial conditions
double dv0 = 0.1;

int ind = 0;

int k=0;
int l=0;
int m=0;

//Add objects
Signal mySignal;
//////////////////////////////////////////////////////////////////////////////
////////////////////////// M A I N   P R O G R A M //////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv) {
	
	srand((unsigned)time(0));
    ////// MAKE DIRECTORY WITH PARAMETER NAMES ////////////////////////////////////////////////////////////////////////
    // char name[550];
    //comment below for cluster
    //sprintf(name, "/Users/zenah/Dropbox/Signalling_sexes/C++_code/SS_agent_based/sims_N%.d_M%.d_q%.3f_mu%.3f_sigma%.3f_nuL1_%.3f_nuR1_%.3f_nuL2_%.3f_nuR2_%.3f/", N, M, q, mu, sigma, nu_bL1, nu_bR1, nu_bL2, nu_bR2);
    //mkdir(name, S_IRWXU);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int repeats = 0;
    
    // read parameters in
    //kon = atof(argv[1]);
    //dv0 = atof(argv[2]);
    //kb = f*kon/koff;
    //K = atof(argv[3]);
    //n = atof(argv[3]);
    
    // fixed parameters (no cluster)

    // params
    double konList[30];
    double fList[15];
    double rList[7];
    
    set_kon(konList);
    set_f(fList);
    set_r(rList);

    
    int i = 0;
    int j = 0;
    int l = 0;
    int k = 0;
    // dv0 = 0.01;
    for(repeats=17; repeats<18; repeats++){
        for(k=0; k<1; k++){
            //n = 1.0*k;
            r_fixed = 0.0;
            //r_fixed = rList[k];
            dv0 = 0.0;//k*0.1;
            for (l=2; l<3; l++) {
                nu_bR = 0.01;//pow(10.0, -1.0*l);
                nu_bL = nu_bR;
                nu_recomb = 0.0;
                for(i=0; i<1; i++){
                    f = 1.0;//fList[i];
                    for(j=19; j<20; j++){
                        kon = 5.0;//konList[j];
                        kb = f*kon/koff;
                        // main loop
                        K = 0.5*mySignal.GetMaxK();
                        printf("K = %.13f\n", K);
                        // initiate recombination locus to 0
                        Init_N_rec(mySignal.recomb);
                        // initiate bR1, bL1, bR2, bL2 kc, kt, gR[N], gL[N] in all N cells
                        Init_productions(mySignal.bL, mySignal.bR, mySignal.bl, mySignal.br);
                        //Init_N_1(mySignal.bL);
                        //Init_N_1(mySignal.bR);
                        //Init_N_0(mySignal.bl);
                        //Init_N_0(mySignal.br);
                        // initiate evolution matrices
                        Init_numt(mySignal.s_Evolution);
                        Init_numt(mySignal.r_Evolution);
                        
                        //int
                        ind = 0;
                        stabCheck = 0;
                        countStab = 0;
                        mySignal.s_Evolution[ind] = mySignal.Get_S(mySignal.bL, mySignal.bR, mySignal.bl, mySignal.br);
                        mySignal.r_Evolution[ind] = mySignal.get_mean_r(mySignal.recomb);

                        while((stabCheck!=1 || ind < 10000) && ind< 100000){
                            if(ind>1000){
                                nu_bR = 0.01;//pow(10.0, -1.0*l);
                                nu_bL = nu_bR;
                            }
                            //mySignal.Update_si_ri(mySignal.si_ri_Evo, mySignal.bL1, mySignal.bR1, mySignal.recomb, ind);
                            /////
                            mySignal.MUTATE(mySignal.bL, mySignal.bR, mySignal.bl, mySignal.br, mySignal.recomb);
                            // mate
                            mySignal.MATE(mySignal.bL, mySignal.bR, mySignal.bl, mySignal.br, mySignal.recomb);
                            // check stability by evaluating total |bR - bL|
                            if(countStab > numSave){
                                countStab = 0;
                                s_now = mySignal.Get_S(mySignal.bL, mySignal.bR, mySignal.bl, mySignal.br);
                                dr = Absolute(r_before, mySignal.get_mean_r(mySignal.recomb));
                                if(Absolute(s_now, s_before) < epsilon && ind>10 && dr<epsilonR) stabCheck = 1;
                                //printf("ds:%.20f      dr:%.10f  ing:%.d\n", Absolute(s_now,s_before), dr, ind);
                                s_before = s_now;
                                r_before = mySignal.get_mean_r(mySignal.recomb);
                                // save profiles dynamically
                            //    Write_Profile(mySignal.bL1, ind, 1);
                            //    Write_Profile(mySignal.bR1, ind, 2);
                            }//if
                            ind++;
                            //printf("here: %.d\n", ind);
                            countStab++;
                            mySignal.s_Evolution[ind] = mySignal.Get_S(mySignal.bL, mySignal.bR, mySignal.bl, mySignal.br);
                            mySignal.r_Evolution[ind] = mySignal.get_mean_r(mySignal.recomb);
                            
                        }//while
                    

                    Write_Profile(mySignal.bL, ind, 1, repeats);
                    Write_Profile(mySignal.bR, ind, 2, repeats);
                    Write_Profile(mySignal.bl, ind, 3, repeats);
                    Write_Profile(mySignal.br, ind, 4, repeats);
//                    Write_recomb(mySignal.recomb, ind, repeats);
            //           Write_Evolution(mySignal.s_Evolution, ind, 1, repeats);
                 //     Write_Evolution(mySignal.r_Evolution, ind, 2, repeats);
                    }//i
                }//j
            }//l
        }//k
    }
	return 0;
}



