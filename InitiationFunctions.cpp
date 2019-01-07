#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"


///////////////////////////////////
/////////Initiating matirces//////
/////////////////////////////////

// initiate N profiels
void Init_N_1(double mat[N]){
    int i=0;
    for (i=0; i<N; i++) mat[i]= 1.0;
}//Init_N

// initiate N profiels
void Init_N_0(double mat[N]){
    int i=0;
    for (i=0; i<N; i++) mat[i]= 0.0;
}//Init_N

void Init_mutants(double R[N], double L[N]){
    int i = 0;
    
    for (i=0; i<100; i++) {
        R[i] = 0.0;
        L[N-i] = 0.0;
    }
}

void Init_N_rec(double mat[N]){
    int i=0;
    for (i=0; i<N; i++) mat[i]= r_fixed;
}//Init_N


void set_kon(double mat[30]){
    
//    mat[0] = 0.1;
//    mat[1] = 0.25;
//    mat[2] = 0.5;
//    mat[3] = 0.75;
//    mat[4] = 1.0;
//    mat[5] = 1.25;
//    mat[6] = 1.5;
//    mat[7] = 2.0;
//    mat[8] = 2.5;
//    mat[9] = 3.0;
//    mat[10] = 3.5;
//    mat[11] = 4.0;
//    mat[12] = 4.5;
//    mat[13] = 5.0;
//    mat[14] = 6.0;
//    mat[15] = 7.0;
//    mat[16] = 8.0;
    
    mat[0] = 0.01;
    mat[1] = 0.25;
    mat[2] = 0.5;
    mat[3] = 0.75;
    mat[4] = 1.0;
    mat[5] = 1.25;
    mat[6] = 1.5;
    mat[7] = 1.75;
    mat[8] = 2.0;
    mat[9] = 2.25;
    mat[10] = 2.5;
    mat[11] = 2.75;
    mat[12] = 3.0;
    mat[13] = 3.25;
    mat[14] = 3.5;
    mat[15] = 3.75;
    mat[16] = 4.0;
    mat[17] = 4.25;
    mat[18] = 4.5;
    mat[19] = 4.75;
    mat[20] = 5.0;
}

void set_koff(double mat[8]){
    
    mat[0] = 0.1;
    mat[1] = 0.25;
    mat[2] = 0.5;
    mat[3] = 0.75;
    mat[4] = 1.0;
    mat[5] = 2.0;
    mat[6] = 5.0;
    mat[7] = 10.0;
}

void set_r(double mat[7]){
    
    mat[0] = 0.01;
    mat[1] = 0.05;
    mat[2] = 0.1;
    mat[3] = 0.2;
    mat[4] = 0.3;
    mat[5] = 0.4;
    mat[6] = 0.5;
}


void set_f(double mat[15]){
    
    mat[0] = 1.0;
    mat[1] = 10.0;
    mat[2] = 100.0;
    mat[3] = 4.0;
    mat[4] = 5.0;
    mat[5] = 6.0;
    mat[6] = 7.0;
    mat[7] = 8.0;
    mat[8] = 9.0;
    mat[9] = 10.0;
    mat[10] = 15.0;
    mat[11] = 20.0;
    mat[12] = 30.0;
    mat[13] = 40.0;
    mat[14] = 50.0;
    
//    mat[0] = 0.5;
//    mat[1] = 1.0;
//    mat[2] = 2.0;
//    mat[3] = 3.0;
//    mat[4] = 4.0;
//    mat[5] = 5.0;
//    mat[6] = 6.0;
//    mat[7] = 7.0;
//    mat[8] = 8.0;
//    mat[9] = 9.0;
//    mat[10] = 10.0;
}


void Init_numt(double mat[numt]){
    int i=0;
    
    for(i=0; i<numt; i++) mat[i] = 0.0;
}

void Init_N_numt(double mat[N][numt]){
    int i=0;
    int j=0;
    
    for(i=0; i<numt; i++)
        for(j=0; j<N; j++)
            mat[j][i] = 0.0;
}

void Init_cells_SS(double mat[N][6]){
    int i=0;
    int j=0;
    
    for (i=0; i<N; i++)
        for(j=0; j<6; j++)
            mat[i][j]= 0.0;
}//Init_N

void Init_productions(double bL[N], double bR[N], double bl[N], double br[N]){
    int i =0;
    
    for(i=0; i<N; i++){
        if(unifRand()<0.5){
            bL[i] = 1.0 - dv0;
            bl[i] = dv0;
            bR[i] = 1.0;
            br[i] = 0.0;
        }
        else{
            bL[i] = 1.0;
            bl[i] = 0.0;
            bR[i] = 1.0 -dv0;
            br[i] = dv0;
            
        }

        

    }
    
}






