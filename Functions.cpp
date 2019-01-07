#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"

//////////////
//FUNTIONS //
/////////////

//rand_Normal
double rand_Normal (double mean, double sigma){
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;
    
    if (call == 1)
    {
        call = !call;
        return (mean + sigma * (double) X2);
    }
    
    do
    {
        U1 = -1 + ((double) rand () / RAND_MAX) * 2;
        U2 = -1 + ((double) rand () / RAND_MAX) * 2;
        W = pow (U1, 2) + pow (U2, 2);
    }
    while (W >= 1 || W == 0);
    
    mult = sqrt ((-2 * log (W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;
    
    call = !call;
    
    return (mean + sigma * (double) X1);
}//rand_Normal

//create random number from 0 to 1
double unifRand()
{
    return rand() / double(RAND_MAX);
}//unifRand

//random integer from rangeLow to rangeHight
int uniform_distribution(int rangeLow, int rangeHigh)
{
    int myRand = (int)rand();
    int range = rangeHigh - rangeLow + 1; //+1 makes it [rangeLow, rangeHigh], inclusive.
    int myRand_scaled = (myRand % range) + rangeLow;
    return myRand_scaled;
}//uniform_distribution

//absolute value of difference between two doubles
double Absolute(double a, double b){
    if(a>b) return(a-b);
    else return(b-a);
}//Absolute

double min(double a, double b){
    if(a>b) return(b);
    else return(a);
}//min

double max(double a, double b){
    if(a<b) return(b);
    else return(a);
}//max


double delta_ij_kl(int i, int j, int k, int l){
    if(i==k && j==l) return(1.0);
    else return(0.0);
}//

void copyNNtoNN(double matOriginal[N][N], double matCopy[N][N]){
    int i = 0;
    int j = 0;
    
    for(i=0; i<N; i++) {
        for(j=0; j<N; j++){
            matCopy[i][j] = matOriginal[i][j];
        }
    }
}


//shuffle list
void shuffleList(double L1[N], double R1[N], double L2[N], double R2[N], double r[N]){
    int i=0;
    double temp1=0.0;
    double temp2=0.0;
    double temp3=0.0;
    double temp4=0.0;
    double temp5 = 0.0;
    int rval = 0;
    
    for (i=0; i<N; i++) {
        rval = rand()%N;
        //printf("rval %.d \n", rval);
        temp1 = L1[i];
        temp2 = R1[i];
        temp3 = L2[i];
        temp4 = R2[i];
        temp5 = r[i];
        
        
        L1[i] = L1[rval];
        R1[i] = R1[rval];
        L2[i] = L2[rval];
        R2[i] = R2[rval];
        r[i] = r[rval];
        
        L1[rval] = temp1;
        R1[rval] = temp2;
        L2[rval] = temp3;
        R2[rval] = temp4;
        r[rval] = temp5;
    }
}//shuffleList









