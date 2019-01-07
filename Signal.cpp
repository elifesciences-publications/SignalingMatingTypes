#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"
#include "Signal.hpp"


/////////////////////////////////////////////////////////////////
//Here all the functions used in the object Signal are defined //
///////////////////////////////////////////////////////////////

// function for free receptors at steady state
double Signal::FreeR(double b_R, double b_L, double k_on, double k_off, double g_L, double g_R, double g_LR){
    double r = 0.0;
    double d=0.0;
    
    d = sqrt(pow(k_off*g_L*g_R + g_LR*(g_L*g_R + k_on*(b_L-b_R)),2.0)+4.0*k_on*g_L*g_R*g_LR*(k_off+g_LR)*b_R);
    
    r = (k_on*g_LR*(b_R-b_L) - k_off*g_L*g_R - g_L*g_R*g_LR + d)/(2.0*k_on*g_R*g_LR);
    
    return(r);
}//

// function for free ligand at steady state
double Signal::FreeL(double b_R, double b_L, double k_on, double k_off, double g_L, double g_R, double g_LR){
    double r = 0.0;
    double d=0.0;
    
    d = sqrt(pow(k_off*g_L*g_R + g_LR*(g_L*g_R + k_on*(b_L-b_R)),2.0)+4.0*k_on*g_L*g_R*g_LR*(k_off+g_LR)*b_R);
    
    r = (k_on*g_LR*(b_L-b_R) - k_off*g_L*g_R - g_L*g_R*g_LR + d)/(2.0*k_on*g_L*g_LR);
    
    return(r);
}//

// function for bound ligand-receptor at steady state
double Signal::FreeLR(double b_R, double b_L, double k_on, double k_off, double g_L, double g_R, double g_LR){
    double r = 0.0;
    double d=0.0;
    
    d = sqrt(pow(k_off*g_L*g_R + g_LR*(g_L*g_R + k_on*(b_L-b_R)),2.0)+4.0*k_on*g_L*g_R*g_LR*(k_off+g_LR)*b_R);
    
    r = (k_on*g_LR*(b_L+b_R) + k_off*g_L*g_R + g_L*g_R*g_LR - d)/(2.0*k_on*g_LR*g_LR);
    
    return(r);
}//

double Signal::GetMaxK(){
    double R1;
    double L1;
    double L1R1;
    double r1;
    double l1;
    double l1r1;
    double R2;
    double L2;
    double L2R2;
    double r2;
    double l2;
    double l2r2;
    double out = -1.0;
    
    R1 = FreeR(1, 0, kon, koff, gL, gR, gLR);
    L1 = FreeL(1, 0, kon, koff, gL, gR, gLR);
    L1R1 = FreeLR(1, 0, kon, koff, gL, gR, gLR);
    
    r1 = FreeR(0, 1, kon, koff, gL, gR, gLR);
    l1 = FreeL(0, 1, kon, koff, gL, gR, gLR);
    l1r1 = FreeLR(0, 1, kon, koff, gL, gR, gLR);
    
    R2 = FreeR(0, 1, kon, koff, gL, gR, gLR);
    L2 = FreeL(0, 1, kon, koff, gL, gR, gLR);
    L2R2 = FreeLR(0, 1, kon, koff, gL, gR, gLR);
    
    r2 = FreeR(1, 0, kon, koff, gL, gR, gLR);
    l2 = FreeL(1, 0, kon, koff, gL, gR, gLR);
    l2r2 = FreeLR(1, 0, kon, koff, gL, gR, gLR);
    
   // printf("R1: %.5f    L1: %.5f    R2: %.5f    L2: %.5f  \n", R1, L1, R2, L2);
   // printf("r1: %.5f    l1: %.5f    r2: %.5f    l2: %.5f  \n", r1, l1, r2, l2);
    
   // printf("w12: %.5f\n",W12(R1, L1, R2, L2, L1R1, kb, n) + W12(r1, l1, r2, l2, l1r1, kb, n));
   // printf("w12: %.5f\n",(W12(R2, L2, R1, L1, L2R2, kb, n))+W12(r2, l2, r1, l1, l2r2, kb, n));

    
    out = (W12(R1, L1, R2, L2, L1R1, kb, n)+W12(r1, l1, r2, l2, l1r1, kb, n))*(W12(R2, L2, R1, L1, L2R2, kb, n)+W12(r2, l2, r1, l1, l2r2, kb, n));
    
    return(out);
    
}


// computes the bidning strength between cells 1 and 2 with free receptor/ligands Ri, Li
double Signal::W12(double R1, double L1, double R2, double L2, double L1R1, double k_b, double n){
    double r = 0.0;
    
 
    if(L1R1 + k_b*R1*L2 < 0.0000001) r = 0.0;
    else r = (k_b*R1*L2)*pow(1.0-L1R1/(L1R1 + k_b*R1*L2),n);
    
    //printf("R1: %.5f    L1: %.5f    R2: %.5f    L2: %.5f  L1R1:%.5f     W12:%.3f \n\n", R1, L1, R2, L2, L1R1, r);
    
    return(r);
}

// returns the probability of mating given a binding strength bindStr
double Signal::p_mate(double W,  double K){
    //return(W*W/(K*K+W*W));
    return(W/(K+W));
}

void Signal::MATE(double bL[N], double bR[N], double bl[N], double br[N], double r[N]){
    int num1=0;
    int num2=0;
    int i1=0;
    int i2=0;
    int i = 0;
    int j = 0;
    int numMated = 0;
    double W = 0.0;
    double prob_mate = 0.0;
    double prob_r = 0.0;
    double rM = 0.0;
    
    int mat[N];
    double RLrl[2][6];
    double bsMated[M][4];   // 0-3 for bL1, bR1, bL2, bR2 and M for number of mated cells
    double rMated[M];       // coupled to bsMated; ith entry for recombination of ith mated cell
    
    //initiate matrices
    for(i=0; i<N; i++) mat[i] = i;
    for(i=0; i<M; i++) rMated[i] = -1.0;
    for(i=0; i<2; i++) for(j=0; j<4; j++) RLrl[i][j] =0.0;
    for(i=0; i<M; i++) for(j=0; j<4; j++) bsMated[i][j] =0.0;
    
    //randomly choose two numbers from 0 to N - 1 (must be different)
    while(numMated < M){
        // randomly choose two cells that have not mated yet and that are different from one another
        num1=-1;
        num2=-1;
        while(num1<0){
            i1 = uniform_distribution(0, N-1);
            num1 = mat[i1];
        }
        while(num2<0 || num2==num1){
            i2 = uniform_distribution(0, N-1);
            num2 = mat[i2];
        }
        //find R, L, LR, r, l, lr at SS for two cells
        // cell 1
        RLrl[0][0] = FreeR(bR[num1], bL[num1], kon, koff, gL, gR, gLR); //R
        RLrl[0][1] = FreeL(bR[num1], bL[num1], kon, koff, gL, gR, gLR); //L
        RLrl[0][2] = FreeLR(bR[num1], bL[num1], kon, koff, gL, gR, gLR); //LR
        RLrl[0][3] = FreeR(br[num1], bl[num1], kon, koff, gL, gR, gLR); //r
        RLrl[0][4] = FreeL(br[num1], bl[num1], kon, koff, gL, gR, gLR); //l
        RLrl[0][5] = FreeLR(br[num1], bl[num1], kon, koff, gL, gR, gLR); //lr
       // printf("Cell 1: R:%.3f L%.3f  RL:%.3f \n", RLrl[0][0], RLrl[0][1], RLrl[0][2]);
        // cell 2
        RLrl[1][0] = FreeR(bR[num2], bL[num2], kon, koff, gL, gR, gLR); //R
        RLrl[1][1] = FreeL(bR[num2], bL[num2], kon, koff, gL, gR, gLR); //L
        RLrl[1][2] = FreeLR(bR[num2], bL[num2], kon, koff, gL, gR, gLR); //LR
        RLrl[1][3] = FreeR(br[num2], bl[num2], kon, koff, gL, gR, gLR); //r
        RLrl[1][4] = FreeL(br[num2], bl[num2], kon, koff, gL, gR, gLR); //l
        RLrl[1][5] = FreeLR(br[num2], bl[num2], kon, koff, gL, gR, gLR); //lr
        ///////////////////////////
        //find the prob they mate//
        //////////////////////////
        // first find bdinding strength
        W = (W12(RLrl[0][0], RLrl[0][1], RLrl[1][0], RLrl[1][1], RLrl[0][2], kb, n)+ W12(RLrl[0][3], RLrl[0][4], RLrl[1][3], RLrl[1][4], RLrl[0][5], kb, n))* (W12(RLrl[1][0], RLrl[1][1], RLrl[0][0], RLrl[0][1], RLrl[1][2], kb, n)+ W12(RLrl[1][3], RLrl[1][4], RLrl[0][3], RLrl[0][4], RLrl[1][5], kb, n));
        // the pmate
        prob_mate = p_mate(W, K);
        if(unifRand()<prob_mate){ //then mate
         //       printf("i:%.d   j:%.d   W:%.5f  pMate:%.5f\n\n\n", num1, num2, W, prob_mate);
            //set cells to -1
            mat[i1] = -1;
            mat[i2] = -1;
            /////////////////////////////////////////////////
            // RECOMBINATION STEP //
            // store mated pairs and implement recombination
//            // PAIR 1 //
//            // offspring 1 pair 1
//            bsMated[numMated][0] = bL[num1];
//            bsMated[numMated][1] = bR[num1];
//            // offspring 2 pair 1
//            bsMated[numMated+1][0] = bL[num2];
//            bsMated[numMated+1][1] = bR[num2];
//            // PAIR 2 //
//            // offpsirng 1 pair 2
//            bsMated[numMated][2] = bl[num1];
//            bsMated[numMated][3] = br[num1];
//            // offpsirng 2 pair 2
//            bsMated[numMated+1][2] = bl[num2];
//            bsMated[numMated+1][3] = br[num2];
            recombination(bsMated, bL, bR, bl, br, rMated, r, 0.5*(r[num1]+r[num2]), numMated, num1, num2);
            //increment number of mated cells
            numMated = numMated + 2;
        }//if
    }

    // once M cells have mated sample with replacement from M mated cells back to N
    int rval = 0;
    for(i=0; i<N; i++){
        rval = rand()%M;
        bL[i] = bsMated[rval][0];
        bR[i] = bsMated[rval][1];
        bl[i] = bsMated[rval][2];
        br[i] = bsMated[rval][3];
        r[i] = rMated[rval];
    }
}// MATE

void Signal::recombination(double bsMated1[M][4], double bL[N], double bR[N], double bl[N], double br[N], double rMated[M], double r[N], double rM, int numMated, int num1, int num2){
    double p=0.0;
    
    p = unifRand();
    
    if(p<= (1-rM)*(1-rM)){ // no recombination
        // offspring 1 RL
        bsMated1[numMated][0] = bL[num1];
        bsMated1[numMated][1] = bR[num1];
        // offpsirng 1 rl
        bsMated1[numMated][2] = bl[num1];
        bsMated1[numMated][3] = br[num1];
        // offspring 1 recombination rate
        rMated[numMated] = r[num1];
        //////////////////////////////////////////////
        // offspring 2 rl
        bsMated1[numMated+1][0] = bL[num2];
        bsMated1[numMated+1][1] = bR[num2];
        // offpsirng 2 rl
        bsMated1[numMated+1][2] = bl[num2];
        bsMated1[numMated+1][3] = br[num2];
        // offspring 1 recombination rate
        rMated[numMated+1] = r[num2];
        
    }
    else if(p > (1-rM)*(1-rM) && p <= (1-rM)*(1-rM) + rM*rM){ // two recombination events, swap rMs
        // offspring 1 RL
        bsMated1[numMated][0] = bL[num1];
        bsMated1[numMated][1] = bR[num1];
        // offpsirng 1 rl
        bsMated1[numMated][2] = bl[num1];
        bsMated1[numMated][3] = br[num1];
        // offspring 1 recombination rate
        rMated[numMated] = r[num2];
        //////////////////////////////////////////////
        // offspring 2 rl
        bsMated1[numMated+1][0] = bL[num2];
        bsMated1[numMated+1][1] = bR[num2];
        // offpsirng 2 rl
        bsMated1[numMated+1][2] = bl[num2];
        bsMated1[numMated+1][3] = br[num2];
        // offspring 1 recombination rate
        rMated[numMated+1] = r[num1];
    }
    else if(p > (1-rM)*(1-rM) + rM*rM && p <= (1-rM)*(1-rM) + rM*rM + rM*(1-rM) ){ // one recombination event, swap rM-Ll
        // offspring 1 RL
        bsMated1[numMated][0] = bL[num2];
        bsMated1[numMated][1] = bR[num1];
        // offpsirng 1 rl
        bsMated1[numMated][2] = bl[num2];
        bsMated1[numMated][3] = br[num1];
        // offspring 1 recombination rate
        rMated[numMated] = r[num2];
        //////////////////////////////////////////////
        // offspring 2 rl
        bsMated1[numMated+1][0] = bL[num1];
        bsMated1[numMated+1][1] = bR[num2];
        // offpsirng 2 rl
        bsMated1[numMated+1][2] = bl[num1];
        bsMated1[numMated+1][3] = br[num2];
        // offspring 1 recombination rate
        rMated[numMated+1] = r[num1];
    }
    else{ // one recombination event,swap rM-Rr
        // offspring 1 RL
        bsMated1[numMated][0] = bL[num1];
        bsMated1[numMated][1] = bR[num2];
        // offpsirng 1 rl
        bsMated1[numMated][2] = bl[num1];
        bsMated1[numMated][3] = br[num2];
        // offspring 1 recombination rate
        rMated[numMated] = r[num2];
        //////////////////////////////////////////////
        // offspring 2 rl
        bsMated1[numMated+1][0] = bL[num2];
        bsMated1[numMated+1][1] = bR[num1];
        // offpsirng 2 rl
        bsMated1[numMated+1][2] = bl[num2];
        bsMated1[numMated+1][3] = br[num1];
        // offspring 1 recombination rate
        rMated[numMated+1] = r[num1];
    }
}//recombination

//void Signal::MUTATE(double bL[N], double bR[N], double bl[N], double br[N], double r[N]){
//    int i=0;
//    double e = 0.0;
//    
//    for(i=0; i<N; i++){
//        if(unifRand()<nu_bL){
//            e = unifRand();
//            if(unifRand()>0.5) e = e;
//            else e = -e;
//            bL[i] = bL[i] - e;
//            bl[i] = bl[i] + e;
//            if(bL[i]<0.0) bL[i] = 0.0;
//            if(bL[i]>1.0) bL[i] = 1.0;
//            if(bl[i]<0.0) bl[i] = 0.0;
//            if(bl[i]>1.0) bl[i] = 1.0;
//        }
//         if(unifRand()<nu_bR){
//             e = unifRand();
//             if(unifRand()>0.5) e = e;
//             else e = -e;
//             bR[i] = bR[i] - e;
//             br[i] = br[i] + e;
//             if(bR[i]<0.0) bR[i] = 0.0;
//             if(bR[i]>1.0) bR[i] = 1.0;
//             if(br[i]<0.0) br[i] = 0.0;
//             if(br[i]>1.0) br[i] = 1.0;
//         }
//        // mutate recombindation locues
//        if(unifRand()<nu_recomb){
//            r[i] = r[i] + rand_Normal(mu_r, sigma_r);
//            if(r[i]<0.0) r[i] = 0.0;
//            if(r[i]>0.5) r[i] = 0.5;
//        }
//    }
//}

void Signal::MUTATE(double bL[N], double bR[N], double bl[N], double br[N], double r[N]){
    int i=0;
    double e = 0.0;
    
    for(i=0; i<N; i++){
        // mutate L
        if(unifRand()<nu_bL){
            e = rand_Normal(0.0, sigma);
            bL[i] = bL[i] + e;
            if(bL[i]<0.0) bL[i] = bL[i] - e;
            if(bL[i]>1.0) bL[i] = bL[i] - e;
        }
        // mutate l
        if(unifRand()<nu_bL){
            e = rand_Normal(0.0, sigma);
            bl[i] = bl[i] + e;
            if(bl[i]<0.0) bl[i] = bl[i] - e;
            if(bl[i]>1.0) bl[i] = bl[i] - e;
        }
        // max at maxAlpha
        if (bL[i] + bl[i]>maxAlpha) {
            bL[i] = maxAlpha*bL[i]/(bL[i]+bl[i]);
            bl[i] = maxAlpha*bl[i]/(bL[i]+bl[i]);
        }
        
        // mutate R
        if(unifRand()<nu_bR){
            e = rand_Normal(0.0, sigma);
            bR[i] = bR[i] + e;
            if(bR[i]<0.0) bR[i] = bR[i] - e;
            if(bR[i]>1.0) bR[i] = bR[i] - e;
        }
        // mutate r
        if(unifRand()<nu_bR){
            e = rand_Normal(0.0, sigma);
            br[i] = br[i] + e;
            if(br[i]<0.0) br[i] = br[i] - e;
            if(br[i]>1.0) br[i] = br[i] - e;
        }
        // max at maxAlpha
        if (bR[i] + br[i]>maxAlpha) {
            bR[i] = maxAlpha*bR[i]/(bR[i]+br[i]);
            br[i] = maxAlpha*br[i]/(bR[i]+br[i]);
        }
        // mutate recombindation locues
        if(unifRand()<nu_recomb){
            e = rand_Normal(0.0, sigma_r);
            r[i] = r[i] + e;
            if(r[i]<0.0) r[i] = 0.0 - e;
            if(r[i]>0.5) r[i] = 0.5 - e;
        }
    }
}

void Signal::Normalise(double bL[N], double bR[N]){
    int i=0;
    
    for(i=0; i<N; i++){
        //printf("divide by sum: %.5f\n", bL[i]+bR[i]);
        if(bL[i]+bR[i]>normal_base){
            bL[i] = normal_base*bL[i]/(bL[i]+bR[i]);
            bR[i] = normal_base*bR[i]/(bL[i]+bR[i]);
        }
    }
    
}

double Signal::Get_S(double bL1[N], double bR1[N], double bL2[N], double bR2[N]){
    int i = 0;
    double sum = 0.0;
    
    for(i=0; i<N; i++){
        sum = sum + Absolute(bL1[i], bR1[i]) + Absolute(bL2[i], bR2[i]);
    }
    return(sum/(2.0*double(N)));
}

double Signal::get_mean_r(double r_mat[N]){
    int i=0;
    double sum = 0.0;
    
    for(i=0; i<N; i++) sum = sum + r_mat[i];
    
    return(sum/(double N));
}




