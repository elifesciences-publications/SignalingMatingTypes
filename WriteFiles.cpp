#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"
#include "WriteFiles.hpp"


void  Write_Evolution(double mat[numt], int time,  int index, int rep){
    char name[550];
    int i=0;
    FILE *pfile=NULL;
    
    // use for cluster
    sprintf(name, "/Users/zenah/Dropbox/Signalling_sexes/C++_code/MT_self_signal/MT_self_signal/maxAlpha/Evo_N%.d_M%.d_mu%.3f_nuL_%.3f_nuR_%.3f_nur%.3f_rGlob%.2f_kon%.2f_koff%.2f_f%.2f_g%.2f_K%.1f_n%.1f_sigma%.2f_dv0%.2f_time%.d_%.d_rep%.d_alpha%.2f.txt", N, M, mu, nu_bL, nu_bR, nu_recomb, r_fixed, kon, koff, f,  gR, log10(K), n, sigma, dv0, time, index, rep, maxAlpha);
    //
    //sprintf(name, "/Users/zenah/Dropbox/Signalling_sexes/C++_code/MT_self_signal/MT_self_signal/Fig4b/Evo_N%.d_M%.d_mu%.3f_nuL_%.3f_nuR_%.3f_nur%.3f_rGlob%.2f_kon%.2f_koff%.2f_f%.2f_g%.2f_K%.1f_n%.1f_sigma%.2f_dv0%.2f_time%.d_%.d_rep%.d.txt", N, M, mu, nu_bL, nu_bR, nu_recomb, r_fixed, kon, koff, f,  gR, log10(K), n, sigma, dv0, time, index, rep);
    
    
    pfile = fopen(name, "w+");
    
    int min = time-20000;
    if(min < 0) min =0;
    
    for(i=min; i<time; i++){
        fprintf(pfile, "%.20f\n", mat[i]);
    }
    
    fclose(pfile);
}//WriteDelta

void  Write_Profile(double mat[N],  int time, int index, int rep){
    char name[550];
    int i=0;
    FILE *pfile=NULL;
    
    // use for desktop
    // sprintf(name, "/Users/zenah/Dropbox/Signalling_sexes/C++_code/SS_agent_based/sims_N%.d_M%.d_q%.3f_mu%.3f_sigma%.3f_nuL1_%.3f_nuR1_%.3f_nuL2_%.3f_nuR2_%.3f/Profile_kc%.5f_f%.3f_gR1%.5f_gL1%.5f_gR2%.5f_gL2%.5f_gg%.7f_%.d_time%.d.txt", N, M, q, mu, sigma, nu_bL1, nu_bR1, nu_bL2, nu_bR2, kc_global, f,  gR1, gL1, gR2, gL2, gg, index, time);
    // use for cluster
        sprintf(name, "/Users/zenah/Dropbox/Signalling_sexes/C++_code/MT_self_signal/MT_self_signal/maxAlpha/Profile_N%.d_M%.d_mu%.3f_nuL_%.3f_nuR_%.3f_nur%.3f_rGlob%.2f_kon%.2f_koff%.2f_f%.2f_g%.2f_K%.1f_n%.1f_time%.d_%.d_rep%.d_alpha%.2f.txt", N, M, mu, nu_bL, nu_bR, nu_recomb, r_fixed, kon, koff, f,  gR, log10(K), n, time, index, rep,maxAlpha);
    
    pfile = fopen(name, "w+");
    
    
    for(i=0; i<N; i++){
        fprintf(pfile, "%.20f\n", mat[i]);
    }
    
    fclose(pfile);
}//WriteDelta

void  Write_recomb(double mat[N],  int time, int rep){
    char name[550];
    int i=0;
    FILE *pfile=NULL;
    
    // use for desktop
    // sprintf(name, "/Users/zenah/Dropbox/Signalling_sexes/C++_code/SS_agent_based/sims_N%.d_M%.d_q%.3f_mu%.3f_sigma%.3f_nuL1_%.3f_nuR1_%.3f_nuL2_%.3f_nuR2_%.3f/Profile_kc%.5f_f%.3f_gR1%.5f_gL1%.5f_gR2%.5f_gL2%.5f_gg%.7f_%.d_time%.d.txt", N, M, q, mu, sigma, nu_bL1, nu_bR1, nu_bL2, nu_bR2, kc_global, f,  gR1, gL1, gR2, gL2, gg, index, time);
    // use for cluster
    sprintf(name, "r_N%.d_M%.d_mu%.3f_sigma%.2f_nuL_%.3f_nuR_%.3f_nur%.3f_rGlob%.2f_kon%.2f_koff%.2f_kb%.2f_f%.2f_g%.2f_norm%.2f_time%.d_rep%.d.txt", N, M, mu, sigma, nu_bL, nu_bR, nu_recomb, r_fixed, kon, koff, kb, f,  gR, normal_base, time, rep);
    
    pfile = fopen(name, "w+");
    
    
    for(i=0; i<N; i++){
        fprintf(pfile, "%.20f\n", mat[i]);
    }
    
    fclose(pfile);
}//WriteDelta
