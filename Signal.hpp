class Signal
{
private:
	int x;
public:

    double bR[N];                       // receptor 1 production rate vector
    double bL[N];                       // ligand 1 production rate vector
    double br[N];                       // receptor 2 production rate vector
    double bl[N];                       // ligand 2 production rate vector
    double cells_SS[N][6];              // R*, L*, LR*, r*, l*, lr*
    double recomb[N];                    // recombination rate vector; shared recombination rate for all
    int mated_index[N];
    double s_Evolution[numt];           // evolution of s over time
    double r_Evolution[numt];           // evolution of recombination over time
    double si_ri_Evo[N][numt];          // evolution of si and ri for individual cells

    double FreeR(double b_R, double b_L, double k_on, double k_off, double g_L, double g_R, double g_LR);
    double FreeL(double b_R, double b_L, double k_on, double k_off, double g_L, double g_R, double g_LR);
    double FreeLR(double b_R, double b_L, double k_on, double k_off, double g_L, double g_R, double g_LR);
    double GetMaxK();
    double W12(double R1, double L1, double R2, double L2, double L1R1, double k_b, double n);
    double p_mate(double W, double K);
    void MATE(double bL1[N], double bR1[N], double bL2[N], double bR2[N], double r[N]);
    void MUTATE(double bL1[N], double bR1[N], double bL2[N], double bR2[N], double r[N]);
    void Normalise(double bL[N], double bR[N]);
    double Get_S(double bL1[N], double bR1[N], double bL2[N], double bR2[N]);
    double get_mean_r(double r_mat[N]);
    void recombination(double bsMated1[M][4], double bL[N], double bR[N], double bl[N], double br[N], double rMated[M], double r[N], double rM, int numMated, int num1, int num2);
    
    // int StabilityCheck(double fnow[N][N], double fbefore[N][N]);
};


