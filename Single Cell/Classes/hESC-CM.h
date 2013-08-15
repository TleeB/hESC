#ifndef HESCM_H
#define HESCM_H
#include <string>
#define numOfAPDs 100
class hESCM{

  public:
    hESCM(char);
    ~hESCM();
    void solve(double,double);
    double getV();
    std::string getName();
    void print_vars(bool, FILE*, double, int);

  private:
//Cell type identifier for switch statement
 char DevelopmentalStage;

//Scaling paramters and shifting parameters
 double RaINa;
 double myCoefTauH;
 double myCoefTauJ;
 double RaICaL;
 double Vh_dCa;
 double kCa;
 double ampICaLDinf;
 double KtaufCa;
 double myVhfCaL;
 double myKfCaL;
 double ampICaLFinf;
 double myKtauF;
 double myTauFShift;
 double myShiftFCaInf;
 double Vth_ICaL;
 double RaICaT;
 double RaINaCa;
 double alfa;
 double RaINaK;
 double RaIK1;
 double myShiftK1;
 double RaIKr;
 double mySlope;
 double myShift;
 double RaIKs;
 double RaIto;
 double myShiftItoR;
 double mySlopeItoR;
 double myShiftItoS;
 double mySlopeItoS;
 double myConstTauR;
 double myConstTauS;
 double RaIup;
 double RaIrel;
 double RaIleak;
 double RaIf;
 double myIfconst;
 double RaCm;
 double RaICap;
 double RaIKp;
 double RaIback;


//// Constants
 double R;
 double F;
 double T;
//// Buffering
 double Bufc;
 double Kbufc;
 double Bufsr;
 double Kbufsr;
//// Extracellular Ionic concentrations
//Commented because parameters are defined in Environment.h header
 //double Ko;
 //double Cao;
 //double Nao;
 double Nao_naca;
//// Intracellular Ionic concentrations
 double Ki;
 double Nai;
 double Vc;
 double Vsr;
 double capacitance;

//// Ionic Current Parameters

//// Fast Na+ Current 
 double Vh_h;
 double k_h;
 double Gnamax;

//// If Current 
 double Gf;

//// L-type Ca2+ Current 
 double GCaL;

//// T-type Ca2+ Current
 double GICaT;

//// Transient Outward Current 
 double GItoepi;
 double GItoendo;
 double GItoMcell;

//// IKs 
 double GKsepi;
 double GKsendo;
 double GKsMcell;
 double pKNa;

//// IKr 
 double GKr;
 double Q;
 double L0;
 double Kc;
 double Ka;

//// Inward Rectifier K+ Current
 double gK1max;

//// Na+/Ca2+ Exchanger Current 
 double knaca;
 double KmNai;
 double KmCa;
 double ksat;
 double n;

//// Na+/K+ Pump Current 
 double knak;
 double KmK;
 double KmNa;

//// IpCa
 double GpCa;
 double kpca;

//// IpK
 double GpK;

//// Background Currents 
 double GbNa;
 double GbCa;

///////////////////////////Calcium Dynamics////////////////////

//// Ileak
 double Vleak;

 double arel;
 double brel;
 double crel;

//// Iup
 double Vmaxup;
 double Kup;

 //Declaration of initial parameters
 double V;
 double Cai;
 double CaSR;
 double m;
 double h;
 double j;
 double xr1;
 double xr2;
 double xs;
 double r;
 double s;
 double d;
 double f;
 double f_ca;
 double g;
 double xf;
 double dCaT;
 double fCaT;

 //Gating parameters
 double alpha_K1, beta_K1, x_K1_inf;
 double alpha_m, beta_m, tau_m, m_inf, alpha_h, beta_h, tau_h, h_inf, alpha_j, beta_j, tau_j, j_inf;
 double alpha_d, beta_d, gamma_d, tau_d, d_inf, f_inf, tau_f, g_inf, tau_g, constg;
 double f_ca_inf, tau_f_ca, constf_ca, r_inf, tau_r, s_inf, tau_s, alpha_xs, beta_xs, xs_inf, tau_xs;
 double alpha_xr1, beta_xr1, xr1_inf, tau_xr1, alpha_xr2, beta_xr2, xr2_inf, tau_xr2, xf_inf, tau_xf;
 double dCaT_inf, tau_dCaT, fCaT_inf, tau_fCaT;
 double INa, ICaL, Ito, IKr, IKs, IK1, INaCa, INaK, IpCa, IpK, IbNa, IbCa, Istim,If,ICaT, Itotal;
 double Ileak, Iup, Irel, CaBuf, CaCSQN, CaCurrent, CaSRCurrent, dCaSR, bjsr, cjsr, dCai, bc, cc;

 //Nerst potentials
 double RTONF, Ek, Ena, Eks, Eca, Ef;

 //Coupling parameters
 double Igap_hes;
 int print_header;

 //// Flag to choose between epi, endo and M cell types
int epi;
int endo;
int Mcell;

int counter;
int start_output;
double vmin[numOfAPDs];
double tvmin[numOfAPDs];
double vmax[numOfAPDs];
double dvdtmax[numOfAPDs];
//double vdvdtmax[numOfAPDs];
double apd_start[numOfAPDs];
double apd_end[numOfAPDs];
double ddr[numOfAPDs];
double top[numOfAPDs];
double top_slope[numOfAPDs];
double apd30[numOfAPDs];
double apd50[numOfAPDs];
double apd70[numOfAPDs];
double apd90[numOfAPDs];
double bcl[numOfAPDs];

//APD paramters
double dvdt,dvdtold;
double vnew;

};

#endif

