#include <cmath>

//#define EARLY
#define LATE
struct ode_params
{
 char celltype;
 double RTONF,R,T,F, Ek, Ena, Eks, Eca, Ef;
 double V, Cai, CaSR, alpha_K1, beta_K1, x_K1_inf, shift;
 double alpha_m, beta_m, tau_m, m_inf, alpha_h, beta_h, tau_h, h_inf, alpha_j, beta_j, tau_j, j_inf, alpha_d, beta_d, gamma_d, tau_d, d_inf, f_inf, tau_f, g_inf, tau_g, constg, f_ca_inf, tau_f_ca, constf_ca, r_inf, tau_r, s_inf, tau_s, alpha_xs, beta_xs, xs_inf, tau_xs, alpha_xr1, beta_xr1, xr1_inf, tau_xr1, alpha_xr2, beta_xr2, xr2_inf, tau_xr2, xf_inf, tau_xf, dCaT_inf, tau_dCaT, fCaT_inf, tau_fCaT;
 double INa, m, h, j, ICaL, d, f, g, f_ca, Ito, r, s, IKr, xr1, v_half, xr2, IKs, xs, IK1, INaCa, INaK, IpCa, IpK, IbNa, IbCa, Istim, If, xf, ICaT, dCaT, fCaT, Itotal;
 double Ileak, Iup, Irel,CaCurrent, CaSRCurrent, Caibufc, CaSRbufsr;
 double Ko, Nao, Cao, Nai, Ki, Vc, Vsr, Bufc, Kbufc,Bufsr,Kbufsr,Vmaxup,Kup,CAP,GNa,	GCaL,	pKNa,Gks ,GK1, knaca, KmNai, KmCa, ksat, n, KmK, KmNa, knak, GpCa, GpK, GbNa, GbCa, Vleak, Q, Gf,arel,brel,crel; 
 double RaINa,myCoefTauH , myCoefTauJ,RaICaL,Vh_dCa,kCa, ampICaLDinf , KtaufCa, myVhfCaL , myKfCaL, ampICaLFinf , myKtauF , myTauFShift ;
 double myShiftFCaInf ,Vth_ICaL , RaICaT , RaINaCa, alfa, RaINaK,RaIK1, myShiftK1 , RaIKr,poffi , mySlope , myShift , RaIKs;
  double RaIto;
  double myShiftItoR ;
  double mySlopeItoR ;
  double myShiftItoS ;
  double mySlopeItoS ;
  double myConstTauR;
  double myConstTauS;
  double RaIup;
  double RaIrel;
  double RaIleak;
  double RaIf;
  double Rax0;
  double Radx;
  double myIfconst ;
  double RaCm;
  double RaICap;
  double RaIKp;
  double RaIback;
  
 double Nao_naca;
 double capacitance;
 double svolt;
 double Vh_h;
 double k_h;
 double Gnamax;
 double GICaT ;
 double GItoepi;
 double GItoendo;
 double GItoMcell;
 double GKsepi   ;
 double GKsendo  ;
 double GKsMcell ;
 double GKr;
 double gK1max;
 double L0;
 double Kc;
 double Ka;
 double kpca;

};
/*
#ifdef EARLY
char celltype = 'e';
double RaINa=0.038;
double myCoefTauH =2.8;
double myCoefTauJ = 1;
double RaICaL=0.25;
double Vh_dCa=12.5;
double kCa=12.5;
double ampICaLDinf = 1;
double KtaufCa=1433;   // [1/mM]  Altamirano & Bers 2007
double myVhfCaL = 20;
double myKfCaL = 7;//7
double ampICaLFinf = 1;//to accelerate the tauF and to inactivate more rapidly Ical
double myKtauF = 1;
double myTauFShift = 0;
double myShiftFCaInf = -0.11e-3;  //[mM] to shift a litte bit the sigmoidal of fCa
double Vth_ICaL = -0.060;
double RaICaT =0.25;
double RaINaCa=1.750e1;
double alfa=0.8;
double RaINaK=0.7;
double RaIK1=0.05*2.67/3; 
double myShiftK1 = -0.015;
double RaIKr=3;
double poffi = 0;
double mySlope = 1;
double myShift = 0;
double RaIKs=0.1;
double RaIto=0.1673*0.4903*0.8;
double myShiftItoR = -25;
double mySlopeItoR = 0.3;
double myShiftItoS = 0;
double mySlopeItoS = 1;
double myConstTauR= 1;
double myConstTauS= 1;
double RaIup=0.4/3;
double RaIrel=0.2/18;
double RaIleak=0.1/18;
double RaIf=0.5389;
double Rax0=1.1061;
double Radx=0.6537;
double myIfconst = 1;
double RaCm=0.22162;
double RaICap=1;
double RaIKp=0;
double RaIback=0.2;
#endif

#ifdef LATE
char celltype='l';
double RaINa=1;    // Itoh
double myCoefTauH = 2.8;
double myCoefTauJ = 1;
double RaICaL=0.422;
double Vh_dCa=16;
double kCa=12.8;
double ampICaLDinf = 1;
double KtaufCa=1433;   // [1/mM]  Altamirano & Bers 2007
double myVhfCaL = 20;
double myKfCaL = 7;
double ampICaLFinf = 1;
double myKtauF = 1;
double myTauFShift = 0;
double myShiftFCaInf = -0.12e-3; //[mM] to move a bit left the sigmoid  fCa
double Vth_ICaL = 0;
double RaICaT =0.05;
double RaINaCa=1.824e+01;
double alfa=0.38;
double RaINaK=0.83; 
double RaIK1=0.4*0.05*2.67*4;
double myShiftK1 = -0.015;
double RaIKr=1.4;
double poffi = 0;
double mySlope = 1;
double myShift = 0;
double RaIKs=0.1;
double RaIto=0.3754*0.4903*0.9;
double myShiftItoR = -25;
double mySlopeItoR = 0.3;
double myShiftItoS = 0;
double mySlopeItoS = 1;
double myConstTauR= 1;
double myConstTauS= 1;
double RaIup=0.33;
double RaIrel=0.4;
double RaIleak=0.3*1;
double RaIf=0.23;
double Rax0=1.1061;
double Radx=0.6537;
double RaCm=0.17838;
double RaICap=1;
double RaIKp=0;
double RaIback=1;
#endif

//// Constants
double R=8.314472;   // [J/mol/K] Gas constant
double F=96485.3415; // [C/mol]	  Faraday constant
double T=310.0;      // [K]       Temperature 
//// Buffering
double Bufc=0.25;   // [mM] total cytoplasmic buffer concentration
double Kbufc=0.001;  // [mM] Cai half saturation constant
double Bufsr=10;     // [mM] total sarcoplasmic buffer concentration
double Kbufsr=0.3;   // [mM] CaSR half saturation constant
//// Extracellular Ionic concentrations
double Ko=4;         // [mM]
double Cao=2.7;      // [mM]
double Nao=150.5;    // [mM]
double Nao_naca = Nao;
//// Intracellular Ionic concentrations
double Ki=140;      // [mM]  & 140 in TT04
double Nai=7;    // [mM]
double Vc=0.016404*RaCm;
double Vsr=0.001094*RaCm;
double capacitance=0.185*1000*RaCm;
//// Resting Potential
double svolt=-70;  // [mV]
//// Ionic Currents

//// Fast Na+ Current 
double Vh_h=-73;
double k_h=5.6;
double Gnamax=14838*RaINa; // [S/F] maximal INa conductance

//// If Current 
double Gf=90.926*RaIf;
double x0=-89.73015*Rax0;
double dx=11.7335*Radx;

//// L-type Ca2+ Current 
double GCaL=0.000175*RaICaL;  // [m^3/F/s] maximal ICaL conductance

//// T-type Ca2+ Current
double GICaT = 183.2*RaICaT; //[S/F] 

//// Transient Outward Current 
double GItoepi=294*RaIto;   // [S/F] maximal ITo conductance
double GItoendo=73;   // [S/F] maximal ITo conductance
double GItoMcell=294; // [S/F] maximal ITo conductance

//// IKs 
double GKsepi   =157*RaIKs; //245; //[S/F] maximal IKs conductance
double GKsendo  =157; //245;// [S/F] maximal IKs conductance
double GKsMcell =40; //62;// [S/F] maximal IKs conductance
double pKNa=0.03;   // [ ]

//// IKr 
double GKr=96*sqrt(Ko/5.4)*RaIKr; //GKr=96 S/F maximal IKr conductance
double xr1o=0;
double xr2o=1;
double Q=2.3;
double L0=0.025;
double Kc=0.58e-3;
double Ka=2.6e-3;

//// Inward Rectifier K+ Current
double gK1max=5405*sqrt(Ko/5.4)*RaIK1; //GK1=5405 S/F maximal IK1 conductance

//// Na+/Ca2+ Exchanger Current 
double knaca=1000*RaINaCa;  // [A/F] maximal INaCa
double KmNai=87.5;  // [mM] Nai half saturation constant 
double KmCa=1.38;   // [mM] Cai half saturation constant 
double ksat=0.1;    // [ ]  saturation factor for INaCa 
double n=0.35;      // [ ]  voltage dependence parameter

//// Na+/K+ Pump Current 
double knak=1.362*RaINaK;  // [A/F] maximal INaK
double KmK=1;       // [mM] Ko half saturation constant 
double KmNa=40;     // [mM] Nai half saturation constant 

//// IpCa
double GpCa=0.825*RaICap;    // [A/F] maximal IpCa 
double kpca=0.0005;   // [mM]  Cai half saturation constant 

//// IpK
double GpK=14.6*RaIKp;    // [S/F] maximal IpK conductance

//// Background Currents 
double GbNa=0*0.29*RaIback;   // [S/F] maximal IbNa conductance
double GbCa=0.592*RaIback;  // [S/F] maximal IbCa conductance

///////////////////////////Calcium Dynamics////////////////////

//// Ileak
double Vleak=0.08*RaIleak; // [1/s] maximal Ileak =0.00008/s

//// Irel
double arel=16.464;   // [mM/s] maximal CaSR-dependent Irel
double brel=0.25;     // [mM] CaSR half saturation constant
double crel=8.232;    // [mM/s] maximal CaSR-independent Irel

//// Iup
double Vmaxup=0.425*RaIup; //[mM/s]   // 0.000425;   // [mM/ms] maximal Iup 
double Kup=0.00025;//0.00025;    // [mM] half saturation constant
*/
