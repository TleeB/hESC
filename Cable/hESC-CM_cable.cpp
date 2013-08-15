#include <stdio.h>
#include <math.h>
#include <iostream>

#define FILENAME 20 //Array size holding file name
//Parameters currently only represent epicardial cells, other cell types will be implemented at
//later date

using namespace std;

//timestep
double DT = 0.01;//Time converted to ms
double tmin = 0.0;
double tmax = 10000;

char output_file_name[ FILENAME ]= "cbl_adult.dat";

//Cable parameters
int NODES = 100;//number of nodes
double D = 0.001;//cm^2/ms
double DX  =  0.025;//cm

//Char to define the developmental stage of the hESC of the cable
char DevelopmentalStage = 'a';

double RaINa;
double myCoefTauH;
double myCoefTauJ;
double RaICaL;
double Vh_dCa;
double kCa;
double ampICaLDinf;
double KtaufCa;   // [1/mM]  Altamirano & Bers 2007
double myVhfCaL;
double myKfCaL;
double ampICaLFinf;
double myKtauF ;
double myTauFShift;
double myShiftFCaInf ;
double Vth_ICaL;
double RaICaT;
double RaINaCa;
double alfa;
double RaINaK;
double RaIK1;
double myShiftK1;
double RaIKr;
double poffi;
double mySlope;
double myShift ;
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
double Rax0;
double Radx;
double myIfconst;
double RaCm;
double RaICap;
double RaIKp;
double RaIback;

int main(){
//Switch case for the developmental stage 'e' = early, 'l' = late and 'a' = adult, where adult represents the modifier tentusscher model
  switch (DevelopmentalStage){
    case 'e':
       RaINa = 0.038;
       myCoefTauH =2.8;
       myCoefTauJ = 1;
       RaICaL = 0.25;
       Vh_dCa = 12.5;
       kCa = 12.5;
       ampICaLDinf = 1;
       KtaufCa = 1433;   // [1/mM]  Altamirano & Bers 2007
       myVhfCaL = 20;
       myKfCaL = 7;
       ampICaLFinf = 1;
       myKtauF = 1;
       myTauFShift = 0;
       myShiftFCaInf = -0.11e-3;
       Vth_ICaL = -0.060;
       RaICaT =0.25;
       RaINaCa = 1.750e1;
       alfa = 0.8;
       RaINaK = 0.7;
       RaIK1 = 0.05*2.67/3;
       myShiftK1  =  -15;
       RaIKr = 3;
       poffi  =  0;
       mySlope  =  1;
       myShift  =  0;
       RaIKs = 0.1;
       RaIto = 0.1673*0.4903*0.8;
       myShiftItoR  =  -25;
       mySlopeItoR  =  0.3;
       myShiftItoS  =  0;
       mySlopeItoS  =  1;
       myConstTauR =  1;
       myConstTauS =  1;
       RaIup = 0.4/3;
       RaIrel = 0.2/18;
       RaIleak = 0.1/18;
       RaIf = 0.5389;
       Rax0 = 1.1061;
       Radx = 0.6537;
       myIfconst  =  1;
       RaCm = 0.22162;
       RaICap = 1;
       RaIKp = 0;//In adult cells only
       RaIback = 0.2;
      break;
    case 'l':
       RaINa = 1;    // Itoh
       myCoefTauH  =  2.8;
       myCoefTauJ  =  1;
       RaICaL = 0.422;
       Vh_dCa = 16;
       kCa = 12.8;
       ampICaLDinf  =  1;
       KtaufCa = 1433;   // [1/mM]  Altamirano & Bers 2007
       myVhfCaL  =  20;
       myKfCaL  =  7;
       ampICaLFinf  =  1;
       myKtauF  =  1;
       myTauFShift  =  0;
       myShiftFCaInf  =  -0.12e-3; //[mM] per spostare un po'  a sx la signmoide di fCa
       Vth_ICaL  =  0;
       RaICaT  = 0.05;
       RaINaCa = 1.824e+01;
       alfa = 0.38;
       RaINaK = 0.83;
       RaIK1 = 0.4*0.05*2.67*4;
       myShiftK1  =  -15;
       RaIKr = 1.4;
       poffi  =  0;
       mySlope  =  1;
       myShift  =  0;
       RaIKs = 0.1;
       RaIto = 0.3754*0.4903*0.9;
       myShiftItoR  =  -25;
       mySlopeItoR  =  0.3;
       myShiftItoS  =  0;
       mySlopeItoS  =  1;
       myConstTauR =  1;
       myConstTauS =  1;
       RaIup = 0.33;
       RaIrel = 0.4;
       RaIleak = 0.3*1;
       RaIf = 0.23;
       Rax0 = 1.1061;
       Radx = 0.6537;
       RaCm = 0.17838;
       RaICap = 1;
       RaIKp = 0;//In adult cells only
       RaIback = 1;  // Itoh
      break;
    case 'a':
       RaINa = 1;
       myCoefTauH  =  1;
       myCoefTauJ  =  1;
       RaICaL = 1;
       Vh_dCa = -5;
       kCa = 7.5;
       RaINaCa = 1;
       alfa = 2.5;
       RaIK1 = 1;
       myShiftK1  =  0;
       RaIKr = 1;
       mySlope  =  1;
       myShift  =  0;
       RaIKs = 1;
       RaIto = 1;
       RaIup = 1;
       RaIrel = 1;
       RaIleak = 1;
       RaICaT = 0;//added, this scaling parameter was left out in the Paci model for the adult phenotype
       RaIf = 0;
       Rax0 = 1;
       Radx = 1;
       RaCm = 1;
       RaINaK = 1;
       RaICap = 1;
       RaIKp = 1;
       RaIback = 1;
      break;
    default: cout <<"No developmental stage selected."<<endl;
             break;
  }
  //// Constants
  double R = 8314.472;   // [J/millimole/K] Gas constant
  double F = 96485.3415; // [C/mol]	  Faraday constant
  double T = 310.0;      // [K]       Temperature

  //// Buffering
  double Bufc = 0.25;   // [mM] total cytoplasmic buffer concentration
  double Kbufc = 0.001;  // [mM] Cai half saturation constant
  double Bufsr = 10;     // [mM] total sarcoplasmic buffer concentration
  double Kbufsr = 0.3;   // [mM] CaSR half saturation constant

  //// Extracellular Ionic concentrations
  ////TT
  //Ko = 5.4;      // [mM]
  //Cao = 1.8;     // [mM]
  //Nao = 140;     // [mM]

  ////Sartiani
  double Ko = 4;         // [mM]
  double Cao = 2.7;      // [mM]
  double Nao = 150.5;    // [mM]
  double Nao_naca  =  Nao;

  //// Intracellular Ionic concentrations
  // // Pre-dialysis
  double Ki = 140;      // [mM]  & 140 in TT04
  double Nai = 7;    // [mM]

  //// Intracellular Volumes
  double Vc = 16.404*RaCm;
  double Vsr = 1.094*RaCm;
  double capacitance = 0.185*1000*RaCm;

  //// Flag to choose between epi, endo and M cell types
  int epi = 1;
  int endo = 0;
  int Mcell = 0;

  //// Ionic Currents

  //// Fast Na+ Current
  double Vh_h = -73;
  double k_h = 5.6;
  double Gnamax = 14.838*RaINa; // [nS/pF] maximal INa conductance

  //// If Current
  double Gf = 0.090926*RaIf;
  double x0 = -89.73015*Rax0;
  double dx = 11.7335*Radx;

  //// L-type Ca2+ Current
  double GCaL = 0.000175*RaICaL;  // [m^3/F/s] maximal ICaL conductance

  //// T-type Ca2+ Current
  double GICaT  =  0.1832*RaICaT; //[nS/pF]

  //// Transient Outward Current
  double GItoepi = 0.294*RaIto;   // [nS/pF] maximal ITo conductance
  double GItoendo = 73;   // [nS/pF] maximal ITo conductance
  double GItoMcell = 294; // [nS/pF] maximal ITo conductance
  int soepi = 1;
  int soendo = 1;

  //// IKs
  double GKsepi    = 0.157*RaIKs; //245; //[S/F] maximal IKs conductance
  double GKsendo   = 157; //245;// [nS/pF] maximal IKs conductance
  double GKsMcell  = 40; //62;// [nS/pF] maximal IKs conductance
  double pKNa = 0.03;   // [ ]

  //// IKr
  double GKr = 0.096*sqrt(Ko/5.4)*RaIKr; //GKr = 96 S/F maximal IKr conductance
  double Q = 2.3;
  double L0 = 0.025;
  double Kc = 0.58e-3;
  double Ka = 2.6e-3;

  //// Inward Rectifier K+ Current
  double gK1max = 5.405*sqrt(Ko/5.4)*RaIK1; //GK1 = 5405 S/F maximal IK1 conductance

  //// Na+/Ca2+ Exchanger Current
  double knaca = 1000*RaINaCa;  // [pA/pF] maximal INaCa
  double KmNai = 87.5;  // [mM] Nai half saturation constant
  double KmCa = 1.38;   // [mM] Cai half saturation constant
  double ksat = 0.1;    // [ ]  saturation factor for INaCa
  double n = 0.35;      // [ ]  voltage dependence parameter

  //// Na+/K+ Pump Current
  double knak = 1.362*RaINaK;  // [pA/pF] maximal INaK
  double KmK = 1;       // [mM] Ko half saturation constant
  double KmNa = 40;     // [mM] Nai half saturation constant

  //// IpCa
  double GpCa = 0.000825*RaICap;    // [pA/pF] maximal IpCa
  double kpca = 0.0005;   // [mM]  Cai half saturation constant

  //// IpK
  double GpK = 0.0146*RaIKp;    // [nS/pF] maximal IpK conductance

  //// Background Currents
  double GbNa = 0*0.00029*RaIback;   // [nS/pF] maximal IbNa conductance
  double GbCa = 0.000592*RaIback;  // [nS/pF] maximal IbCa conductance

  //// Calcium Dynamics
  //// Ileak
  double Vleak = 0.00008*RaIleak; // [1/s] maximal Ileak  = 0.00008/s

  //// Irel
  double arel = 0.016464;   // [mM/s] maximal CaSR-dependent Irel
  double brel = 0.25;     // [mM] CaSR half saturation constant
  double crel = 0.008232;    // [mM/s] maximal CaSR-independent Irel

  //// Iup
  double Vmaxup = 0.000425*RaIup; //[mM/s]   // 0.000425;   // [mM/ms] maximal Iup
  double Kup = 0.00025;//0.00025;    // [mM] half saturation constant



  FILE *output;
  FILE *outputCV;
  double CVstart;
  outputCV  =  fopen("CV.dat","w");
  fclose(outputCV);
  int file_output_counter = 0;
  double RTONF, Ek, Ena, Eks, Eca, Ef;

  double V[NODES], Vn[NODES], Vp[NODES], Cai[NODES], CaSR[NODES], alpha_K1, beta_K1, x_K1_inf;
  double alpha_m, beta_m, tau_m, m_inf, alpha_h, beta_h, tau_h, h_inf, alpha_j, beta_j, tau_j, j_inf, alpha_d, beta_d, gamma_d, tau_d, d_inf, f_inf, tau_f, g_inf, tau_g, constg, f_ca_inf, tau_f_ca, constf_ca, r_inf, tau_r, s_inf, tau_s, alpha_xs, beta_xs, xs_inf, tau_xs, alpha_xr1, beta_xr1, xr1_inf, tau_xr1, alpha_xr2, beta_xr2, xr2_inf, tau_xr2, xf_inf, tau_xf, dCaT_inf, tau_dCaT, fCaT_inf, tau_fCaT;
  double INa, m[NODES], h[NODES], j[NODES], ICaL, d[NODES], df, tempf, f[NODES], g[NODES], f_ca[NODES], Ito, r[NODES], s[NODES], IKr, xr1[NODES], xr2[NODES], IKs, xs[NODES], IK1, INaCa, INaK, IpCa, IpK, IbNa, IbCa, Istim, If,xf[NODES],ICaT,dCaT[NODES],fCaT[NODES], Itotal;
  double Ileak, Iup, Irel, CaBuf, CaCSQN, CaCurrent, CaSRCurrent, dCaSR[NODES], bjsr, cjsr, dCai[NODES], bc, cc;

  output  =  fopen( output_file_name, "w" );

  double t;
  int count;

  //Initial values
  int N;
  for (N = 0; N<NODES; N++){
    V[N] =  -70;
    Vn[N]  =  -70;
    Vp[N]  =  -70;
    Cai[N]  =  0.0002;
    CaSR[N]  =  0.2;
    m[N]  =  0.;
    h[N]  =  0.75;
    j[N]  =  0.75;
    xr1[N]  =  0.;
    xr2[N]  =  1.;
    xs[N]  =  0.;
    r[N]  =  0.;
    s[N]  =  1.;
    d[N]  =  0.;
    f[N]  =  1.;
    f_ca[N]  =  1.;
    g[N]  =  1.;
    //Initialization parameters added for hESC_Cm model
    xf[N] = 0.1;
    dCaT[N] = 0.;
    fCaT[N] = 1.;
  }

  t = tmin;
  count  =  0;

  for (t  =  tmin + DT; t <= tmax; t += DT) {
    count++;
    for (N = 0; N<NODES; N++){
      //Reversal potentials
      Eca  =  0.5*(R*T/F)*log(Cao/Cai[N]);
      Ek  =  (R*T/F)*log(Ko/Ki);
      Ena  =  (R*T/F)*log(Nao/Nai);
      Eks  =  (R*T/F)*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
      //INa m gate//////////////////////////////////////////////////////////////
      alpha_m  =  1/(1+exp((-60-V[N])/5));
      beta_m  =  0.1/(1+exp((V[N]+35)/5))+0.1/(1+exp((V[N]-50)/200));
      tau_m  =  alpha_m*beta_m;
      m_inf  =  1/((1+exp((-56.86-V[N])/9.03))*(1+exp((-56.86-V[N])/9.03)));
      //INa, h gate /////////////////////////////////////////////////////////////
      switch(DevelopmentalStage){
        case 'e': h_inf  =  pow((1/(1+exp((V[N]-Vh_h)/k_h))),0.5);
                  break;
        case 'l'://this case proceeds to adult equations
        case 'a': h_inf  =  1/((1+exp((V[N]+71.55)/7.43))*(1+exp((V[N]+71.55)/7.43)));
                  break;
      }
      if ( V[N] >=  -40. ){
        alpha_h  =  0.;
        beta_h  =  ((0.77/(0.13*(1+exp(-(V[N]+10.66)/11.1)))));
      }
      else{
        alpha_h  =  (0.057*exp(-(V[N]+80)/6.8));
        beta_h  =  (2.7*exp(0.079*V[N])+(3.1e5)*exp(0.3485*V[N]));
      }
      tau_h  =  myCoefTauH/(alpha_h+beta_h);
      //INa j gate///////////////////////////////////////////////////////////////
      //j_inf varies depending on the  cell type////////////////////////////////
      switch(DevelopmentalStage){
        case 'e': j_inf  =  pow((1/(1+exp((V[N]-Vh_h)/k_h))),0.5);
                break;
        case 'l'://this case proceeds to adult equations
        case 'a': j_inf  =  1/((1+exp((V[N]+71.55)/7.43))*(1+exp((V[N]+71.55)/7.43)));
                break;
      }
      if( V[N] >=  -40. ){
        alpha_j  =  0.0;
        beta_j  =  ((0.6*exp((0.057)*V[N])/(1+exp(-0.1*(V[N]+32)))));
      }
      else{
        alpha_j  =  (-(25428)*exp(0.2444*V[N])-(0.000006948)*exp(-0.04391*V[N]))*(V[N]+37.78)/(1+exp(0.311*(V[N]+79.23)));
        beta_j  =  ((0.02424*exp(-0.01052*V[N])/(1+exp(-0.1378*(V[N]+40.14)))));
      }
      tau_j  =  myCoefTauJ/(alpha_j+beta_j);
      //ICaL, d gate, and Irel d gate//////////////////////////////////////////
      alpha_d  =  1.4/(1+exp((-35-V[N])/13))+0.25;
      beta_d  =  1.4/(1+exp((V[N]+5)/5));
      gamma_d  =  1/(1+exp((50-V[N])/20));
      tau_d  =  alpha_d*beta_d+gamma_d;
      d_inf  =  1/(1+exp(-(V[N]-Vh_dCa)/kCa));
      //ICaL, f gate///////////////////////////////////////////////////////////
      f_inf  =  1/(1+exp((V[N]+myVhfCaL)/myKfCaL));
      //Equations for this parameter is a bit confusing must email author
      switch(DevelopmentalStage){
        case 'e': tau_f  =  100.;
                 break;
        case 'l'://this case proceeds to adult equations
        case 'a':
                 if(f_inf > f[N]) {tau_f  =  (1125*exp(-((V[N]-myTauFShift)+27)*((V[N]-myTauFShift)+27)/240)+80+165/(1+exp((25-(V[N]-myTauFShift))/10)))*(1+KtaufCa*(Cai[N]-.5e-4));}
                 else tau_f  =  (1125*exp(-((V[N]-myTauFShift)+27)*((V[N]-myTauFShift)+27)/240)+80+165/(1+exp((25-(V[N]-myTauFShift))/10)));
                 break;
      }
      //ICaL, fCa gate///////////////////////////////////////////////////////
      switch(DevelopmentalStage){
        case 'e': f_ca_inf  =  (1/(1+(pow(((Cai[N]-myShiftFCaInf)/0.000325),8)))+0.1/(1+exp(((Cai[N]-myShiftFCaInf)-0.0005)/0.0001))+0.2/(1+exp(((Cai[N]-myShiftFCaInf)-0.00075)/0.0008))+0.23)/1.46;
                  break;
        case 'l': //this case proceeds to adult equations
        case 'a': f_ca_inf  =  (1/(1+(pow((Cai[N]/0.0006),8)))+0.1/(1+exp((Cai[N]-0.0009)/0.0001))+0.3/(1+exp((Cai[N]-0.00075)/0.0008)))/1.3156;
                 break;
      }
      tau_f_ca  =  2.0;//ms
      if ( V[N] > -60.0 ){
        if ( f_ca_inf > f_ca[N] ) f_ca[N]  =  f_ca[N];
        else f_ca[N]  =  f_ca_inf-(f_ca_inf-f_ca[N])*exp(-DT/tau_f_ca);
      }
      else   f_ca[N]  =  f_ca_inf-(f_ca_inf-f_ca[N])*exp(-DT/tau_f_ca);
      //Irel, g gate//////////////////////////////////////////////////////////
      tau_g  =  2.0;//units in ms
      if (Cai[N] <= 0.00035) g_inf  =  (1/(1+pow((Cai[N]/0.00035),6)));
      else g_inf  =  (1/(1+pow((Cai[N]/0.00035),16)));
      if ( V[N] > -60.0 ){
        if ( g_inf > g[N] ) g[N]  =  g[N];
        else g[N]  =  g_inf-(g_inf-g[N])*exp(-DT/tau_g);
      }
      else g[N]  =  g_inf-(g_inf-g[N])*exp(-DT/tau_g);
      //Ito, r gate///////////////////////////////////////////////////////
      r_inf  =  1/(1+exp((-V[N]+20+myShiftItoR)/(6*mySlopeItoR)));
      tau_r  =  myConstTauR*(9.5*exp(-pow((V[N]+40),2)/1800)+0.8);
      //Ito, s gate///////////////////////////////////////////////////////
      s_inf  =  1/(1+exp((V[N]+20+myShiftItoS)/(5*mySlopeItoS)));
      tau_s  =  myConstTauS*(85*exp(-(V[N]+45)*(V[N]+45)/320)+5/(1+exp((V[N]-20)/5))+3);
      //IKs, Xs gate///////////////////////////////////////////////////////
      xs_inf  =  1/(1+exp((-5-V[N])/14));
      alpha_xs  =  1100/(sqrt(1+exp((-10-V[N])/6)));
      beta_xs  =  1/(1+exp((V[N]-60)/20));
      tau_xs  =  alpha_xs*beta_xs;
      //IKr, Xr1 gate///////////////////////////////////////////////////////
      //parameter added for rapid delayed rectifier current
      xr1_inf  =  1/(1+exp((((-R*T/F/Q*log(1/pow(((1+Cao*0.001/Kc)/(1+Cao*0.001/Ka)),4)/L0))-26)-(V[N]-myShift))/7));
      alpha_xr1  =  450/(1+exp((-45-(V[N]-myShift))/(10)));
      beta_xr1  =  6/(1+exp(((V[N]-myShift)-(-30))/11.5));
      tau_xr1  =  alpha_xr1*beta_xr1;
      //IKr, Xr2 gate//////////////////////////////////////////////////////
      xr2_inf  =  1/(1+exp(((V[N]-myShift)-(-88))/24));
      alpha_xr2  =  3/(1+exp((-60-(V[N]-myShift))/20));
      beta_xr2  =  1.12/(1+exp(((V[N]-myShift)-60)/20));
      tau_xr2  =  alpha_xr2*beta_xr2;
      //If, Xf gate////////////////////////////////////////////////////////
      tau_xf  =  1900;//ms
      xf_inf  =  1/(1+exp((V[N]-(-102.4))/(7.6)));
      //ICaT, dCaT gate////////////////////////////////////////////////////
      dCaT_inf  =  1/(1+exp(-(V[N]+26.3)/(6)));
      tau_dCaT  =  1/(1.068*exp((V[N]+26.3)/(30))+1.068*exp(-(V[N]+26.3)/(30)));
      //ICaT, fCaT gate////////////////////////////////////////////////////
      fCaT_inf  =  1/(1+exp((V[N]+61.7)/(5.6)));
      tau_fCaT  =  1/(0.0153*exp(-(V[N]+61.7)/(83.3))+ 0.015*exp((V[N]+61.7)/(15.38)));
      //IK1, alphas and betas//////////////////////////////////////////////
      alpha_K1  =  0.1/(1+exp(0.06*((V[N]-myShiftK1)-Ek-200)));
      beta_K1  =  (3*exp(0.0002*((V[N]-myShiftK1)-Ek+100))+exp(0.1*((V[N]-myShiftK1)-Ek-10)))/(1+exp(-0.5*((V[N]-myShiftK1)-Ek)));
      x_K1_inf  =  alpha_K1/(alpha_K1+beta_K1);

      //Na+ current, INa
      INa  =  Gnamax*m[N]*m[N]*m[N]*h[N]*j[N]*(V[N]-Ena);

      //L-type Ca2+ current, ICaL
      ICaL  =  GCaL*d[N]*f[N]*f_ca[N]*4*V[N]*(F*F)/R/T*(Cai[N]*exp(2*V[N]*F/R/T)-0.341*Cao)/(exp(2*V[N]*F/R/T)-1);

      //Transient outward current, Ito
      Ito  =  GItoepi*r[N]*s[N]*(V[N]-Ek);

      //Rapid delayed rectifier K+ current, IKr
      IKr  =  GKr*xr1[N]*xr2[N]*(V[N]-Ek);

      //Slow delayed rectifier K+ current, IKs, added a scaling factor (1+0.6/(1+pow((3.8e-5/Cai),1.4))
      IKs  =  GKsepi*(1+.6/pow((1+(3.8e-5/Cai[N])),1.4))*xs[N]*xs[N]*(V[N]-Eks);

      //Inward rectifier K+ current, IK1
      IK1 = x_K1_inf*gK1max*(V[N]-Ek);

      //Na+/Ca2+ exchanger current, INaCa, not than n represents lambda
      INaCa  =  knaca*(1./(KmNai*KmNai*KmNai+Nao_naca*Nao_naca*Nao_naca))*(1/(KmCa+Cao))*(1/(1+ksat*exp((n-1)*V[N]*F/(R*T))))*(exp(n*V[N]*F/(R*T))*Nai*Nai*Nai*Cao-exp((n-1)*V[N]*F/(R*T))*Nao_naca*Nao_naca*Nao_naca*Cai[N]*alfa);

      //Na+/K+ pump current, INaK, note that code uses knak variable for pnak
      INaK  =  (1./(1.+0.1245*exp(-0.1*V[N]*F/(R*T))+0.0353*exp(-V[N]*F/(R*T))))*(knak*Ko/(Ko+KmK)*Nai/(Nai+KmNa));

      //Calcium pump current, IpCa
      IpCa  =  GpCa*Cai[N]/(kpca+Cai[N]);

      //Plateau K+ current
      IpK  =  GpK*(V[N]-Ek)*(1./(1.+exp((25-V[N])/5.98)));//Added for adult total

      //Background sodium current
      IbNa  =  GbNa*(V[N]-Ena);

      //Calcium background current, IbCa
      IbCa  =  GbCa*(V[N]-Eca);

      //Currents added for the hESC_CM model
      ////Hyperpolarization activated funny current
      If  =  Gf*xf[N]*(V[N]+17);

      //T-type Ca2+ Current, ICaT
      ICaT  =  dCaT[N]*fCaT[N]*GICaT*(V[N]-Eca);

      /* Current injection */
      if ( N > - 1 && N < 5 ){//First 5 nodes stimulus
        if ( t >=  500 && t <=  500.5 ){
          Istim  =  -52.0;
        }
        else{
          Istim  =  0.0;
        }
      }
      else{
        Istim =  0.0;
      }

      Itotal  =  INa + ICaL + Ito + IKr + IKs + IK1 + INaCa + INaK + IpCa + IbNa + IpK + IbCa + IpK + If + ICaT + Istim;
      
      //Nai = Nai+dt*(-(INa+IbNa+3.*INaK+3.*INaCa)*(1.0/(Vc*F))*CAP);
      //Ki = Ki+dt*(-(Istim+IK1+Ito+IKr+IKs-2.*INaK+IpK)*(1.0/(Vc*F))*CAP);
      //Ki = Ki+dt*(-(IK1+Ito+IKr+IKs-2.*INaK+IpK)*(1.0/(Vc*F))*CAP);

      V[N]  =  V[N] - DT*Itotal;

      //Calcium dynamics
      Ileak  =  Vleak*(CaSR[N]-Cai[N]);
      Iup  =  Vmaxup/(1.+((Kup*Kup)/(Cai[N]*Cai[N])));

      //modified by RaIrel
      Irel  =  (d[N]*g[N]*(crel+arel*(CaSR[N]*CaSR[N])/((brel*brel)+(CaSR[N]*CaSR[N]))))*RaIrel;

      //Analytic solution from the original ten tusscher 2004 model
      CaBuf  =  Bufc*Cai[N]/(Cai[N]+Kbufc);
      CaCSQN  =  Bufsr*CaSR[N]/(CaSR[N]+Kbufsr);
      CaSRCurrent  =  Iup-Irel-Ileak;

      //Added ICaT to CaCurrent
      CaCurrent  =  -(ICaL+IbCa+ICaT+IpCa-2.0*INaCa)*(1.0/(2.0*Vc*F))*capacitance;

      dCaSR[N]  =  DT*(Vc/Vsr)*CaSRCurrent;
      bjsr  =  Bufsr-CaCSQN-dCaSR[N]-CaSR[N]+Kbufsr;
      cjsr  =  Kbufsr*(CaCSQN+dCaSR[N]+CaSR[N]);
      CaSR[N]  =  (sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;

      dCai[N]  =  DT*(CaCurrent-CaSRCurrent);
      bc  =  Bufc-CaBuf-dCai[N]-Cai[N]+Kbufc;
      cc  =  Kbufc*(CaBuf+dCai[N]+Cai[N]);
      Cai[N]  =  (sqrt(bc*bc+4*cc)-bc)/2;

      //Update gates
      m[N]  =  m_inf-(m_inf-m[N])*exp(-DT/(tau_m));
      h[N]  =  h_inf-(h_inf-h[N])*exp(-DT/(tau_h));
      j[N]  =  j_inf-(j_inf-j[N])*exp(-DT/tau_j);
      d[N]  =  d_inf-(d_inf-d[N])*exp(-DT/tau_d);
      f[N]  =  f_inf-(f_inf-f[N])*exp(-DT/tau_f);
      r[N]  =  r_inf-(r_inf-r[N])*exp(-DT/tau_r);
      s[N]  =  s_inf-(s_inf-s[N])*exp(-DT/tau_s);
      xs[N]  =  xs_inf-(xs_inf-xs[N])*exp(-DT/tau_xs);
      xr1[N]  =  xr1_inf-(xr1_inf-xr1[N])*exp(-DT/tau_xr1);
      xr2[N]  =  xr2_inf-(xr2_inf-xr2[N])*exp(-DT/tau_xr2);
      xf[N]  =  xf_inf-(xf_inf-xf[N])*exp(-DT/tau_xf);
      dCaT[N]  =  dCaT_inf-(dCaT_inf-dCaT[N])*exp(-DT/tau_dCaT);
      fCaT[N]  =  fCaT_inf-(fCaT_inf-fCaT[N])*exp(-DT/tau_fCaT);


      /*{
      //fprintf( output, "%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.105f\t%4.10f\t%4.10f\t%4.10f\n", time+(number-1)*basic_cycle_length, Ileak, Iup, Irel, CaBuf, CaSR, Cai, Nai, Ki, V );
      fprintf( output, "%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\n", t, IKr, IKs, IK1, Ito, INa, If, INaK, ICaL, IbCa, INaCa, V,ICaT, Cai, IpCa, Irel, Itotal );
      }
      */
      // printf( "\nV = %f;\nCai = %f;\nCaSR = %f;\nm = %f;\nh = %f;\nj = %f;\nxr1 = %f;\nxr2 = %f;\nxs = %f;\nr = %f;\ns = %f;\nd = %f;\nf = %f;\nf_ca = %f;\ng = %f;\ndCaT = %f;\nfCaT = %f;\nxf = %f;\n", V, Cai, CaSR, m, h, j, xr1, xr2, xs, r, s, d, f, f_ca,g, dCaT, fCaT, xf );
    }//End of Nodes loop
    for (N=0;N<NODES; N++){
      if (N==0)  Vn[N] = DT/(DX*DX)*D*(2.0*V[N+1]-2.0*V[N]) + V[N];
      else if (N==NODES-1)  Vn[N] = DT/(DX*DX)*D*(2.0*V[N-1]-2.0*V[N]) + V[N];
      else  Vn[N] = DT/(DX*DX)*D*(V[N+1]+V[N-1]-2.0*V[N]) + V[N];
    }
    //Calculate the conduction velocity using Trine's method, ask for reasoning, margo deduced
    //that points were chosen to not be two close to the edges of the cable but also far enough
    //apart to get a good estimate of conduction velocity.
    if (Vn[19] >= -40. && Vp[19] < -40.) CVstart = t;
    if (Vn[79] >= -40. && Vp[79] < -40.){
      outputCV=fopen("CV.dat","a");
      fprintf(outputCV,"%.2f\t%.2f\n",t,60.*dx/(t-CVstart)*1000.);//60*dx is the total distance units should be cm/s
      fclose(outputCV);
    }

    for (N=0; N<NODES; N++) {
      V[N] = Vn[N];
      Vp[N] = Vn[N];

    }
    if (!(count % 25))fprintf( output, "%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n", t, Vn[10], Vn[20], Vn[30], Vn[40], Vn[50], Vn[60], Vn[70], Vn[80], Vn[90], Vn[100]);
  }//End of Time loop
  fclose( output );
  return 0;
}

