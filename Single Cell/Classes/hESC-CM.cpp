#include <cmath>
#include <iostream>
#include <cstdlib>
#include <string>
//#include "../include/Environment_LivRudy.h"
#include "Environment_hESCM.h"
#include "hESC-CM.h"
using namespace std;

//Note model was converted from V/s to mV/ms by converting:
// myShiftK1 *1000 to convert from V to mV
// R constant  *1000 to convert the  reversal potentials in mV 
// Maximal conductances from S/F to mS/F or nS/pF
// Time dependent only parameters
// Vleak *(1/1000) to convert 1/s to 1/ms units
// arel and crel *(1/1000) to convert mM/s to mM/ms units
// Vc and Vsr converted from um^3 to picoliters

hESCM::hESCM(char stageChoice){
  DevelopmentalStage = stageChoice;
  switch(DevelopmentalStage){
    case 'e':
      RaINa=0.038;
      myCoefTauH =2.8;
      myCoefTauJ = 1;
      RaICaL=0.25;
      Vh_dCa=12.5;
      kCa=12.5;
      ampICaLDinf = 1;
      KtaufCa=1433;   // [1/mM]  Altamirano & Bers 2007
      myVhfCaL = 20;
      myKfCaL = 7;
      ampICaLFinf = 1;
      myKtauF = 1;
      myTauFShift = 0;
      myShiftFCaInf = -0.11e-3;  //[mM] 
      Vth_ICaL = -0.060;
      RaICaT =0.25;
      RaINaCa=1.750e1;
      alfa=0.8;
      RaINaK=0.7;
      RaIK1=0.05*2.67/3;
      myShiftK1 = -15;//Converted to mV
      RaIKr=3;
      mySlope = 1;
      myShift = 0;
      RaIKs=0.1;
      RaIto=0.1673*0.4903*0.8;
      myShiftItoR = -25;
      mySlopeItoR = 0.3;
      myShiftItoS = 0;
      mySlopeItoS = 1;
      myConstTauR= 1;
      myConstTauS= 1;
      RaIup=0.4/3;
      RaIrel=0.2/18;
      RaIleak=0.1/18;
      RaIf=0.5389;
      myIfconst = 1;
      RaCm=0.22162;
      RaICap=1;
      RaIKp=0;
      RaIback=0.2;
/*
      V = -70.709774;
      Cai = 0.000026;
      CaSR = 0.123913;
      m = 0.031486;
      h = 0.644831;
      j = 0.637497;
      xr1 = 0.002078;
      xr2 = 0.327344;
      xs = 0.009054;
      r = 0.000000;
      s = 0.999961;
      d = 0.001283;
      f = 0.999295;
      f_ca = 1.002155;
      g = 1.000000;
      dCaT = 0.000610;
      fCaT = 0.836858;
      xf = 0.014176;
*/
      break;
    case 'l':
      RaINa=1;    // Itoh
      myCoefTauH = 2.8;
      myCoefTauJ = 1;
      RaICaL=0.422;
      Vh_dCa=16;
      kCa=12.8;
      ampICaLDinf = 1;
      KtaufCa=1433;   // [1/mM]  Altamirano & Bers 2007
      myVhfCaL = 20;
      myKfCaL = 7;
      ampICaLFinf = 1;
      myKtauF = 1;
      myTauFShift = 0;
      myShiftFCaInf = -0.12e-3;
      Vth_ICaL = 0;
      RaICaT =0.05;
      RaINaCa=1.824e+01;
      alfa=0.38;
      RaINaK=0.83;
      RaIK1=0.4*0.05*2.67*4;
      myShiftK1 = -15;//Converted to mV 
      RaIKr=1.4;
      mySlope = 1;
      myShift = 0;
      RaIKs=0.1;
      RaIto=0.3754*0.4903*0.9;
      myShiftItoR = -25;
      mySlopeItoR = 0.3;
      myShiftItoS = 0;
      mySlopeItoS = 1;
      myConstTauR= 1;
      myConstTauS= 1;
      RaIup=0.33;
      RaIrel=0.4;
      RaIleak=0.3*1;
      RaIf=0.23;
      RaCm=0.17838;
      RaICap=1;
      RaIKp=0;
      RaIback=1;
/*
      V = -70.350507;
      Cai = 0.000072;
      CaSR = 0.095304;
      m = 0.033609;
      h = 0.226380;
      j = 0.194417;
      xr1 = 0.034771;
      xr2 = 0.324050;
      xs = 0.009288;
      r = 0.000000;
      s = 0.999958;
      d = 0.001174;
      f = 0.998545;
      f_ca = 0.995698;
      g = 0.999920;
      dCaT = 0.000647;
      fCaT = 0.828093;
      xf = 0.010919;
*/
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
      RaCm = 1;
      RaINaK = 1;
      RaICap = 1;
      RaIKp = 1;
      RaIback = 1;
/*
      //Initial values
      V= -70;
      Cai = 0.0002;
      CaSR = 0.2;
      m = 0.;
      h = 0.75;
      j = 0.75;
      xr1 = 0.;
      xr2 = 1.;
      xs = 0.;
      r = 0.;
      s = 1.;
      d = 0.;
      f = 1.;
      f_ca = 1.;
      g = 1.;
      //Initialization parameters added for hESC_Cm model
      xf=0.1;
      dCaT=0.;
      fCaT=1.;
*/
      break;
    default: cout <<"No developmental stage selected."<<endl;
             break;

  }
  //Flag to choose between epi, endo and M cell types
  epi=1;
  endo=0;
  Mcell=0;

  //Constants
  R=8314.472;   // [J/millimoles/K] Gas constant
  F=96485.3415; // [C/mol]	  Faraday constant
  T=310.0;      // [K]       Temperature

  //Buffering
  Bufc=0.25;   // [mM] total cytoplasmic buffer concentration
  Kbufc=0.001;  // [mM] Cai half saturation constant
  Bufsr=10;     // [mM] total sarcoplasmic buffer concentration
  Kbufsr=0.3;   // [mM] CaSR half saturation constant

  //Extracellular Ionic concentrations
  //Commented because parameters are defined in Environment.h header
  //Ko=4;         // [mM]
  //Cao=2.7;      // [mM]
  //Nao=150.5;    // [mM]
  Nao_naca = Nao;

  //Intracellular Ionic concentrations
  //Pre-dialysis
  Ki=140;      // [mM]  & 140 in TT04
  Nai=7;    // [mM]

  // Intracellular Volumes
  Vc=16.404*RaCm;//picoliters
  Vsr=1.094*RaCm;//picoliters
  capacitance=0.185*1000*RaCm;//pF

  //Ionic Currents
  //Fast Na+ Current
  Vh_h=-73;
  k_h=5.6;
  Gnamax=14.838*RaINa; // [nS/pF] maximal INa conductance

  //If Current
  Gf=0.090926*RaIf; // [nS/pF] maximal If conductance

  //L-type Ca2+ Current, not sure why this stays the same
  GCaL=0.000175*RaICaL;  // [m^3/F/s] maximal ICaL conductance

  //T-type Ca2+ Current
  GICaT = 0.1832*RaICaT; //[nS/pF]

  //Transient Outward Current
  GItoepi=0.294*RaIto;   // [nS/pF] maximal ITo conductance
  GItoendo=0.073;   // [nS/pF] maximal ITo conductance
  GItoMcell=0.294; // [nS/pF] maximal ITo conductance

  //IKs
  GKsepi   =0.157*RaIKs; //[nS/pF] maximal IKs conductance
  GKsendo  =0.157; // [nS/pF] maximal IKs conductance
  GKsMcell =0.040; // [nS/pF] maximal IKs conductance
  pKNa=0.03;   // [dimensionless]

  //IKr
  GKr=0.096*sqrt(Ko/5.4)*RaIKr; //[nS/pF] maximal IKr conductance
  Q=2.3;
  L0=0.025;
  Kc=0.58e-3;
  Ka=2.6e-3;

  //Inward Rectifier K+ Current
  gK1max=5.405*sqrt(Ko/5.4)*RaIK1; // [nS/pF] maximal IK1 conductance

  //Na+/Ca2+ Exchanger Current
  knaca=1000*RaINaCa;  // [pA/pF] maximal INaCa
  KmNai=87.5;  // [mM] Nai half saturation constant
  KmCa=1.38;   // [mM] Cai half saturation constant
  ksat=0.1;    // [dimensionless]  saturation factor for INaCa
  n=0.35;      // [dimensionless]  voltage dependence parameter

  //Na+/K+ Pump Current
  knak=1.362*RaINaK;  // [pA/pF] maximal INaK
  KmK=1;       // [mM] Ko half saturation constant
  KmNa=40;     // [mM] Nai half saturation constant

  //IpCa
  GpCa=0.825*RaICap;    // [pA/pF] maximal IpCa
  kpca=0.0005;   // [mM]  Cai half saturation constant

  //IpK
  GpK=0.0146*RaIKp;    // [nS/pF] maximal IpK conductance

  //Background Currents
  GbNa=0*0.29*RaIback;   // [nS/pF] maximal IbNa conductance
  GbCa=0.000592*RaIback;  // [nS/pF] maximal IbCa conductance

  //Calcium Dynamics
  //Ileak
  Vleak=0.00008*RaIleak; // [1/ms] maximal Ileak =0.00008/s

  //Irel
  arel=0.016464;   // [mM/ms] maximal CaSR-dependent Irel
  brel=0.25;     // [mM] CaSR half saturation constant
  crel=0.008232;    // [mM/ms] maximal CaSR-independent Irel

  //Iup
  Vmaxup=0.000425*RaIup; // [mM/ms] maximal Iup
  Kup=0.00025;    // [mM] half saturation constant
  
  //Initial values
  for (int node = 0; node<CELLS_fib; node++){
        yh[node][0]= -0.069160;//V
        yh[node][1]= 0.000027;//Cai
        yh[node][2]= 0.086301;//CaSR
        yh[node][3]= 0.041569;//m
        yh[node][4]= 0.600505;//h
        yh[node][5]= 0.618972;//j
        yh[node][6]= 0.001208;//xr1
        yh[node][7]= 0.313321;//xr2
        yh[node][8]= 0.010087;//xs
        yh[node][9]= 0.000000;//r
        yh[node][10]= 0.999948;//s
        yh[node][11]= 0.001452;//d
        yh[node][12]= 0.999173;//f
        yh[node][13]= 1.002065;//f_ca
        yh[node][14]= 1.000000;//g
        yh[node][15]= 0.000789;//dCaT
        yh[node][16]= 0.797842;//fCat
        yh[node][17]= 0.014191;//xf
  }

  Igap_hes = 0;
  //Initialize parameter calculation
  for(counter=0;counter<numOfAPDs;counter++){
    vmin[counter] = 100000.0;
    vmax[counter] = -10000.0;
    dvdtmax[counter]       = -10000.0;
    ddr[counter]           = -10000.0;
    top[counter]           =  10000.0;
    top_slope[counter]     = -10000.0;
    apd30[counter]         = -10000.0;
    apd50[counter]         = -10000.0;
    apd70[counter]         = -10000.0;
    apd90[counter]         = -10000.0;
    bcl[counter]  = -10000.0;
  }
  //printf("Vmin(mV)\tVmax(mV)\tdVdtmax(mV/s)\tAPD30(ms)\tAPD50(ms)\tAPD70(ms)\tAPD90(ms)\tFreq.(bpm)\tDDR(mV/s)\tTOP(mV)\n");
  counter =0;
  start_output = 0;
  dvdt = 1000.0;
  dvdtold = 500.0;
  start_output = 0;
}
hESCM::~hESCM(){
}
string hESCM::getName(){
  switch(DevelopmentalStage){
    case 'e': return "hESCM_Early";
              break;
    case 'l': return "hESCM_Late";
              break;
    case 'a': return "Adult";
              break;
    default:
              cout<<"Developmental Stage not Specified!"<<endl;
              break;
  }
}
void hESCM::solve(double Istim, double t){
  //Scale Igap by the fibroblast capacitance
 // Igap_hes = Igap * (1/(capacitance));

  //Start of the hESCM dimensions
          ///////////////////////////////////////////////////////////////////////////////////////////
          for (int x = 0; x<CELLS_fib; x++){
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

            //Reversal potentials
            Eca = 0.5*(R*T/F)*log(Cao/yh[x][1]);
            Ek = (R*T/F)*log(Ko/Ki);
            Ena = (R*T/F)*log(Nao/Nai);
            Eks = (R*T/F)*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
            //INa m gate//////////////////////////////////////////////////////////////
            alpha_m = 1/(1+exp((-60-yh[x][0])/5));
            beta_m = 0.1/(1+exp((yh[x][0]+35)/5))+0.1/(1+exp((yh[x][0]-50)/200));
            tau_m = alpha_m*beta_m;
            m_inf = 1/((1+exp((-56.86-yh[x][0])/9.03))*(1+exp((-56.86-yh[x][0])/9.03)));
            //m = m_inf-(m_inf-m)*exp(-DT/(tau_m));
            yhdot[x][3] = (m_inf - yh[x][3])/tau_m;
            //INa, h gate /////////////////////////////////////////////////////////////
            switch(DevelopmentalStage){
              case 'e': h_inf = pow((1/(1+exp((yh[x][0]-Vh_h)/k_h))),0.5);//yh[x][0]h_h = -73, k_h = 5.6
                        break;
              case 'l': //this case proceeds to the adult phenotype
              case 'a': h_inf = 1/((1+exp((yh[x][0]+71.55)/7.43))*(1+exp((yh[x][0]+71.55)/7.43)));
                        break;
            }
            if( yh[x][0] < -40. ){
              alpha_h = (0.057*exp(-(yh[x][0]+80)/6.8));
              beta_h = (2.7*exp(0.079*yh[x][0])+(3.1e5)*exp(0.3485*yh[x][0]));
            }
            else{
              alpha_h = 0.;
              beta_h = ((0.77/(0.13*(1+exp(-(yh[x][0]+10.66)/11.1)))));
            }
            tau_h = myCoefTauH/(alpha_h+beta_h);//myCoefTauH = 2.8
            //h = h_inf-(h_inf-h)*exp(-DT/(tau_h));
            yhdot[x][4] = (h_inf - yh[x][4])/tau_h;
            //INa j gate///////////////////////////////////////////////////////////////
            //j_inf varies depending on the  cell type////////////////////////////////
            switch(DevelopmentalStage){
              case 'e': j_inf = pow((1/(1+exp((yh[x][0]-Vh_h)/k_h))),0.5);
                        break;
              case 'l'://this case proceeds to adult phenotype
              case 'a': j_inf = 1/((1+exp((yh[x][0]+71.55)/7.43))*(1+exp((yh[x][0]+71.55)/7.43)));
                        break;
            }
            if( yh[x][0] < -40. ){
              alpha_j = (-(25428)*exp(0.2444*yh[x][0])-(0.000006948)*exp(-0.04391*yh[x][0]))*(yh[x][0]+37.78)/(1+exp(0.311*(yh[x][0]+79.23)));
              beta_j = ((0.02424*exp(-0.01052*yh[x][0])/(1+exp(-0.1378*(yh[x][0]+40.14)))));
            }
            else{
              alpha_j = 0.0;
              beta_j = ((0.6*exp((0.057)*yh[x][0])/(1+exp(-0.1*(yh[x][0]+32)))));
            }
            tau_j = myCoefTauJ/(alpha_j+beta_j);
            //j = j_inf-(j_inf-j)*exp(-DT/tau_j);
            yhdot[x][5] = (j_inf - yh[x][5])/tau_j;

            //ICaL, d gate, and Irel d gate//////////////////////////////////////////
            alpha_d = 1.4/(1+exp((-35-yh[x][0])/13))+0.25;
            beta_d = 1.4/(1+exp((yh[x][0]+5)/5));
            gamma_d = 1/(1+exp((50-yh[x][0])/20));
            tau_d = alpha_d*beta_d+gamma_d;

            d_inf = 1/(1+exp(-(yh[x][0]-Vh_dCa)/kCa));

            //d = d_inf-(d_inf-d)*exp(-DT/tau_d);
            yhdot[x][11]=(d_inf - yh[x][11])/d_inf;

            //ICaL, f gate///////////////////////////////////////////////////////////
            f_inf = 1/(1+exp((yh[x][0]+myVhfCaL)/myKfCaL));
            switch(DevelopmentalStage){
              case 'e': tau_f = 100.;
                        break;
              case 'l'://this case proceeds to adult phenotype
              case 'a':
                        if(f_inf > yh[x][12]){ tau_f = (1125*exp(-((yh[x][0]-myTauFShift)+27)*((yh[x][0]-myTauFShift)+27)/240)+80+165/(1+exp((25-(yh[x][0]-myTauFShift))/10)))*(1+KtaufCa*(yh[x][1]-.5e-4));}
                        else {tau_f = (1125*exp(-((yh[x][0]-myTauFShift)+27)*((yh[x][0]-myTauFShift)+27)/240)+80+165/(1+exp((25-(yh[x][0]-myTauFShift))/10)));}
                        break;
            }
            //f = f_inf-(f_inf-f)*exp(-DT/tau_f);
            yhdot[x][12]=(f_inf - yh[x][12])/tau_f;
            //ICaL, fCa gate/////////////////////////////////////////////////////// 
            switch(DevelopmentalStage){
              case 'e':     f_ca_inf = (1/(1+(pow(((yh[x][1]-myShiftFCaInf)/0.000325),8)))+0.1/(1+exp(((yh[x][1]-myShiftFCaInf)-0.0005)/0.0001))+0.2/(1+exp(((yh[x][1]-myShiftFCaInf)-0.00075)/0.0008))+0.23)/1.46;
                            break;
              case 'l'://this case proceeds to adult phenotype    
              case 'a': f_ca_inf = (1/(1+(pow((yh[x][1]/0.0006),8)))+0.1/(1+exp((yh[x][1]-0.0009)/0.0001))+0.3/(1+exp((yh[x][1]-0.00075)/0.0008)))/1.3156;
                        break;
            }
            tau_f_ca = 2.0;//ms
            if ( yh[x][0] > -60.0 ){
              if ( f_ca_inf > yh[x][13] ) { yhdot[x][13] = yh[x][13];}
              else { 
                //f_ca = f_ca_inf-(f_ca_inf-f_ca)*exp(-DT/tau_f_ca);
                yhdot[x][13]=(f_ca_inf - yh[x][13])/tau_f_ca;
              }
            }
            else { 
              //f_ca = f_ca_inf-(f_ca_inf-f_ca)*exp(-DT/tau_f_ca);
              yhdot[x][13]=(f_ca_inf - yh[x][13])/tau_f_ca;
            }
            //Irel, g gate////////////////////////////////////////////////////////// 
            tau_g = 2.0;//units in ms
            if (yh[x][1]<=0.00035) {g_inf = (1/(1+pow((yh[x][1]/0.00035),6)));}
            else { g_inf = (1/(1+pow((yh[x][1]/0.00035),16)));}
            if ( yh[x][0] > -60.0 ){
              if (g_inf > yh[x][14]) {yhdot[x][14]= yh[x][14];}
              else { 
                //g = g_inf-(g_inf-g)*exp(-DT/tau_g);
                yhdot[x][14]=(g_inf - yh[x][14])/tau_g;
              }
            }
            else { 
              //g = g_inf-(g_inf-g)*exp(-DT/tau_g);
              yhdot[x][14]=(g_inf - yh[x][14])/tau_g;
            }
            //IKs, Xs gate/////////////////////////////////////////////////////// 
            xs_inf = 1/(1+exp((-5-yh[x][0])/14));
            alpha_xs = 1100/(sqrt(1+exp((-10-yh[x][0])/6)));
            beta_xs = 1/(1+exp((yh[x][0]-60)/20));
            tau_xs = alpha_xs*beta_xs;
            //xs = xs_inf-(xs_inf-xs)*exp(-DT/tau_xs);
            yhdot[x][8] = (xs_inf - yh[x][8])/tau_xs;
            //Ito, r gate/////////////////////////////////////////////////////// 
            r_inf = 1/(1+exp((-yh[x][0]+20+myShiftItoR)/(6*mySlopeItoR)));
            tau_r = myConstTauR*(9.5*exp(-pow((yh[x][0]+40),2)/1800)+0.8);
            //r = r_inf-(r_inf-r)*exp(-DT/tau_r);
            yhdot[x][9] = (r_inf - yh[x][9])/tau_r;
            //Ito, s gate/////////////////////////////////////////////////////// 
            s_inf = 1/(1+exp((yh[x][0]+20+myShiftItoS)/(5*mySlopeItoS)));
            tau_s = myConstTauS*(85*exp(-(yh[x][0]+45)*(yh[x][0]+45)/320)+5/(1+exp((yh[x][0]-20)/5))+3);
            //s = s_inf-(s_inf-s)*exp(-DT/tau_s);
            yhdot[x][10] = (s_inf - yh[x][10])/tau_s;
            //If, Xf gate////////////////////////////////////////////////////////
            tau_xf = 1900;//ms
            xf_inf = 1/(1+exp((yh[x][0]-(-102.4))/(7.6)));
            //xf = xf_inf-(xf_inf-xf)*exp(-DT/tau_xf);
            yhdot[x][17] = (xf_inf - yh[x][17])/tau_xf;
            //ICaT, dCaT gate//////////////////////////////////////////////////// 
            dCaT_inf = 1/(1+exp(-(yh[x][0]+26.3)/(6)));
            tau_dCaT = 1/(1.068*exp((yh[x][0]+26.3)/(30))+1.068*exp(-(yh[x][0]+26.3)/(30)));
            //dCaT = dCaT_inf-(dCaT_inf-dCaT)*exp(-DT/tau_dCaT);
            yhdot[x][15] = (dCaT_inf - yh[x][15])/tau_dCaT;
            //ICaT, fCaT gate//////////////////////////////////////////////////// 
            fCaT_inf = 1/(1+exp((yh[x][0]+61.7)/(5.6)));
            tau_fCaT = 1/(0.0153*exp(-(yh[x][0]+61.7)/(83.3))+ 0.015*exp((yh[x][0]+61.7)/(15.38)));
            //fCaT = fCaT_inf-(fCaT_inf-fCaT)*exp(-DT/tau_fCaT);
            yhdot[x][16] = (fCaT_inf - yh[x][16])/tau_fCaT;
            //IK1, alphas and betas////////////////////////////////////////////// 
            alpha_K1 = 0.1/(1+exp(0.06*((yh[x][0]-myShiftK1)-Ek-200)));
            beta_K1 = (3*exp(0.0002*((yh[x][0]-myShiftK1)-Ek+100))+exp(0.1*((yh[x][0]-myShiftK1)-Ek-10)))/(1+exp(-0.5*((yh[x][0]-myShiftK1)-Ek)));
            x_K1_inf = alpha_K1/(alpha_K1+beta_K1);


            //Na+ current, INa 
            INa = Gnamax*yh[x][3]*yh[x][3]*yh[x][3]*yh[x][4]*yh[x][5]*(yh[x][0]-Ena);

            //L-type Ca2+ current, ICaL 
            ICaL = GCaL*yh[x][11]*yh[x][12]*yh[x][13]*4*yh[x][0]*(F*F)/R/T*(yh[x][1]*exp(2*yh[x][0]*F/R/T)-0.341*Cao)/(exp(2*yh[x][0]*F/R/T)-1);

            //Transient outward current, Ito 
            Ito = GItoepi*yh[x][9]*yh[x][10]*(yh[x][0]-Ek);

            //Rapid delayed rectifier K+ current, IKr 
            IKr = GKr*yh[x][6]*yh[x][7]*(yh[x][0]-Ek);

            //Slow delayed rectifier K+ current, IKs, added a scaling factor (1+0.6/(1+pow((3.8e-5/yh[x][1]),1.4))
            IKs = GKsepi*(1+.6/pow((1+(3.8e-5/yh[x][1])),1.4))*yh[x][8]*yh[x][8]*(yh[x][0]-Eks);

            //Inward rectifier K+ current, IK1
            IK1=x_K1_inf*gK1max*(yh[x][0]-Ek);

            //Na+/Ca2+ exchanger current, INaCa, not than n represents lambda 
            INaCa = knaca*(1/(KmNai_hesc*KmNai_hesc*KmNai_hesc+Nao_naca*Nao_naca*Nao_naca))*(1/(KmCa+Cao))*(1/(1+ksat_hesc*exp((n-1)*yh[x][0]*F/(R*T))))*(exp(n*yh[x][0]*F/(R*T))*Nai*Nai*Nai*Cao-exp((n-1)*yh[x][0]*F/(R*T))*Nao_naca*Nao_naca*Nao_naca*yh[x][1]*alfa);

            //Na+/K+ pump current, INaK, note that code uses knak variable for pnak
            INaK = (1/(1+0.1245*exp(-0.1*yh[x][0]*F/(R*T))+0.0353*exp(-yh[x][0]*F/(R*T))))*(knak*Ko/(Ko+KmK)*Nai/(Nai+KmNa));

            //Calcium pump current, IpCa 
            IpCa = GpCa*yh[x][1]/(kpca+yh[x][1]);

            //Plateau K+ current 
            IpK = GpK*(yh[x][0]-Ek)*(1/(1+exp((25-yh[x][0])/5.98)));

            //Background sodium current 
            IbNa = GbNa*(yh[x][0]-Ena);

            //Calcium background current, IbCa 
            IbCa = GbCa*(yh[x][0]-Eca);

            //Currents added for the hESC_CM model
            ////Hyperpolarization activated funny current, If
            If = Gf*yh[x][17]*(yh[x][0]+17);

            //T-type Ca2+ Current, ICaT, not originally in the tentusscher model 
            ICaT = yh[x][15]*yh[x][16]*GICaT*(yh[x][0]-Eca);


            //dvdtold = dvdt;

            Itotal = INa + ICaL + Ito + IKr + IKs + IK1 + INaCa + INaK + IpCa + IbNa + IpK + IbCa + If + ICaT + I_stim_fib;

            Ileak = Vleak*(yh[x][2]-yh[x][1]);
            Iup = Vmaxup/(1.+((Kup*Kup)/(yh[x][1]*yh[x][1])));

            //modified by RaIrel
            Irel = (yh[x][11]*yh[x][14]*(crel+arel*(yh[x][2]*yh[x][2])/((brel*brel)+(yh[x][2]*yh[x][2]))))*RaIrel;

            CaBuf = Bufc*yh[x][1]/(yh[x][1]+Kbufc);
            CaCSQN = Bufsr*yh[x][2]/(yh[x][2]+Kbufsr);
            CaSRCurrent = Iup-Irel-Ileak;

            //Added ICaT to CaCurrent
            CaCurrent = -(ICaL+IbCa+ICaT+IpCa-2.0*INaCa)*(1.0/(2.0*Vc*F))*capacitance;

            //dCaSR = DT*(Vc/Vsr)*CaSRCurrent;
            yhdot[x][2] = (Vc/Vsr_hesc)*CaSRCurrent;
            bjsr = Bufsr-CaCSQN-dCaSR-yh[x][2]+Kbufsr;
            cjsr = Kbufsr*(CaCSQN+dCaSR+yh[x][2]);
            yh[x][2] = (sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;

            //dCai = DT*(CaCurrent-CaSRCurrent);
            yhdot[x][1] = (CaCurrent-CaSRCurrent);
            bc = Bufc-CaBuf-dCai-yh[x][1]+Kbufc;
            cc = Kbufc*(CaBuf+dCai+yh[x][1]);
            yh[x][1] = (sqrt(bc*bc+4*cc)-bc)/2;

            // Current injection
            if ( (( count <= stim_dur_int)&& (stim_fib[x] == 1))&& (num> stim_equil)) {
              I_stim_fib = -stim_mag;
            }
            else {
              I_stim_fib = 0.0;

            }


            yhdot[x][0] = -Itotal;

            //Update membrane potential and gates
            for (int cc=0; cc<18; cc++) yh[x][cc] = yh[x][cc] + DT * yhdot[x][cc];

            if(count%file_filter==0) fprintf(output2,"%.12f\t",yh[x][0]);
          }
/*
//Method for parameters calculations for spontaneously beating cells taken from Kharche et al 2010
  if(dvdt>=0.0 && dvdtold<0.0){
    vmin[counter] = V;
    tvmin[counter] = t;
    start_output = 1;
  }
  if(dvdt>dvdtmax[counter]&&start_output>0){
    dvdtmax[counter] = dvdt;
    apd_start[counter] = t;
  }
  if(dvdtold>0.0&&dvdt<=0.0){
    vmax[counter] = V;
    top_slope[counter] = (vmax[counter]-vmin[counter])/(t - tvmin[counter]);
  }
  if((counter>0)&&(dvdtold<=top_slope[counter-1])&&(dvdt>top_slope[counter-1])){
    top[counter] = V;
    ddr[counter] = (V - vmin[counter])/(t - tvmin[counter]);
  }
  if(vnew<=0.3*vmin[counter] && V>0.3*vmin[counter]){
    if(apd_start[counter]>0.0)
      apd30[counter] = t - apd_start[counter];
  }
  if(vnew<=0.5*vmin[counter] && V>0.5*vmin[counter]){
    if(apd_start[counter]>0.0)
      apd50[counter] = t - apd_start[counter];
  }
  if(vnew<=0.7*vmin[counter] && V>0.7*vmin[counter]){
    if(apd_start[counter]>0.0)
      apd70[counter] = t - apd_start[counter];
  }
  if(vnew<=0.9*vmin[counter]&&V>0.9*vmin[counter]){
    if(apd_start[counter]>0.0){
      apd_end[counter] = t;
      apd90[counter] = t - apd_start[counter];
      bcl[counter] = apd_start[counter]-apd_start[counter-1];
      printf("%10.5f\t", vmin[counter]);
      printf("%10.5f\t", vmax[counter]);
      printf("%10.5f\t", dvdtmax[counter]*1000);
      printf("%10.5f\t", apd30[counter]);
      printf("%10.5f\t", apd50[counter]);
      printf("%10.5f\t", apd70[counter]);
      printf("%10.5f\t", apd90[counter]);
      printf("%10.5f\t", (60000/bcl[counter]));//convert cycle length to beats per minute
      printf("%10.5f\t", ddr[counter]*1000);
      printf("%10.5f\t", top[counter]);//take off potential
      printf("\n");
      counter++;
    }
  }
  */
 // V = vnew;
}
double hESCM::getV(){
  return V;//Converts voltage to mV
}
void hESCM::print_vars(bool FLAG, FILE *output,double t, int stimulus_number)
{
  if (FLAG == true) {
    printf( "\nV = %f;\nCai = %f;\nCaSR = %f;\nm = %f;\nh = %f;\nj = %f;\nxr1 = %f;\nxr2 = %f;\nxs = %f;\nr = %f;\ns = %f;\nd = %f;\nf = %f;\nf_ca = %f;\ng = %f;\ndCaT = %f;\nfCaT = %f;\nxf = %f;\n", V, Cai, CaSR, m, h, j, xr1, xr2, xs, r, s, d, f, f_ca,g, dCaT, fCaT, xf );
  }
  else {
    if(print_header==0){
      // fprintf(output,"Time(ms)\tVoltage(mV)\tI_K_1(pA/pF)\tI_K_v(pA/pF)\tI_N_a_K(pA/pF)\tI_b_N_a(pA/pF)\tI_g_a_p(pA/pF)\n");
      print_header++;
    }
    fprintf( output, "%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\n", t, IKr, IKs, IK1, Ito, INa, If, INaK, ICaL, IbCa, INaCa, V,ICaT, Cai, IpCa, Irel, Itotal );
  }
}

