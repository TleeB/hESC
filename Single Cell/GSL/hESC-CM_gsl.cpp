#include <cstdio>
#include <cmath>
#include <iostream>
//gsl libraries for ode solver
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#define NAME 20 //characters in file name
#define EPI
//#include "Parameters.h"//contains global variables of model parameters
//#define EARLY
#define LATE

using namespace std;
//double DT = 0.0000001;
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
 double xr1o;
 double xr2o;
 double x0;
 double dx;
};
//ODE variables
const int dimension = 17; //number of ODE derivatives
double y[dimension];
double dy[dimension];
double t, t_next;		// current and next independent variable
double tmin, tmax, delta_t;	// range of t and step size for output
double dh;//note that gsl example uses h, but h is one of the state variables!
double dfdttemp;

// ODE function prototype
static int nstates(double , const double [], double [], void *);


int main(){
  FILE *output;
  char output_file_name[ NAME ];
  output = fopen( "hescm.dat", "w" );
  int file_output_counter=0;
////////////////////////////////////////////////////////////////
    ode_params multipliers;
    ode_params *odeParams;
    double eps_abs = 1.e-6;	// absolute error requested
    double eps_rel = 1.e-6;	// relative error requested
    
// define the type of routine for making steps:
    const gsl_odeiv2_step_type *type_ptr = gsl_odeiv2_step_rk4;
    //rk2
    //rkf45
    //gsl_odeiv2_step_msadams
    //   = gsl_odeiv2_step_rk4;
    //   = gsl_odeiv2_step_rkck;
    //   = gsl_odeiv2_step_rk8pd;
    //   = gsl_odeiv2_step_rk4imp;
    
// allocate/initialize the stepper, the control function, and the evolution function.
    gsl_odeiv2_step *step_ptr = gsl_odeiv2_step_alloc (type_ptr, dimension);
    gsl_odeiv2_control *control_ptr = gsl_odeiv2_control_y_new (eps_abs, eps_rel);
    gsl_odeiv2_evolve *evolve_ptr = gsl_odeiv2_evolve_alloc (dimension);

    gsl_odeiv2_system my_system;	// structure with the diffeq function, etc.
    
    
// load values into the my_system structure
    my_system.function = &nstates;	// the right-hand-side functions dy[i]/dt
    my_system.jacobian = NULL;
    my_system.dimension = dimension;	// number of diffeq's
    my_system.params = &multipliers;	// parameters to pass to function and jacobian
    
//need to use driver function for some of the algorithms
    gsl_odeiv2_driver * driver_ptr = gsl_odeiv2_driver_alloc_y_new (&my_system,type_ptr,1e-7,eps_rel,eps_abs);
    gsl_odeiv2_step_set_driver (step_ptr, driver_ptr);
    gsl_odeiv2_control_set_driver (control_ptr, driver_ptr);
    gsl_odeiv2_evolve_set_driver (evolve_ptr, driver_ptr);

//Define a pointer for the defines, segmentation fault without
odeParams = &multipliers;

#include "odedefines.h"

#ifdef EARLY
 celltype = 'e';
 RaINa=0.038;
 myCoefTauH =2.8;
 myCoefTauJ = 1;
 RaICaL=0.25;
 Vh_dCa=12.5;
 kCa=12.5;
 ampICaLDinf = 1;
 KtaufCa=1433;   // [1/mM]  Altamirano & Bers 2007
 myVhfCaL = 20;
 myKfCaL = 7;//7
 ampICaLFinf = 1;//to accelerate the tauF and to inactivate more rapidly Ical
 myKtauF = 1;
 myTauFShift = 0;
 myShiftFCaInf = -0.11e-3;  //[mM] to shift a litte bit the sigmoidal of fCa
 Vth_ICaL = -0.060;
 RaICaT =0.25;
 RaINaCa=1.750e1;
 alfa=0.8;
 RaINaK=0.7;
 RaIK1=0.05*2.67/3; 
 myShiftK1 = -0.015;
 RaIKr=3;
 poffi = 0;
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
 Rax0=1.1061;
 Radx=0.6537;
 myIfconst = 1;
 RaCm=0.22162;
 RaICap=1;
 RaIKp=0;
 RaIback=0.2;
#endif

#ifdef LATE
char celltype='l';
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
 myShiftFCaInf = -0.12e-3; //[mM] to move a bit left the sigmoid  fCa
 Vth_ICaL = 0;
 RaICaT =0.05;
 RaINaCa=1.824e+01;
 alfa=0.38;
 RaINaK=0.83; 
 RaIK1=0.4*0.05*2.67*4;
 myShiftK1 = -0.015;
 RaIKr=1.4;
 poffi = 0;
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
 Rax0=1.1061;
 Radx=0.6537;
 RaCm=0.17838;
 RaICap=1;
 RaIKp=0;
 RaIback=1;
#endif

//// Constants
 R=8.314472;   // [J/mol/K] Gas constant
 F=96485.3415; // [C/mol]	  Faraday constant
 T=310.0;      // [K]       Temperature 37 degrees celsius
 RTONF = (R*T)/F;
//// Buffering
 Bufc=0.25;   // [mM] total cytoplasmic buffer concentration
 Kbufc=0.001;  // [mM] Cai half saturation constant
 Bufsr=10;     // [mM] total sarcoplasmic buffer concentration
 Kbufsr=0.3;   // [mM] CaSR half saturation constant
//// Extracellular Ionic concentrations
 Ko=4;         // [mM]
 Cao=2.7;      // [mM]
 Nao=150.5;    // [mM]
 Nao_naca = Nao;
//// Intracellular Ionic concentrations
 Ki=140;      // [mM]  & 140 in TT04
 Nai=7;    // [mM]
 Vc=0.016404*RaCm;
 Vsr=0.001094*RaCm;
 capacitance=0.185*1000*RaCm;
//// Resting Potential
 svolt=-70;  // [mV]
//// Ionic Currents

//// Fast Na+ Current 
 Vh_h=-73;
 k_h=5.6;
 Gnamax=14838*RaINa; // [S/F] maximal INa conductance

//// If Current 
 Gf=90.926*RaIf;
 //x0=-89.73015*Rax0;
 //dx=11.7335*Radx;

//// L-type Ca2+ Current 
 GCaL=0.000175*RaICaL;  // [m^3/F/s] maximal ICaL conductance

//// T-type Ca2+ Current
 GICaT = 183.2*RaICaT; //[S/F] 

//// Transient Outward Current 
 GItoepi=294*RaIto;   // [S/F] maximal ITo conductance
 GItoendo=73;   // [S/F] maximal ITo conductance
 GItoMcell=294; // [S/F] maximal ITo conductance

//// IKs 
 GKsepi   =157*RaIKs; //245; //[S/F] maximal IKs conductance
 GKsendo  =157; //245;// [S/F] maximal IKs conductance
 GKsMcell =40; //62;// [S/F] maximal IKs conductance
 pKNa=0.03;   // [ ]

//// IKr 
 GKr=96*sqrt(Ko/5.4)*RaIKr; //GKr=96 S/F maximal IKr conductance
 //xr1o=0;
 //xr2o=1;
 Q=2.3;
 L0=0.025;
 Kc=0.58e-3;
 Ka=2.6e-3;

//// Inward Rectifier K+ Current
 gK1max=5405*sqrt(Ko/5.4)*RaIK1; //GK1=5405 S/F maximal IK1 conductance

//// Na+/Ca2+ Exchanger Current 
 knaca=1000*RaINaCa;  // [A/F] maximal INaCa
 KmNai=87.5;  // [mM] Nai half saturation constant 
 KmCa=1.38;   // [mM] Cai half saturation constant 
 ksat=0.1;    // [ ]  saturation factor for INaCa 
 n=0.35;      // [ ]  voltage dependence parameter

//// Na+/K+ Pump Current 
 knak=1.362*RaINaK;  // [A/F] maximal INaK
 KmK=1;       // [mM] Ko half saturation constant 
 KmNa=40;     // [mM] Nai half saturation constant 

//// IpCa
 GpCa=0.825*RaICap;    // [A/F] maximal IpCa 
 kpca=0.0005;   // [mM]  Cai half saturation constant 

//// IpK
 GpK=14.6*RaIKp;    // [S/F] maximal IpK conductance

//// Background Currents 
 GbNa=0*0.29*RaIback;   // [S/F] maximal IbNa conductance
 GbCa=0.592*RaIback;  // [S/F] maximal IbCa conductance

///////////////////////////Calcium Dynamics////////////////////

//// Ileak
 Vleak=0.08*RaIleak; // [1/s] maximal Ileak =0.00008/s

//// Irel
 arel=16.464;   // [mM/s] maximal CaSR-dependent Irel
 brel=0.25;     // [mM] CaSR half saturation constant
 crel=8.232;    // [mM/s] maximal CaSR-independent Irel

//// Iup
 Vmaxup=0.425*RaIup; //[mM/s]   // 0.000425;   // [mM/ms] maximal Iup 
 Kup=0.00025;//0.00025;    // [mM] half saturation constant

//Initial values
 V= -0.07;
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

//ODE solving loop/////////////////////////////////////////////////
//cout<<"Loading time variables for ode"<<endl;
dh = 1e-3;		// starting step size for ode solver
tmin = 0.0;			// starting t value (s)
tmax = 10;			// final t simulation value
delta_t = 1e-7;  //time step or max step size, in ms use this parameter to downsample data
int count = 0;
t=tmin;

//for (t = tmin + DT; t<=tmax; t+=DT) {
for (t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
{
    while (t < t_next)	// evolve from t to t_next
    {
        //Reversal potentials
        Eca = 0.5*(R*T/F)*log(Cao/Cai);
        Ek = (R*T/F)*log(Ko/Ki);
        Ena = (R*T/F)*log(Nao/Nai);
        Eks = (R*T/F)*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));

//Gating variables/////////////////////////////////////////////////////////////////////////////////////////////
//INa m gate//////////////////////////////////////////////////////////////
        alpha_m = 1/(1+exp((-60-V*1000)/5));
        beta_m = 0.1/(1+exp((V*1000+35)/5))+0.1/(1+exp((V*1000-50)/200));
        tau_m = alpha_m*beta_m;
        m_inf = 1/((1+exp((-56.86-V*1000)/9.03))*(1+exp((-56.86-V*1000)/9.03)));
        //m = m_inf-(m_inf-m)*exp(-DT*1000/(tau_m));
//INa, h gate /////////////////////////////////////////////////////////////
        switch(celltype){
          case 'e':  h_inf = pow((1/(1+exp((1000*V-Vh_h)/k_h))),0.5);break;
          case 'l': h_inf = 1/((1+exp((V*1000+71.55)/7.43))*(1+exp((V*1000+71.55)/7.43)));break;
        }
                if ( V*1000 >= -40. ){
            alpha_h = 0.;
            beta_h = ((0.77/(0.13*(1+exp(-(V*1000+10.66)/11.1)))));
                }
                else{
            alpha_h = (0.057*exp(-(V*1000+80)/6.8));
            beta_h = (2.7*exp(0.079*V*1000)+(3.1e5)*exp(0.3485*V*1000));
                }
        tau_h = 2.8/(alpha_h+beta_h);
        //h = h_inf-(h_inf-h)*exp(-DT*1000/(tau_h));
        // dhdt = (h_inf-h)/(myCoefTauH/(beta_h+alpha_h))*1000*slo;
//INa j gate///////////////////////////////////////////////////////////////
        //j_inf varies depending on the  cell type////////////////////////////////
        switch(celltype){
          case 'e':   j_inf = pow((1/(1+exp((1000*V-Vh_h)/k_h))),0.5);break;
          case 'l':   j_inf = 1/((1+exp((V*1000+71.55)/7.43))*(1+exp((V*1000+71.55)/7.43)));break;
        }
        if( V*1000 >= -40. ){
            alpha_j = 0.0;
            beta_j = ((0.6*exp((0.057)*V*1000)/(1+exp(-0.1*(V*1000+32)))));
                }
                else{
            alpha_j = (-(25428)*exp(0.2444*V*1000)-(0.000006948)*exp(-0.04391*V*1000))*(V*1000+37.78)/(1+exp(0.311*(V*1000+79.23)));
            beta_j = ((0.02424*exp(-0.01052*V*1000)/(1+exp(-0.1378*(V*1000+40.14)))));
        }
        tau_j = 1.0/(alpha_j+beta_j);
        //j = j_inf-(j_inf-j)*exp(-DT*1000/tau_j);
        //djdt = (j_inf-j)/(myCoefTauJ/(alpha_j+beta_j))*1000*slo;
//ICaL, d gate, and Irel d gate//////////////////////////////////////////
        alpha_d = 1.4/(1+exp((-35-V*1000)/13))+0.25;
        beta_d = 1.4/(1+exp((V*1000+5)/5));
        gamma_d = 1/(1+exp((50-V*1000)/20));
        tau_d = alpha_d*beta_d+gamma_d;
        d_inf = 1/(1+exp(-(V*1000-Vh_dCa)/kCa));
        //dddt = (d_inf-d)/tau_d*1000*slo;
        //d = d_inf-(d_inf-d)*exp(-DT*1000/tau_d);
//ICaL, f gate///////////////////////////////////////////////////////////
        f_inf = 1/(1+exp((V*1000+myVhfCaL)/myKfCaL));
        switch(celltype){
          //Equations for this parameter is a bit confusing must email author
          case 'e': tau_f = 100.;break;
          case 'l':
            if(dfdttemp > 0){
            tau_f = (1125*exp(-((V-myTauFShift)*1000+27)*((V-myTauFShift)*1000+27)/240)+80+165/(1+exp((25-(V-myTauFShift)*1000)/10)))*(1+KtaufCa*(Cai-.5e-4));
                }
                else{
            tau_f = (1125*exp(-((V-myTauFShift)*1000+27)*((V-myTauFShift)*1000+27)/240)+80+165/(1+exp((25-(V-myTauFShift)*1000)/10)));
                }
          break;
        }
        //f = f_inf-(f_inf-f)*exp(-DT*1000/tau_f);
        //dfdt = (f_inf-f)/tau_f*1000*slo;
//ICaL, fCa gate/////////////////////////////////////////////////////// done
        switch(celltype){
          case 'e':   f_ca_inf = (1/(1+(pow(((Cai-myShiftFCaInf)/0.000325),8)))+0.1/(1+exp(((Cai-myShiftFCaInf)-0.0005)/0.0001))+0.2/(1+exp(((Cai-myShiftFCaInf)-0.00075)/0.0008))+0.23)/1.46;break;
          case 'l':   f_ca_inf = (1/(1+(pow((Cai/0.0006),8)))+0.1/(1+exp((Cai-0.0009)/0.0001))+0.3/(1+exp((Cai-0.00075)/0.0008)))/1.3156;break;
        }
        tau_f_ca = 2.0;//ms
         if ( V*1000 > -60.0 ){
            if ( f_ca_inf > f_ca ){constf_ca = 0;}//Note the integral of dfCa/dt = 0 is fCa
            else{constf_ca = 1;}
        }
        else{constf_ca = 1;}

        /*if ( V > -60.0 ){
            if ( f_ca_inf > f_ca ){
                f_ca = f_ca;//Note the integral of dfCa/dt = 0 is fCa
            }
            else{
                f_ca = f_ca_inf-(f_ca_inf-f_ca)*exp(-DT*1000/tau_f_ca);
            }
        }
        else{
            f_ca = f_ca_inf-(f_ca_inf-f_ca)*exp(-DT*1000/tau_f_ca);
        }
        */
        //df_cadt = (f_ca_inf-f_ca)/tau_f_ca*1000*slo*(1-(f_ca_inf>f_ca)*(V>-0.06));
        //cout <<"Calculated gates for ICaL"<<endl;
//Irel, g gate////////////////////////////////////////////////////////// done
        tau_g = 2.0;//units in ms
        if (Cai<=0.00035){
            g_inf = (1/(1+pow((Cai/0.00035),6)));
        }
        else{
            g_inf = (1/(1+pow((Cai/0.00035),16)));
        }
        if ( V*1000 > -60.0 ){
            if ( g_inf > g ){
              constg = 0;
            }
            else{
              constg = 1;
            }
        }
        else{
          constg = 1;
        }//changed to differential expression

       /* if ( V > -60.0 ){
            if ( g_inf > g ){
                g = g;
            }
            else{
                g = g_inf-(g_inf-g)*exp(-DT*1000/tau_g);
            }
        }
        else{
            g = g_inf-(g_inf-g)*exp(-DT*1000/tau_g);
        }
        */
        //dgdt = (g_inf-g)/tau_g*1000*slo*(1-(g_inf>g)*(V>-0.06));
//Ito, r gate/////////////////////////////////////////////////////// done
        r_inf = 1/(1+exp((-1000*V+20+myShiftItoR)/(6*mySlopeItoR)));
        tau_r = myConstTauR*(9.5*exp(-pow((1000*V+40),2)/1800)+0.8);
        //drdt = (r_inf-r)/tau_r*1000*slo;
        //r = r_inf-(r_inf-r)*exp(-DT*1000/tau_r);
//Ito, s gate/////////////////////////////////////////////////////// done
        s_inf = 1/(1+exp((1000*V+20+myShiftItoS)/(5*mySlopeItoS)));
        tau_s = myConstTauS*(85*exp(-(1000*V+45)*(1000*V+45)/320)+5/(1+exp((1000*V-20)/5))+3);
        //dsdt = (s_inf - s)/tau_s *1000*slo;
        //s = s_inf-(s_inf-s)*exp(-DT*1000/tau_s);
        //cout <<"Calculated gates for Ito"<<endl;
//IKs, Xs gate/////////////////////////////////////////////////////// done
        xs_inf = 1/(1+exp((-5-V*1000)/14));
        alpha_xs = 1100/(sqrt(1+exp((-10-V*1000)/6)));
        beta_xs = 1/(1+exp((V*1000-60)/20));
        tau_xs = alpha_xs*beta_xs;
        //dxsdt = ((xs_inf-xs)/(alpha_xs*beta_xs))*slo*1000;
        //xs = xs_inf-(xs_inf-xs)*exp(-DT*1000/tau_xs);
        //cout <<"Calculated gates for IKs"<<endl;
//IKr, Xr1 gate/////////////////////////////////////////////////////// done
        //parameter added for rapid delayed rectifier current
        xr1_inf = 1/(1+exp((((-R*T/F/Q*log(1/pow(((1+Cao*0.001/Kc)/(1+Cao*0.001/Ka)),4)/L0))-0.026)*1000-(V-myShift)*1000)/7));
        alpha_xr1 = 450/(1+exp((-45-(V-myShift)*1000)/(10)));
        beta_xr1 = 6/(1+exp(((V-myShift)*1000-(-30))/11.5));
        tau_xr1 = alpha_xr1*beta_xr1;
        //dxr1dt = ((xr1_inf-xr1)/(myRedTauxr1*alpha_xr1*beta_xr1 ))*1000;
        //xr1 = xr1_inf-(xr1_inf-xr1)*exp(-DT*1000/tau_xr1);
//IKr, Xr2 gate////////////////////////////////////////////////////// done
        xr2_inf = 1/(1+exp(((V-myShift)*1000-(-88))/24));
        alpha_xr2 = 3/(1+exp((-60-(V-myShift)*1000)/20));
        beta_xr2 = 1.12/(1+exp(((V-myShift)*1000-60)/20));
        tau_xr2 = alpha_xr2*beta_xr2;
        //dxr2dt = ((xr2_inf-xr2)/(myRedTauxr2*alpha_xr2*beta_xr2))*1000;
        //xr2 = xr2_inf-(xr2_inf-xr2)*exp(-DT*1000/tau_xr2);
//If, Xf gate////////////////////////////////////////////////////////done
        tau_xf = 1900;//ms
        xf_inf = 1/(1+exp((V*1000-(-102.4))/(7.6)));
        //xf_inf = 1/(1 + exp((u*1000-x0)/dx));
        //dxfdt=(xf_inf-xf)/tauif*1000;
        //xf = xf_inf-(xf_inf-xf)*exp(-DT*1000/tau_xf);
        //xf = xf + DT*((xf_inf-xf)/tau_xf*1000);
//ICaT, dCaT gate//////////////////////////////////////////////////// done
        dCaT_inf = 1/(1+exp(-(1000*V+26.3)/(6)));
        tau_dCaT = 1/(1.068*exp((1000*V+26.3)/(30))+1.068*exp(-(1000*V+26.3)/(30)));
        //ddCaTdt = (dCaT_inf-dCaT)/tau_dCaT*1000*slo;
        //dCaT = dCaT_inf-(dCaT_inf-dCaT)*exp(-DT*1000/tau_dCaT);
        //dCaT = dCaT + DT *((dCaT_inf - dCaT)/tau_dCaT*1000);
//ICaT, fCaT gate//////////////////////////////////////////////////// done
        fCaT_inf = 1/(1+exp((1000*V+61.7)/(5.6)));
        tau_fCaT = 1/(0.0153*exp(-(1000*V+61.7)/(83.3))+ 0.015*exp((1000*V+61.7)/(15.38)));
        //dfCaTdt = (fCaT_inf-fCaT)/tau_fCaT*1000*slo;
        //fCaT = fCaT_inf-(fCaT_inf-fCaT)*exp(-DT*1000/tau_fCaT);
        //fCaT = fCaT + DT *((fCaT_inf - fCaT)/tau_fCaT*1000);
//IK1, alphas and betas////////////////////////////////////////////// done
        alpha_K1 = 0.1/(1+exp(0.06*((V-myShiftK1)*1000-Ek*1000-200)));
        beta_K1 = (3*exp(0.0002*((V-myShiftK1)*1000-Ek*1000+100))+exp(0.1*((V-myShiftK1)*1000-Ek*1000-10)))/(1+exp(-0.5*((V-myShiftK1)*1000-Ek*1000)));
        x_K1_inf = alpha_K1/(alpha_K1+beta_K1);
///Current equations///////////////////////////////////////////////////////////////////////
        //Na+ current, INa done
        INa = Gnamax*m*m*m*h*j*(V-Ena);
        //L-type Ca2+ current, ICaL done
        ICaL = GCaL*d*f*f_ca*4*V*(F*F)/R/T*(Cai*exp(2*V*F/R/T)-0.341*Cao)/(exp(2*V*F/R/T)-1);
        //Transient outward current, Ito done
        Ito = GItoepi*r*s*(V-Ek);
        //Rapid delayed rectifier K+ current, IKr done
        IKr = GKr*xr1*xr2*(V-Ek);
        //Slow delayed rectifier K+ current, IKs, added a scaling factor (1+0.6/(1+pow((3.8e-5/Cai),1.4))
        IKs = GKsepi*(1+.6/pow((1+(3.8e-5/Cai)),1.4))*xs*xs*(V-Eks);
        //Inward rectifier K+ current, IK1 done
        IK1=x_K1_inf*gK1max*(V-Ek);
        //Na+/Ca2+ exchanger current, INaCa, not than n represents lambda done
        INaCa = knaca*(1/(KmNai*KmNai*KmNai+Nao_naca*Nao_naca*Nao_naca))*(1/(KmCa+Cao))*(1/(1+ksat*exp((n-1)*V*F/(R*T))))*(exp(n*V*F/(R*T))*Nai*Nai*Nai*Cao-exp((n-1)*V*F/(R*T))*Nao_naca*Nao_naca*Nao_naca*Cai*alfa);
        //Na+/K+ pump current, INaK, note that code uses knak variable for pnak
        INaK =(1/(1+0.1245*exp(-0.1*V*F/(R*T))+0.0353*exp(-V*F/(R*T))))*(knak*Ko/(Ko+KmK)*Nai/(Nai+KmNa));
        //Calcium pump current, IpCa done
        IpCa = GpCa*Cai/(kpca+Cai);
        //Plateau K+ current done
        IpK = GpK*(V-Ek)*(1/(1+exp((25-V*1000)/5.98)));
        //Background sodium current done
        IbNa = GbNa*(V-Ena);
        //Calcium background current, IbCa done
        IbCa = GbCa*(V-Eca);
        //Currents added for the hESC_CM model
        ////Hyperpolarization activated funny current, If note this is a new current incorporated into the model
        //not initially in the ten Tusscher model done
        If = Gf*xf*(V+.017);
        //T-type Ca2+ Current, ICaT, not originally in the tentusscher model done
        ICaT = dCaT*fCaT*GICaT*(V-Eca);
        Itotal = INa + ICaL + Ito + IKr + IKs + IK1 + INaCa + INaK + IpCa + IbNa + IpK + IbCa + If + ICaT;
///Calcium dynamics//////////////////////////////////////////////////////////////////////////////////////////
        Ileak = Vleak*(CaSR-Cai);
        Iup = Vmaxup/(1.+((Kup*Kup)/(Cai*Cai)));
        //modified by RaIrel
        Irel = (d*g*(crel+arel*(CaSR*CaSR)/((brel*brel)+(CaSR*CaSR))))*RaIrel;
      /*  CaBuf = Bufc*Cai/(Cai+Kbufc);
        CaCSQN = Bufsr*CaSR/(CaSR+Kbufsr);
        CaSRCurrent = Iup-Irel-Ileak;
        //Added ICaT to CaCurrent
        CaCurrent = -(ICaL+IbCa+ICaT+IpCa-2.0*INaCa)*(1.0/(2.0*Vc*F))*capacitance;
        dCaSR = DT*(Vc/Vsr)*CaSRCurrent;
        bjsr = Bufsr-CaCSQN-dCaSR-CaSR+Kbufsr;
        cjsr = Kbufsr*(CaCSQN+dCaSR+CaSR);
        CaSR = (sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
        dCai = DT*(CaCurrent-CaSRCurrent);
        bc = Bufc-CaBuf-dCai-Cai+Kbufc;
        cc = Kbufc*(CaBuf+dCai+Cai);
        Cai = (sqrt(bc*bc+4*cc)-bc)/2;
        */
         //adding Caibufc variable
        Caibufc = 1./(1.+Bufc*Kbufc/((Cai + Kbufc)*(Cai + Kbufc)));
        //adding CaSRbufsr variable
        CaSRbufsr = 1./(1.+Bufsr*Kbufsr/((CaSR + Kbufsr)*(CaSR + Kbufsr)));

///////////////////////////////////////////////////////////////////////////////////////////
        //V = V+DT*(-Itotal);
        int status = gsl_odeiv2_evolve_apply (evolve_ptr, control_ptr, step_ptr, &my_system, &t, t_next, &dh, y);
        
        //int status = gsl_odeiv2_driver_apply (driver_ptr, &t, t_next, y);
        if (status != GSL_SUCCESS)
        {
            printf ("error, return value=%d\n", status);
            break;//note only breaks a single loop
        }

///Print screen and print to file///////////////////////////////////////////////////////
        file_output_counter++;
        if(file_output_counter%10000==0)
        {
            //printf( "\nV = %f;\nCai = %f;\nCaSR = %f;\nm = %f;\nh = %f;\nj = %f;\nxr1 = %f;\nxr2 = %f;\nxs = %f;\nr = %f;\ns = %f;\nd = %f;\nf = %f;\nf_ca = %f;\ng = %f;\ndCaT = %f;\nfCaT = %f;\nxf = %f;\n", V, Cai, CaSR, m, h, j, xr1, xr2, xs, r, s, d, f, f_ca,g, dCaT, fCaT, xf );
            //fprintf( output, "%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\n", time+(number-1)*basic_cycle_length, V, m, h, j, d, f, g, f_ca, r, s, xs, xr1, xr2 );
            //fprintf( output, "%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.105f\t%4.10f\t%4.10f\t%4.10f\n", time+(number-1)*basic_cycle_length, Ileak, Iup, Irel, CaBuf, CaSR, Cai, Nai, Ki, V );
            fprintf( output, "%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\n", t, IKr, IKs, IK1, Ito, INa, If, INaK, ICaL, IbCa, INaCa, V,ICaT, Cai, IpCa, Irel, Itotal );
        }
        // printf( "\nV = %f;\nCai = %f;\nCaSR = %f;\nm = %f;\nh = %f;\nj = %f;\nxr1 = %f;\nxr2 = %f;\nxs = %f;\nr = %f;\ns = %f;\nd = %f;\nf = %f;\nf_ca = %f;\ng = %f;\ndCaT = %f;\nfCaT = %f;\nxf = %f;\n", V, Cai, CaSR, m, h, j, xr1, xr2, xs, r, s, d, f, f_ca,g, dCaT, fCaT, xf );
    }
}
    printf( "\nV = %f;\nCai = %f;\nCaSR = %f;\nm = %f;\nh = %f;\nj = %f;\nxr1 = %f;\nxr2 = %f;\nxs = %f;\nr = %f;\ns = %f;\nd = %f;\nf = %f;\nf_ca = %f;\ng = %f;\ndCaT = %f;\nfCaT = %f;\nxf = %f;\n", V, Cai, CaSR, m, h, j, xr1, xr2, xs, r, s, d, f, f_ca,g, dCaT, fCaT, xf );
    

//free up the gsl_odeiv stuff
gsl_odeiv2_evolve_free (evolve_ptr);
gsl_odeiv2_control_free (control_ptr);
gsl_odeiv2_step_free (step_ptr);
gsl_odeiv2_driver_free (driver_ptr);
fclose( output );
    
    return 0;

}
static int nstates(double t, const double y[], double dy[], void *params_ptr){
    struct ode_params *odeParams = (struct ode_params *) params_ptr;
    // //cout<<"Load params structure"<<endl;
    dVdt = -Itotal;
    dmdt = (m_inf-m)/tau_m*1000;
    dhdt = (h_inf-h)/tau_h*1000;
    djdt = (j_inf-j)/tau_j*1000;
    dxr1dt = (xr1_inf-xr1)/tau_xr1*1000;
    dxr2dt = (xr2_inf-xr2)/tau_xr2*1000;
    dxsdt = (xs_inf-xs)/tau_xs*1000;
    drdt = (r_inf-r)/tau_r*1000;
    dsdt = (s_inf-s)/tau_s*1000;
    dddt = (d_inf-d)/tau_d*1000;
    dfdt = (f_inf-f)/tau_f*1000;
    dfdttemp = dfdt;
    df_cadt = constf_ca*(f_ca_inf - f_ca)/tau_f_ca*1000;
    dgdt = constg*(g_inf-g)/tau_g*1000;
    ////cout<<"constg"<<constg<<endl;
    dxfdt = (xf_inf-xf)/tau_xf*1000;
    ddCaTdt = (dCaT_inf-dCaT)/tau_dCaT*1000;
    dfCaTdt = (fCaT_inf-fCaT)/tau_fCaT*1000;
    dCaidt = Caibufc*(Ileak-Iup+Irel-((ICaL+IbCa+ICaT+IpCa-2.0*INaCa)*(1.0/(2.0*Vc*F)))*CAP);
    dCaSRdt = ((CaSRbufsr*Vc)/Vsr)*(Iup-(Irel +Ileak));

    return GSL_SUCCESS;		// GSL_SUCCESS defined in gsl/errno.h as 0
}

