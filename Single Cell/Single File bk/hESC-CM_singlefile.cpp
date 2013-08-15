#include <cmath>
#include <iostream>
#include <cstdio>

//Note model was converted from V/s to mV/ms by converting:
// myShiftK1 *1000 to convert from V to mV
// R constant  *1000 to convert the  reversal potentials in mV 
// Maximal conductances from S/F to mS/F or nS/pF
// Time dependent only parameters
// Vleak *(1/1000) to convert 1/s to 1/ms units
// arel and crel *(1/1000) to convert mM/s to mM/ms units
// Vc and Vsr converted from um^3 to picoliters
// To reproduce simulations in paper experiments were run for 350 seconds,
// drug block was administered at 300 seconds.
#define beats 1000
#define EPI
#define INITIALVALS
#define EARLY
//#define LATE


using namespace std;


//Variables to measure action potential parameters
int counter = 0;
int start_output = 0;
double vmin[beats];
double tvmin[beats];
double vmax[beats];
double dvdtmax[beats];
//double vdvdtmax[beats];
double apd_start[beats];
double apd30_start[beats];
double apd50_start[beats];
double apd70_start[beats];
double apd90_start[beats];
double APA[beats];
double apd_end[beats];
double ddr[beats];
double top[beats];
double ttop[beats];
double top_slope[beats];
double apd30[beats];
double apd50[beats];
double apd70[beats];
double apd90[beats];
double bcl[beats];

//timestep
double DT = 0.001;//ms

bool Ikrblock = true;
double blocktime = 300000;
double E4031 = 1;
double tmin = 0.0;
double tmax = 350000;
double t;
char  output_file_name[20] = "Ikrblock_early.dat";


#ifdef EARLY
int CellTypeFlag= 2;
double cardiom_C = 41e-12; //[F] spostata dal file parametri, è il valore approssimato della C del cardiom per le equazioni dei fibro.
double EarlySRTypeFlag = 2;
double lock_Nai=0;
double lock_Ki = 0;
double Cap_div=1;
double RaINa=0.038;
double myCoefTauH =2.8;//3.7;
double myCoefTauJ = 1;
double RaICaL=0.25;
double Vh_dCa=12.5;
double kCa=12.5;
double ampICaLDinf = 1;
double KtaufCa=1433;   // [1/mM]  Altamirano & Bers 2007
double myVhfCaL = 20;
double myKfCaL = 7;//7
double ampICaLFinf = 1;//per accelerare la tauF e far inattivare piu' rapidamente la ICaL
double myKtauF = 1;
double myTauFShift = 0;
double myShiftFCaInf = -0.11e-3;  //[mM] per spostare un po'  a sx la signmoide di fCa
double Vth_ICaL = -60;
double RaICaT =0.25;
double RaINaCa=1.750e1;
double alfa=0.8;
double RaINaK=0.7;//0.79;    ////////RETICOLO ON
double RaIK1=0.05*2.67/3;  //// Sartiani 20100412 modificato per prove
double myShiftK1 = -0.015*1000;
double RaIKr=3;//3;//2.6;//2;//4;//2*2;          //// RETICOLO ON
double poffi = 0;
double mySlope = 1;//0.5;
double myShift = 0;
double RaIKs=0.1;
double RaIto=0.1673*0.4903*0.8;
double myShiftItoR = -25;//-20;
double mySlopeItoR = 0.3;
double myShiftItoS = 0;//-5;
double mySlopeItoS = 1;
double myConstTauR= 1;
double myConstTauS= 1;
double RaIup=0.4/3;//0.7/3;//0.3/18;//0.003*2; //0.03*10;  ////RETICOLO ON
double RaIrel=0.2/18;//0.2/18;//0.005*2;//0.05*10;
double RaIleak=0.1/18;//0.4/18;//0.004*2;//0.04*10;
double RaIf=0.5389;//0.5389;     ////Sartiani
double Rax0=1.1061;
double Radx=0.6537;
double myIfconst = 1;
double RaCm=0.22162;
double RaICap=1;
double RaIKp=0;
double RaIback=0.2;//0.2;//0.45;
#endif

#ifdef LATE
int CellTypeFlag= 3;
double ryanSR = 1;
double cardiom_C = 33e-12;// [F] spostata dal file parametri, è il valore approssimato della C del cardiom per le equazioni dei fibro.
double lock_Nai=0;
double lock_Ki = 0;
double Cap_div=1;
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
double myShiftFCaInf = -0.12e-3; //[mM] per spostare un po'  a sx la signmoide di fCa
double Vth_ICaL = 0;
double RaICaT =0.05;
double RaINaCa=1.824e+01;     //////////////////// dai dati Sartiani
double alfa=0.38;
double RaINaK=0.83;
double RaIK1=0.4*0.05*2.67*4;
double myShiftK1 = -0.015*1000;
double RaIKr=1.4;
double poffi = 0;
double mySlope = 1;
double myShift = 0;
double RaIKs=0.1;
double RaIto=0.3754*0.4903*0.9;         // dai dati Sartiani
double myShiftItoR = -25;
double mySlopeItoR = 0.3;
double myShiftItoS = 0;
double mySlopeItoS = 1;
double myConstTauR= 1;
double myConstTauS= 1;
double RaIup=0.33;
double RaIrel=0.4;
double RaIleak=0.3*1;
double RaIf=0.23;            // dati Sartiani
double Rax0=1.1061;
double Radx=0.6537;
double RaCm=0.17838; // dati Sartiani
double RaICap=1;
double RaIKp=0;
double RaIback=1;  // Itoh
#endif

//// Constants
double R=8314.472;   // [J/millimoles/K] Gas constant
double F=96485.3415; // [C/mol]	  Faraday constant
double T=310.0;      // [K]       Temperature

//// Buffering
double Bufc=0.25;   // [mM] total cytoplasmic buffer concentration
double Kbufc=0.001;  // [mM] Cai half saturation constant
double Bufsr=10;     // [mM] total sarcoplasmic buffer concentration
double Kbufsr=0.3;   // [mM] CaSR half saturation constant

//// Extracellular Ionic concentrations
////TT
//Ko=5.4;      // [mM]
//Cao=1.8;     // [mM]
//Nao=140;     // [mM]

////Sartiani
double Ko=4;         // [mM]
double Cao=2.7;      // [mM]
double Nao=150.5;    // [mM]
double Nao_naca = Nao;



//// Intracellular Ionic concentrations
// // Pre-dialysis
double Ki=140;      // [mM]  & 140 in TT04
double Nai=7;    // [mM]

//// Intracellular Volumes
double Vc=16.404*RaCm;
double Vsr=1.094*RaCm;
double capacitance=0.185*1000*RaCm;


//// Flag to choose between epi, endo and M cell types
int epi=1;
int endo=0;
int Mcell=0;

//// Ionic Currents

//// Fast Na+ Current
double Vh_h=-73;
double k_h=5.6;
double Gnamax=14.838*RaINa; // [nS/pF] maximal INa conductance


//// If Current
double Gf=0.090926*RaIf;
double x0=-89.73015*Rax0;
double dx=11.7335*Radx;


//// L-type Ca2+ Current

double GCaL=0.000175*RaICaL;  // [m^3/F/s] maximal ICaL conductance

//// T-type Ca2+ Current
double GICaT = 0.1832*RaICaT; //[S/F]


//// Transient Outward Current
double GItoepi=0.294*RaIto;   // [S/F] maximal ITo conductance
double GItoendo=73;   // [S/F] maximal ITo conductance
double GItoMcell=294; // [S/F] maximal ITo conductance
int soepi=1;
int soendo=1;


//// IKs
double GKsepi   =0.157*RaIKs; //245; //[S/F] maximal IKs conductance
double GKsendo  =157; //245;// [S/F] maximal IKs conductance
double GKsMcell =40; //62;// [S/F] maximal IKs conductance
double pKNa=0.03;   // [ ]

//// IKr

double GKr=0.096*sqrt(Ko/5.4)*RaIKr; //GKr=96 nS/pF maximal IKr conductance
double Q=2.3;
double L0=0.025;
double Kc=0.58e-3;
double Ka=2.6e-3;


//// Inward Rectifier K+ Current
double gK1max=5.405*sqrt(Ko/5.4)*RaIK1; // maximal IK1 conductance

//// Na+/Ca2+ Exchanger Current
double knaca=1000*RaINaCa;  // [pA/pF] maximal INaCa
double KmNai=87.5;  // [mM] Nai half saturation constant
double KmCa=1.38;   // [mM] Cai half saturation constant
double ksat=0.1;    // [dimensionless]  saturation factor for INaCa
double n=0.35;      // [dimensionless]  voltage dependence parameter

//// Na+/K+ Pump Current
double knak=1.362*RaINaK;  // [pA/pF] maximal INaK
double KmK=1;       // [mM] Ko half saturation constant
double KmNa=40;     // [mM] Nai half saturation constant

//// IpCa
double GpCa=0.825*RaICap;    // [pA/pF] maximal IpCa
double kpca=0.0005;   // [mM]  Cai half saturation constant

//// IpK
double GpK=0.0146*RaIKp;    // [nS/F] maximal IpK conductance

//// Background Currents
double GbNa=0*0.29*RaIback;   // [nS/pF] maximal IbNa conductance
double GbCa=0.000592*RaIback;  // [nS/pF] maximal IbCa conductance

//// Calcium Dynamics
//// Ileak
double Vleak=0.00008*RaIleak; // [1/ms] maximal Ileak =0.00008/s

//// Irel
double arel=0.016464;   // [mM/ms] maximal CaSR-dependent Irel
double brel=0.25;     // [mM] CaSR half saturation constant
double crel=0.008232;    // [mM/ms] maximal CaSR-independent Irel

//// Iup
double Vmaxup=0.000425*RaIup; //[mM/ms]   // 0.000425;   // [mM/ms] maximal Iup
double Kup=0.00025;//0.00025;    // [mM] half saturation constant

//APD paramters
double dvdt,dvdtold;
double vnew;

int main(){
    FILE *output;
    int file_output_counter=0;
    double RTONF, Ek, Ena, Eks, Eca, Ef;
    
    double V, Cai, CaSR, alpha_K1, beta_K1, x_K1_inf;
    double alpha_m, beta_m, tau_m, m_inf, alpha_h, beta_h, tau_h, h_inf, alpha_j, beta_j, tau_j, j_inf, alpha_d, beta_d, gamma_d, tau_d, d_inf, f_inf, tau_f, g_inf, tau_g, constg, f_ca_inf, tau_f_ca, constf_ca, r_inf, tau_r, s_inf, tau_s, alpha_xs, beta_xs, xs_inf, tau_xs, alpha_xr1, beta_xr1, xr1_inf, tau_xr1, alpha_xr2, beta_xr2, xr2_inf, tau_xr2, xf_inf, tau_xf, dCaT_inf, tau_dCaT, fCaT_inf, tau_fCaT;
    double INa, m, h, j, ICaL, d, df, tempf, f, g, f_ca, Ito, r, s, IKr, xr1, xr2, IKs, xs, IK1, INaCa, INaK, IpCa, IpK, IbNa, IbCa, Istim,If,xf,ICaT,dCaT,fCaT, Itotal;
    double Ileak, Iup, Irel, CaBuf, CaCSQN, CaCurrent, CaSRCurrent, dCaSR, bjsr, cjsr, dCai, bc, cc;
    
    output = fopen( output_file_name , "w" );
    

//Initial values
#ifdef INITIALVALS
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
#endif

#ifndef INITIALVALS
//Initial Conditions after 250 seconds ~100 beats 
#ifdef EARLY
V = -0.069160;
Cai = 0.000027;
CaSR = 0.086301;
m = 0.041569;
h = 0.600505;
j = 0.618972;
xr1 = 0.001208;
xr2 = 0.313321;
xs = 0.010087;
r = 0.000000;
s = 0.999948;
d = 0.001452;
f = 0.999173;
f_ca = 1.002065;
g = 1.000000;
dCaT = 0.000789;
fCaT = 0.797842;
xf = 0.014191;
#endif

#ifdef LATE
//Initial Conditions after 100 seconds
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
#endif
#endif

RTONF = (R*T)/F;
for(counter=0;counter<beats;counter++){
    vmin[counter] = 100000.0;
    vmax[counter] = -10000.0;
    dvdtmax[counter]       = -10000.0;
    ddr[counter]           = -10000.0;
    APA[counter]           =  10000.0;
    top[counter]           =  10000.0;
    top_slope[counter]     = -10000.0;
    apd30[counter]         = -10000.0;
    apd50[counter]         = -10000.0;
    apd70[counter]         = -10000.0;
    apd90[counter]         = -10000.0;
    bcl[counter]  = -10000.0;
}
printf("Vmin(mV)\tVmax(mV)\tdVdtmax(mV/s)\tAPD30(ms)\tAPD50(ms)\tAPD70(ms)\tAPD90(ms)\tFreq.(bpm)\tDDR(mV/s)\tTOP(mV)\tAPA(mV)\n");
counter =0;
start_output = 0;
dvdt = 1000.0;
dvdtold = 500.0;
start_output = 0;
t=tmin;

    for (t = tmin + DT; t<=tmax; t+=DT) {

        if ((Ikrblock) && (t > (blocktime))){
          #ifdef EARLY
            E4031 = 0.5;
          #endif
          #ifdef LATE
            E4031 = 0.4;
          #endif
        }

        //Reversal potentials
        Eca = 0.5*(R*T/F)*log(Cao/Cai);
        Ek = (R*T/F)*log(Ko/Ki);
        Ena = (R*T/F)*log(Nao/Nai);
        Eks = (R*T/F)*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));

///////////////////////////////////////////////////////////////////////////////////////////////
//INa m gate//////////////////////////////////////////////////////////////
        alpha_m = 1/(1+exp((-60-V)/5));
        beta_m = 0.1/(1+exp((V+35)/5))+0.1/(1+exp((V-50)/200));
        tau_m = alpha_m*beta_m;
        m_inf = 1/((1+exp((-56.86-V)/9.03))*(1+exp((-56.86-V)/9.03)));
        m = m_inf-(m_inf-m)*exp(-DT/(tau_m));
//INa, h gate /////////////////////////////////////////////////////////////
        #ifdef EARLY
                h_inf = pow((1/(1+exp((V-Vh_h)/k_h))),0.5);
        #endif
        #ifdef LATE
                h_inf = 1/((1+exp((V+71.55)/7.43))*(1+exp((V+71.55)/7.43)));
        #endif
            if ( V >= -40. ){
            alpha_h = 0.;
            beta_h = ((0.77/(0.13*(1+exp(-(V+10.66)/11.1)))));
            }
            else{
            alpha_h = (0.057*exp(-(V+80)/6.8));
            beta_h = (2.7*exp(0.079*V)+(3.1e5)*exp(0.3485*V));
            }
        tau_h = 2.8/(alpha_h+beta_h);
        h = h_inf-(h_inf-h)*exp(-DT/(tau_h));
        // dhdt = (h_inf-h)/(myCoefTauH/(beta_h+alpha_h))*1000*slo;
//INa j gate///////////////////////////////////////////////////////////////
        //j_inf varies depending on the  cell type////////////////////////////////
        #ifdef EARLY
                j_inf = pow((1/(1+exp((V-Vh_h)/k_h))),0.5);
        #endif
        #ifdef LATE
                j_inf = 1/((1+exp((V+71.55)/7.43))*(1+exp((V+71.55)/7.43)));
        #endif
            if( V >= -40. ){
            alpha_j = 0.0;
            beta_j = ((0.6*exp((0.057)*V)/(1+exp(-0.1*(V+32)))));
            }
            else{
            alpha_j = (-(25428)*exp(0.2444*V)-(0.000006948)*exp(-0.04391*V))*(V+37.78)/(1+exp(0.311*(V+79.23)));
            beta_j = ((0.02424*exp(-0.01052*V)/(1+exp(-0.1378*(V+40.14)))));
            }
        tau_j = 1.0/(alpha_j+beta_j);
        j = j_inf-(j_inf-j)*exp(-DT/tau_j);
        //djdt = (j_inf-j)/(myCoefTauJ/(alpha_j+beta_j))*1000*slo;
//ICaL, d gate, and Irel d gate//////////////////////////////////////////
        alpha_d = 1.4/(1+exp((-35-V)/13))+0.25;
        beta_d = 1.4/(1+exp((V+5)/5));
        gamma_d = 1/(1+exp((50-V)/20));
        tau_d = alpha_d*beta_d+gamma_d;
        d_inf = 1/(1+exp(-(V-Vh_dCa)/kCa));
        //dddt = (d_inf-d)/tau_d*1000*slo;
        d = d_inf-(d_inf-d)*exp(-DT/tau_d);
//ICaL, f gate///////////////////////////////////////////////////////////
        f_inf = 1/(1+exp((V+myVhfCaL)/myKfCaL));
        #ifdef EARLY
                //Equations for this parameter is a bit confusing must email author
                tau_f = 100.;
        #endif
        #ifdef LATE
            if(f_inf > f){
            tau_f = (1125*exp(-((V-myTauFShift)+27)*((V-myTauFShift)+27)/240)+80+165/(1+exp((25-(V-myTauFShift))/10)))*(1+KtaufCa*(Cai-.5e-4));
            }
            else{
            tau_f = (1125*exp(-((V-myTauFShift)+27)*((V-myTauFShift)+27)/240)+80+165/(1+exp((25-(V-myTauFShift))/10)));
            }
        #endif
        f = f_inf-(f_inf-f)*exp(-DT/tau_f);
        //dfdt = (f_inf-f)/tau_f*1000*slo;
//ICaL, fCa gate/////////////////////////////////////////////////////// done
        #ifdef EARLY
                f_ca_inf = (1/(1+(pow(((Cai-myShiftFCaInf)/0.000325),8)))+0.1/(1+exp(((Cai-myShiftFCaInf)-0.0005)/0.0001))+0.2/(1+exp(((Cai-myShiftFCaInf)-0.00075)/0.0008))+0.23)/1.46;
        #endif
        #ifdef LATE
                f_ca_inf = (1/(1+(pow((Cai/0.0006),8)))+0.1/(1+exp((Cai-0.0009)/0.0001))+0.3/(1+exp((Cai-0.00075)/0.0008)))/1.3156;
        #endif
                tau_f_ca = 2.0;//ms
        if ( V > -60.0 ){
            if ( f_ca_inf > f_ca ){
                f_ca = f_ca;//Note the integral of dfCa/dt = 0 is fCa
            }
            else{
                f_ca = f_ca_inf-(f_ca_inf-f_ca)*exp(-DT/tau_f_ca);
            }
        }
        else{
            f_ca = f_ca_inf-(f_ca_inf-f_ca)*exp(-DT/tau_f_ca);
        }
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
        if ( V > -60.0 ){
            if ( g_inf > g ){
                g = g;
            }
            else{
                g = g_inf-(g_inf-g)*exp(-DT/tau_g);
            }
        }
        else{
            g = g_inf-(g_inf-g)*exp(-DT/tau_g);
        }
        //dgdt = (g_inf-g)/tau_g*1000*slo*(1-(g_inf>g)*(V>-0.06));
//Ito, r gate/////////////////////////////////////////////////////// done
        r_inf = 1/(1+exp((-V+20+myShiftItoR)/(6*mySlopeItoR)));
        tau_r = myConstTauR*(9.5*exp(-pow((V+40),2)/1800)+0.8);
        //drdt = (r_inf-r)/tau_r*1000*slo;
        r = r_inf-(r_inf-r)*exp(-DT/tau_r);
//Ito, s gate/////////////////////////////////////////////////////// done
        s_inf = 1/(1+exp((V+20+myShiftItoS)/(5*mySlopeItoS)));
        tau_s = myConstTauS*(85*exp(-(V+45)*(V+45)/320)+5/(1+exp((V-20)/5))+3);
        //dsdt = (s_inf - s)/tau_s *1000*slo;
        s = s_inf-(s_inf-s)*exp(-DT/tau_s);
        //cout <<"Calculated gates for Ito"<<endl;
//IKs, Xs gate/////////////////////////////////////////////////////// done
        xs_inf = 1/(1+exp((-5-V)/14));
        alpha_xs = 1100/(sqrt(1+exp((-10-V)/6)));
        beta_xs = 1/(1+exp((V-60)/20));
        tau_xs = alpha_xs*beta_xs;
        //dxsdt = ((xs_inf-xs)/(alpha_xs*beta_xs))*slo*1000;
        xs = xs_inf-(xs_inf-xs)*exp(-DT/tau_xs);
        //cout <<"Calculated gates for IKs"<<endl;
//IKr, Xr1 gate/////////////////////////////////////////////////////// done
        //parameter added for rapid delayed rectifier current
        xr1_inf = 1/(1+exp((((-R*T/F/Q*log(1/pow(((1+Cao*0.001/Kc)/(1+Cao*0.001/Ka)),4)/L0))-26)-(V-myShift))/7));
        alpha_xr1 = 450/(1+exp((-45-(V-myShift))/(10)));
        beta_xr1 = 6/(1+exp(((V-myShift)-(-30))/11.5));
        tau_xr1 = alpha_xr1*beta_xr1;
        //dxr1dt = ((xr1_inf-xr1)/(myRedTauxr1*alpha_xr1*beta_xr1 ))*1000;
        xr1 = xr1_inf-(xr1_inf-xr1)*exp(-DT/tau_xr1);
//IKr, Xr2 gate////////////////////////////////////////////////////// done
        xr2_inf = 1/(1+exp(((V-myShift)-(-88))/24));
        alpha_xr2 = 3/(1+exp((-60-(V-myShift))/20));
        beta_xr2 = 1.12/(1+exp(((V-myShift)-60)/20));
        tau_xr2 = alpha_xr2*beta_xr2;
        //dxr2dt = ((xr2_inf-xr2)/(myRedTauxr2*alpha_xr2*beta_xr2))*1000;
        xr2 = xr2_inf-(xr2_inf-xr2)*exp(-DT/tau_xr2);
//If, Xf gate////////////////////////////////////////////////////////done
        tau_xf = 1900;//ms
        xf_inf = 1/(1+exp((V-(-102.4))/(7.6)));
        //xf_inf = 1/(1 + exp((u*1000-x0)/dx));
        //dxfdt=(xf_inf-xf)/tauif*1000;
        xf = xf_inf-(xf_inf-xf)*exp(-DT/tau_xf);
        //xf = xf + DT*((xf_inf-xf)/tau_xf*1000);
//ICaT, dCaT gate//////////////////////////////////////////////////// done
        dCaT_inf = 1/(1+exp(-(V+26.3)/(6)));
        tau_dCaT = 1/(1.068*exp((V+26.3)/(30))+1.068*exp(-(V+26.3)/(30)));
        //ddCaTdt = (dCaT_inf-dCaT)/tau_dCaT*1000*slo;
        dCaT = dCaT_inf-(dCaT_inf-dCaT)*exp(-DT/tau_dCaT);
        //dCaT = dCaT + DT *((dCaT_inf - dCaT)/tau_dCaT*1000);
//ICaT, fCaT gate//////////////////////////////////////////////////// done
        fCaT_inf = 1/(1+exp((V+61.7)/(5.6)));
        tau_fCaT = 1/(0.0153*exp(-(V+61.7)/(83.3))+ 0.015*exp((V+61.7)/(15.38)));
        //dfCaTdt = (fCaT_inf-fCaT)/tau_fCaT*1000*slo;
        fCaT = fCaT_inf-(fCaT_inf-fCaT)*exp(-DT/tau_fCaT);
        //fCaT = fCaT + DT *((fCaT_inf - fCaT)/tau_fCaT*1000);
//IK1, alphas and betas////////////////////////////////////////////// done
        alpha_K1 = 0.1/(1+exp(0.06*((V-myShiftK1)-Ek-200)));
        beta_K1 = (3*exp(0.0002*((V-myShiftK1)-Ek+100))+exp(0.1*((V-myShiftK1)-Ek-10)))/(1+exp(-0.5*((V-myShiftK1)-Ek)));
        x_K1_inf = alpha_K1/(alpha_K1+beta_K1);
        
        //Na+ current, INa done
        INa = Gnamax*m*m*m*h*j*(V-Ena);
        
        //L-type Ca2+ current, ICaL done
        ICaL = GCaL*d*f*f_ca*4*V*(F*F)/R/T*(Cai*exp(2*V*F/R/T)-0.341*Cao)/(exp(2*V*F/R/T)-1);
        
        //Transient outward current, Ito done
        Ito = GItoepi*r*s*(V-Ek);
        
        //Rapid delayed rectifier K+ current, IKr done
        IKr = E4031*GKr*xr1*xr2*(V-Ek);
        
        //Slow delayed rectifier K+ current, IKs, added a scaling factor (1+0.6/(1+pow((3.8e-5/Cai),1.4))
        IKs = GKsepi*(1+.6/pow((1+(3.8e-5/Cai)),1.4))*xs*xs*(V-Eks);
        
        //Inward rectifier K+ current, IK1 done
        IK1=x_K1_inf*gK1max*(V-Ek);
        
        //Na+/Ca2+ exchanger current, INaCa, not than n represents lambda done
        INaCa = knaca*(1/(KmNai*KmNai*KmNai+Nao_naca*Nao_naca*Nao_naca))*(1/(KmCa+Cao))*(1/(1+ksat*exp((n-1)*V*F/(R*T))))*(exp(n*V*F/(R*T))*Nai*Nai*Nai*Cao-exp((n-1)*V*F/(R*T))*Nao_naca*Nao_naca*Nao_naca*Cai*alfa);
        
        //Na+/K+ pump current, INaK, note that code uses knak variable for pnak
        //INaK = knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*(1./(1.+0.1245*exp(-0.1*V*F/(R*T))+0.0353*exp(-V*F/(R*T))));
        INaK = (1/(1+0.1245*exp(-0.1*V*F/(R*T))+0.0353*exp(-V*F/(R*T))))*(knak*Ko/(Ko+KmK)*Nai/(Nai+KmNa));
        
        //Calcium pump current, IpCa done
        IpCa = GpCa*Cai/(kpca+Cai);
        
        //Plateau K+ current done
        IpK = GpK*(V-Ek)*(1/(1+exp((25-V)/5.98)));
        
        //Background sodium current done
        IbNa = GbNa*(V-Ena);
        
        //Calcium background current, IbCa done
        IbCa = GbCa*(V-Eca);
        
        //Currents added for the hESC_CM model
        ////Hyperpolarization activated funny current, If note this is a new current incorporated into the model
        //not initially in the ten Tusscher model done
        If = Gf*xf*(V+17);
        
        //T-type Ca2+ Current, ICaT, not originally in the tentusscher model done
        ICaT = dCaT*fCaT*GICaT*(V-Eca);
        
        dvdtold = dvdt; 
        
        Itotal = INa + ICaL + Ito + IKr + IKs + IK1 + INaCa + INaK + IpCa + IbNa + IpK + IbCa + If + ICaT;
        
        
        Ileak = Vleak*(CaSR-Cai);
        Iup = Vmaxup/(1.+((Kup*Kup)/(Cai*Cai)));
        
        //modified by RaIrel
        Irel = (d*g*(crel+arel*(CaSR*CaSR)/((brel*brel)+(CaSR*CaSR))))*RaIrel;
//Analytic solution of calcium dynamics taken from Pete Jordan's version of the original Ten
//Tusscher model        
        CaBuf = Bufc*Cai/(Cai+Kbufc);
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

        file_output_counter++;
        dvdt = -Itotal;
        
        vnew = V+DT*dvdt;
//Method of calculation of parameters taken from Kharche from modelDB         
        if(dvdt>=0.0 && dvdtold<0.0){
            vmin[counter] = V;
            tvmin[counter] = t;
            start_output = 1;
        }
        if(dvdt>dvdtmax[counter]&&start_output>0){
            dvdtmax[counter] = dvdt;
            apd_start[counter] = t;
        }
        if(dvdtold>0.0 && dvdt<=0.0){
            vmax[counter] = V;
            APA[counter] = vmax[counter] - vmin[counter];
            top_slope[counter] = (vmax[counter]-vmin[counter])/(t - tvmin[counter]);
        }
        if((counter>0)&&(dvdtold <= top_slope[counter-1]) && (dvdt>top_slope[counter-1])){
            top[counter] = V;
            ttop[counter] = t;
            ddr[counter] = (V - vmin[counter])/(t - tvmin[counter]);
        }
        if((vnew <= (top[counter] + 0.7*APA[counter])) && (V > (top[counter] + 0.7*APA[counter])) ){
            if(apd_start[counter]>0.0)
                apd30[counter] = t - apd_start[counter];
        }
        if((vnew <= (top[counter] + 0.5*APA[counter])) && (V > (top[counter] + 0.5*APA[counter])) ){
            if(apd_start[counter]>0.0)
                apd50[counter] = t - apd_start[counter];
        }
        if((vnew <= (top[counter] + 0.3*APA[counter])) && (V > (top[counter] + 0.3*APA[counter])) ){
          //cout << vnew << " " <<  (top[counter] + 0.3*APA[counter])  << " "<< V << endl;
            if(apd_start[counter]>0.0)
                apd70[counter] = t - apd_start[counter];
        }
        if((vnew <= (0.9*vmin[counter])) && (V > (0.9*vmin[counter]))){//Algorith fails for the case of APD90
            if(apd_start[counter]>0.0){
                apd_end[counter] = t;
                apd90[counter] = t - apd_start[counter];
            }
        }
        if (start_output && (apd_end[counter] > 0.0)){
                bcl[counter] = apd_start[counter]-apd_start[counter-1];
                printf("%10.5f\t",vmin[counter]);
                printf("%10.5f\t",vmax[counter]);
                printf("%10.5f\t",dvdtmax[counter]);
                printf("%10.5f\t",apd30[counter]);
                printf("%10.5f\t",apd50[counter]);
                printf("%10.5f\t",apd70[counter]);
                printf("%10.5f\t",apd90[counter]);
                printf("%10.5f\t",(60000/bcl[counter]));//convert cycle length to beats per minute
                printf("%10.5f\t",ddr[counter]);
                printf("%10.5f\t",top[counter]);//take off potential
                printf("%10.5f\t",APA[counter]);
                printf("\n");
                counter++;
        }
        V = vnew;
        if((file_output_counter%1000==0) && (t>200000))
        {
        //  if(print_header){
        //    print_header = false;
         // }
            //fprintf( output, "%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\n", time+(number-1)*basic_bcl, V, m, h, j, d, f, g, f_ca, r, s, xs, xr1, xr2 );
            //fprintf( output, "%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.105f\t%4.10f\t%4.10f\t%4.10f\n", time+(number-1)*basic_bcl, Ileak, Iup, Irel, CaBuf, CaSR, Cai, Nai, Ki, V );
            //fprintf( output, "%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\t%4.10f\n", t, IKr, IKs, IK1, Ito, INa, If, INaK, ICaL, IbCa, INaCa, V, ICaT, Cai, IpCa, Irel, Itotal,dvdt );
            fprintf( output, "%4.10f\t%4.10f\t%4.10f\n", t/1000,V/1000, dvdt);
        }
    }
    printf( "\nV = %f;\nCai = %f;\nCaSR = %f;\nm = %f;\nh = %f;\nj = %f;\nxr1 = %f;\nxr2 = %f;\nxs = %f;\nr = %f;\ns = %f;\nd = %f;\nf = %f;\nf_ca = %f;\ng = %f;\ndCaT = %f;\nfCaT = %f;\nxf = %f;\n", V, Cai, CaSR, m, h, j, xr1, xr2, xs, r, s, d, f, f_ca,g, dCaT, fCaT, xf );

    fclose( output );
    return 0;
}
