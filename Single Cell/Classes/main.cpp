#include <stdio.h>
#include <math.h>
#include <iostream>
#include "hESC-CM.h"
//#include "LivRudy.h"
#include "Environment_hESCM.h"
using namespace std;

//APD paramters
double dvdt,dvdtold;
double vnew;

int main(){
    FILE *output;
    char output_file_name[ 20 ];
    output = fopen( "teste.dat", "w" );
    int file_output_counter=0;

    hESCM hescm = hESCM('a');
    double tmin = 0.0;
    double tmax = 100000;
    double bcl = 1000;
    double stim_dur = 1;
    int bcl_int = bcl/DT;
    int stim_int = stim_dur/DT;
    int n_stim = 0;
    double Istim;
    double t;
    bool stim_active = false;
    t=tmin;
    int count = 0;
    for (t = tmin + DT; t<=tmax; t+=DT) {
      count++;
      Istim = 0.0; 
      if ((count%bcl_int) == 0){
       // cout << count << endl;
        stim_active = true;
      }
      if ((n_stim <= stim_int) && stim_active){
        Istim = -52; 
        n_stim++;
      }
      if (n_stim == stim_int) {
  
        stim_active = false; 
        n_stim = 0;
      }
    //  cout << Istim << endl;
        hescm.solve(Istim,t);
        file_output_counter++;
        if((file_output_counter%1000==0))
        {
          hescm.print_vars( false, output, t, 1);
        }
    }
    hescm.print_vars( true, output, t, 1);

fclose (output);
return 0;
}
