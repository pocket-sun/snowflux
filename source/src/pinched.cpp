#include <iostream>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string> 
#include <sstream>

using namespace std;

/*
 * functions from supernova_mixing.cpp
 */

// Function to write energy and fluxes into fluxfile 
static void write(double a, double B[], ofstream& outfile){

  // No flavor transitions assumed
  // Input:  B[0]: nue, B[1]: nuebar, B[2]: nux (any one of them)
  // Output flux file: first neutrinos, e, mu, tau, then antinus, ebar, mubar, taubar
  // This is basically writing out the flux before mixing

  outfile << setw(8) << a << "\t " ;
  outfile << setw(8) << B[0] << "\t " ;
  outfile << setw(8) << B[2] << "\t " ;
  outfile << setw(8) << B[2] << "\t " ;
  outfile << setw(8) << B[1] << "\t " ;
  outfile << setw(8) << B[2] << "\t " ;
  outfile << setw(8) << B[2] << "\t " ;

  outfile << "\n";

}


// Function to write energy and fluxes into fluxfile 
static void write_nh(double a, double B[], double th12, ofstream& outfile){

  // Input:  B[0]: nue, B[1]: nuebar, B[2]: nux (any one of them, all equal)
  // Output flux file: first neutrinos, e, mu, tau, then antinus, ebar, mubar, taubar

  // Normal ordering:
  // nue=nux0
  // nuebar = cos^2th12 nuebar0 + sin^2th12 nux0
  // Pee = p = 0
  // Peebar = pbar = cos^2th12 

  // Neutrinos:
  // numu+nutau = (1-p)nue0 + (1+p)nux0 
  //  so numu=nutau = (nue0+nux0)/2
   
  // Antineutrinos
  // numubar + nutaubar = ((1-pbar)nuebar0 + (1+pbar)nuxbar0)
  //  so numubar = nutaubar = (sin^2th12 nuebar0 + (1+cos^2th12) nuxbar0)/2
  
  //For  th12 = 0.588336
  //double s2th12 = 0.308;
  //double c2th12 = 0.692;

  double s2th12;
  double c2th12;
  s2th12 = pow(sin(th12),2);
  c2th12 = 1-s2th12;

  outfile << setw(8) << a << "\t " ;
  outfile << setw(8) << B[2] << "\t " ;
  outfile << setw(8) << (B[0]+B[2])/2. << "\t " ;
  outfile << setw(8) << (B[0]+B[2])/2. << "\t " ;
  outfile << setw(8) << c2th12*B[1]+s2th12*B[2] << "\t " ;
  outfile << setw(8) << ((1.-c2th12)*B[1]+(1+c2th12)*B[2])/2. << "\t " ;
  outfile << setw(8) << ((1.-c2th12)*B[1]+(1+c2th12)*B[2])/2. << "\t " ;
  outfile << "\n";

}


// Function to write energy and fluxes into fluxfile 
static void write_ih(double a, double B[], double th12, ofstream& outfile){


  // Input:  B[0]: nue, B[1]: nuebar, B[2]: nux (any one of them, all equal)
  // Output flux file: first neutrinos, e, mu, tau, then antinus, ebar, mubar, taubar

  // Inverted ordering:
  // nuebar=nuxbar0
  // nue = sin^2th12 nue0 + cos^2th12 nux0
  // Pee = p = sin^2th12
  // Peebar = pbar = 0

  // Neutrinos
  // numu + nutau = ((1-p)nue0 + (1+p)nux0)
  //  so numu = nutau = (cos^2th12 nue0 + (1+sin^2th12) nux0)/2

  // Antineutrinos:
  // numubar+nutaubar = (1-pbar)nuebar0 + (1+pbar)nuxbar0 
  //  so numubar=nutaubar = (nuebar0+nuxbar0)/2
   
  
  // th12 = 0.588336

  //  double s2th12 = 0.308;
  //double c2th12 = 0.692;

  double s2th12;
  double c2th12;
  s2th12 = pow(sin(th12),2);
  c2th12 = 1-s2th12;

  // First neutrinos then antinus
  outfile << setw(8) << a << "\t " ;
  outfile << setw(8) << s2th12*B[0]+c2th12*B[2] << "\t " ;
  outfile << setw(8) << ((1-s2th12)*B[0]+(1+s2th12)*B[2])/2. << "\t " ;
  outfile << setw(8) << ((1-s2th12)*B[0]+(1+s2th12)*B[2])/2. << "\t " ;
  outfile << setw(8) << B[2] << "\t " ;
  outfile << setw(8) << (B[1]+B[2])/2. << "\t " ;
  outfile << setw(8) << (B[1]+B[2])/2. << "\t " ;
  outfile << "\n";

}

//-----------------------------------------------------------------------------

// Function to calculate the energy spectrum
static double phi(double E_nu, double E_nu0, double alpha) {
    double N=pow((alpha+1.),(alpha+1.))/(E_nu0*tgamma(alpha+1.));
	double R=N*pow((E_nu/E_nu0),alpha)*exp((-1.)*(alpha+1.)*E_nu/E_nu0); 

	return R;
}


//-----------------------------------------------------------------------------

/*
alpha_nu_e  alpha_nubar_e  alpha_nu_x  <E>_nu_e  <E>_nubar_e  <E>_nu_x  L_nu_e  L_nubar_e  L_nu_x
units: 1, MeV, erg
*/
int GarchingFluence(double alpha[], double E0_0[], double L_0[], double dist){

    const double kpc2cm=3.08568025e21; // [dist]=cm, 1kpc
    const double gevpererg = 624.15; // 1 erg = 624.151 GeV
    double estep=0.0002; // 0.0002 GeV
    double th12=0.; // no oscillation
    dist = dist*kpc2cm;

    ifstream infile;
    ofstream outfile;
    ofstream outfile_nh;
    ofstream outfile_ih;
    double E_nu;
    double F[3];
    double L[3], E0[3];

	// E must be in GeV
	// Convert luminosity to GeV/s
	for(int k=0; k<3; k++){
	    E0[k] = E0_0[k]/1000.;
        L[k] = L_0[k]*gevpererg;
	}

	// create filename for output file and open respective file 
	//		string filename="pinched_";
	outfile.open("./snowglobe/fluxes/pinched_0.dat"); // out and trunc implicit

	if (!outfile.good()) {
        cerr << "fail to open pinched_0.dat" << endl;
	    exit(-1);
	}

    
	if (fabs(th12)>0) {
	  // Make the MSW files

	    string filename_nh;
	    outfile_nh.open(filename_nh.c_str());

	    if (!outfile_nh.good()) {
	      cout << "Outfile "<<filename_nh<<" not opened..."<<endl;
	      cout << "Check that directory pointed to by OUTFLUXDIR environment variable exists"<<endl;
	      exit(-1);
	    }

	    string filename_ih;
	    outfile_ih.open(filename_ih.c_str());

	    if (!outfile_ih.good()) {
	      cout << "Outfile "<<filename_ih<<" not opened..."<<endl;
	      cout << "Check that directory pointed to by OUTFLUXDIR environment variable exists"<<endl;
	      exit(-1);
	    }
	}

	// File should have flux in the 0.0002 GeV bin, so multiply by bin size
	E_nu=0;
	for(int i=0; i<=500; i++) {
	  // Create fluxes F^0
	  for(int j=0; j<3; j++){

	    if (E0[j]>0.) {
	      F[j]=1/(4*M_PI*dist*dist)*L[j]/E0[j]*phi(E_nu,E0[j],alpha[j])*estep;  
	    } else {
	      F[j]=0.;
	    }
	    }
	  // Write data to output file 
	  write(E_nu,F, outfile); // energies in output file need to be in GeV!
	  if (fabs(th12)>0) {
	    write_nh(E_nu,F,th12, outfile_nh);
	    write_ih(E_nu,F,th12, outfile_ih);
	  }

	  E_nu+=estep;
	}	
	outfile.close();

	if (fabs(th12)>0) {
	  outfile_nh.close();
	  outfile_ih.close();
	}
	return 0; 
}

