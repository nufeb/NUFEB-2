// ------------------------------------------------ Master Project ------------------------------------------ 
// -------------------------------------------- Division process ----------------------------------------
// --------------------------------------------Arturo Alvarez-Arenas ----------------------------------------
// ------------------------------------------- Newcastle University -----------------------------------------

#include <iostream>     // std::cout
#include <cstddef>      // std::size_t
#include <valarray>     // std::valarray
#include <ctime>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <stdio.h>      /* printf */
#include <math.h>       /* pow */
#include <algorithm>    // std::min_element, std::max_element
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define PI 3.1415926535897932384626433832795

using namespace std;


void exc_eps_division();
void EPS_excretion(int i);
void Division(int i);

const long double mass_het = 1e-16;
const long double mass_aob = 1e-16;
const long double mass_nob = 1e-16;
const long double mass_eps = 2.6e-17;  // % kg for grid 1e-6 % for having eps density 25 kg/m3
const long double mass_inert = 1.1e-16; //% kg for grid 1e-6 % for having eps density 105 kg/m3

const long double mass_cells[] = {mass_het,mass_aob,mass_nob,mass_eps,mass_inert}; 

long number_bacterium = 5; /* Total number of bacterium, it will be modified during this process*/


const int Nx = 450;       //grid nodes in x direction
const int Ny = 450;       // grid nodes in y direction

const long double dx = 1.1e-6;

const long double dy = 1.1e-6;
const long double dt = 0.1;

const long double Lx = Nx*dx;  //% m, space length in x; defining a grid/matrix for placing the bacteria before sending to COMSOL
const long double Ly = Ny*dy;  //% m, space height in y


const float bac_h    = dx;
const int bac_rho = 100;
const int bac_rhod = 105;
const int bac_rhoe = 25;
// kinetic constants and others




//       1.7) Initializating vectors with cells information

long double bac_x[] = {0, -dx, dx, dx, -dx};  /* x's position of each cell*/
long double bac_y[] = {0, dy, -dy, dy, -dy};  /* y's position of each cell*/
long double bac_m[] = {2.1*mass_het,mass_aob,mass_nob,0,0};  /* mass of each cell*/
long double bac_e_d[] = {mass_eps,0,0,mass_eps,mass_inert}; /* mass of the eps and solid debries*/
int bac_s[] = {0,1,2,3,4};		  /* type of cell 0-HET 1-AOB 2-NOB 3-EPS 4-DEAD*/	
int bac_c[] = {0,1,2,0,0};        /* number of the cluster that the cells belongs to*/



vector<long double> Bac_x (bac_x, bac_x + sizeof(bac_x) / sizeof(bac_x[0]) );
vector<long double> Bac_y (bac_y, bac_y + sizeof(bac_y) / sizeof(bac_y[0]) );
vector<long double> Bac_m (bac_m, bac_m + sizeof(bac_m) / sizeof(bac_m[0]) );
vector<long double> Bac_e_d (bac_e_d, bac_e_d + sizeof(bac_e_d) / sizeof(bac_e_d[0]) );
vector<int> Bac_s (bac_s, bac_s + sizeof(bac_s) / sizeof(bac_s[0]) );
vector<int> Bac_c (bac_c, bac_c + sizeof(bac_c) / sizeof(bac_c[0]) );
vector<long double> Bac_r, Bac_ra;


vector<double> Bac_theta,Bac_R;


//       1.9) Parameters for division

const long double bac_mmax[] = {2*mass_het,2*mass_aob,2*mass_nob,100,100}; /* The last two elements are put in order to run the loop in the division processi without diferentiating types of cells */
const float DivisionFraction = 0.5;

int main(){
	
	
	for (unsigned int i=0; i<Bac_x.size();++i){
		Bac_r.push_back(sqrt((Bac_m[i]/bac_rho + Bac_e_d[i]/bac_rhoe)/PI/bac_h));
		Bac_ra.push_back(sqrt((Bac_m[i]/bac_rho)/PI/bac_h));
	}

	for (unsigned int i=0; i<Bac_x.size();++i){				/* changing from cartesian to polar coordinates*/	
		Bac_theta.push_back(atan2(Bac_y[i],Bac_x[i]));
		Bac_R.push_back (sqrt(pow(Bac_x[i],2) + pow(Bac_y[i],2)));
	}
	


	for(unsigned int i=0; i<Bac_m.size(); ++i){
		cout << "\n i: " << i;
		cout << "\n Bac_m: " << Bac_m[i];
		cout << "\n Bac_e_d: " << Bac_e_d[i];
		cout << "\n Bac_s: " << Bac_s[i];
		cout << "\n Bac_c: " << Bac_c[i];
		cout << "\n Bac_x: " << Bac_x[i];
		cout << "\n Bac_y: " << Bac_y[i];
		cout << "\n Bac_r: " << Bac_r[i];
		cout << "\n Bac_ra: " << Bac_ra[i];
		cout << "\n Bac_R: " << Bac_R[i];
		cout << "\n Bac_theta: " << Bac_theta[i];
	}
	
	exc_eps_division();
	cout << "\n\n\n\n after attachment group of cell:";
	for(unsigned int i=0; i<Bac_m.size(); ++i){		
		cout << "\n i: " << i;
		cout << "\n Bac_m: " << Bac_m[i];
		cout << "\n Bac_e_d: " << Bac_e_d[i];
		cout << "\n Bac_s: " << Bac_s[i];
		cout << "\n Bac_c: " << Bac_c[i];
		cout << "\n Bac_x: " << Bac_x[i];
		cout << "\n Bac_y: " << Bac_y[i];
		cout << "\n Bac_r: " << Bac_r[i];
		cout << "\n Bac_ra: " << Bac_ra[i];
		cout << "\n Bac_R: " << Bac_R[i];
		cout << "\n Bac_theta: " << Bac_theta[i];
		
	}	
	return(0);
}


void exc_eps_division(){
	// -------INPUTS----------
	
	// -------OUTPUTS---------
	/* "Bac_x", "Bac_y" ... are modified inside the function*/
	
	// -------EXPLANATION----------	
	/* the excretion of eps and division processes are carried out */
	
	const int n_old = Bac_x.size(); /* the number of cells is saved because might change during the loop*/
	for (int i=0;i<n_old;++i){  /* Each cell is checked */
		if (Bac_s[i] == 0 && Bac_r[i]>Bac_ra[i]*sqrt(2)){
			EPS_excretion(i);
		}
		if (Bac_m[i] > bac_mmax[Bac_s[i]]){
			Division(i);
		}
	}	
}

void EPS_excretion(int i){
	// -------INPUTS----------
	/* "index" is the index of the cell to remove*/
	
	// -------OUTPUTS---------
	/* "Bac_x", "Bac_y" ... are modified inside the function*/
	
	// -------EXPLANATION----------	
	/* this function eliminates all the properties of a cell */
		
	// create a rng
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
	// initialise the rng
	srand(time(NULL));
	gsl_rng_set (r, rand() + 1);
	// generates a uniform in 0,1
	float fi = gsl_rng_uniform (r)*PI;
	float frac_ra = 0.2*(gsl_rng_uniform (r)-0.5);
	long double fraction = DivisionFraction + frac_ra;
	
	/* Updating data new cell*/
	Bac_x.push_back(Bac_x[i] + 2*Bac_r[i]*cos(fi));
	Bac_y.push_back(Bac_y[i] + 2*Bac_r[i]*sin(fi));
	Bac_s.push_back(3);
	Bac_c.push_back(Bac_c[i]);	
	Bac_m.push_back(0);
	Bac_e_d.push_back(fraction*Bac_e_d[i]);
	Bac_ra.push_back(0);
	Bac_r.push_back(sqrt((Bac_e_d.back()/bac_rhoe)/PI/bac_h));
	Bac_theta.push_back(atan2(Bac_y.back(),Bac_x.back()));
	Bac_R.push_back (sqrt(pow(Bac_x.back(),2) + pow(Bac_y.back(),2)));	
	/* Recalculating some properties previous cell*/
	Bac_e_d[i] *= (1-fraction);
	Bac_r[i] = sqrt((Bac_m[i]/bac_rho + Bac_e_d[i]/bac_rhoe)/PI/bac_h);
	
	/* Updating indexs*/
	//Index_cells[3].push_back(Bac_x.size()-1);  --------- INCLUDE THIS IN THE WASTEWATER
}

void Division(int i){
	// -------INPUTS----------
	/* "index" is the index of the cell that will divide*/
	
	// -------OUTPUTS---------
	/* "Bac_x", "Bac_y" ... are modified inside the function*/
	
	// -------EXPLANATION----------	
	/* the division process is carried out */
	
		
	// create a rng
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
	// initialise the rng
	srand(time(NULL));
	gsl_rng_set (r, rand() + 1);
	// generates a uniform in 0,1
	float fi = gsl_rng_uniform (r)*PI;	
	float frac_ra = 0.2* (gsl_rng_uniform (r)-0.5);
	long double fraction1 = DivisionFraction + frac_ra;

	/* Updating data new cell*/
	Bac_x.push_back(Bac_x[i] + Bac_r[i]*cos(fi));
	Bac_y.push_back(Bac_y[i] + Bac_r[i]*sin(fi));
	Bac_s.push_back(Bac_s[i]);
	Bac_c.push_back(Bac_c[i]);
	Bac_m.push_back(fraction1*Bac_m[i]);
	Bac_e_d.push_back(fraction1*Bac_e_d[i]);
	Bac_ra.push_back(sqrt((Bac_m.back()/bac_rho)/PI/bac_h));
	Bac_r.push_back(sqrt((Bac_m.back()/bac_rho + Bac_e_d.back()/bac_rhoe)/PI/bac_h));	
	Bac_theta.push_back(atan2(Bac_y.back(),Bac_x.back()));
	Bac_R.push_back (sqrt(pow(Bac_x.back(),2) + pow(Bac_y.back(),2)));
	
	/* Moving the previous cell*/
	Bac_x[i] -= Bac_r[i]*cos(fi);
	Bac_y[i] -= Bac_r[i]*sin(fi);
	/* Recalculating mass and size*/
	Bac_m[i] *= (1-fraction1);
	Bac_e_d[i] *= (1-fraction1);
	Bac_ra[i] = sqrt((Bac_m[i]/bac_rho)/PI/bac_h);
	Bac_r[i] = sqrt((Bac_m[i]/bac_rho + Bac_e_d[i]/bac_rhoe)/PI/bac_h);
	
	/* Updating indexs*/
	//Index_cells[Bac_s[i]].push_back(Bac_x.size()-1);  --------- INCLUDE THIS IN THE WASTEWATER
	
		
}
