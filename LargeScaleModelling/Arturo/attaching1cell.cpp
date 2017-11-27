// ------------------------------------------------ Master Project ------------------------------------------ 
// -------------------------------------------- Attachment 1 cell----------------------------------------
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


void attachment_single_cell(vector<long double> &Bac_m,vector<long double> &Bac_e_d,vector<int> &Bac_s,vector<int> &Bac_c, vector<long double> &Bac_r,double Rmax);

// create a rng in order to create random numbers
gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);

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

const long double Lx = Nx*dx;  //% m, space length in x; defining a grid/matrix for placing the bacteria before sending to COMSOL
const long double Ly = Ny*dy;  //% m, space height in y


const float bac_h    = dx;
const int bac_rho = 100;
const int bac_rhod = 105;
const int bac_rhoe = 25;


//       1.7) Initializating vectors with cells information

long double bac_x[] = {0, -dx, dx, dx, -dx};  /* x's position of each cell*/
long double bac_y[] = {0, dy, -dy, dy, -dy};  /* y's position of each cell*/
long double bac_m[] = {mass_het,mass_aob,mass_nob,0,0};  /* mass of each cell*/
long double bac_e_d[] = {mass_eps,0,0,mass_eps,mass_inert}; /* mass of the eps and solid debries*/
int bac_s[] = {0,1,2,3,4};		  /* type of cell 1-HET 2-AOB 3-NOB 4-EPS 5-DEAD*/	
int bac_c[] = {0,1,2,0,0};        /* number of the cluster that the cells belongs to*/



vector<long double> Bac_x (bac_x, bac_x + sizeof(bac_x) / sizeof(bac_x[0]) );
vector<long double> Bac_y (bac_y, bac_y + sizeof(bac_y) / sizeof(bac_y[0]) );
vector<long double> Bac_m (bac_m, bac_m + sizeof(bac_m) / sizeof(bac_m[0]) );
vector<long double> Bac_e_d (bac_e_d, bac_e_d + sizeof(bac_e_d) / sizeof(bac_e_d[0]) );
vector<int> Bac_s (bac_s, bac_s + sizeof(bac_s) / sizeof(bac_s[0]) );
vector<int> Bac_c (bac_c, bac_c + sizeof(bac_c) / sizeof(bac_c[0]) );
vector<long double> Bac_r, Bac_ra;


vector<double> Bac_theta,Bac_R;



//       1.6) Parameters single attach process

				   /* pHET pAOB pNOB */ 
const float prob[] = {0.6, 0.2, 0.2};  /* probabilities for the new cell in single attach process*/
const long double mass_new_cell[] = {mass_het,mass_aob,mass_nob}; 
const long double eps_new_cell[]={mass_eps,0,0};	




int main()
{
	// initialise the rng
	srand(time(NULL));	
	
	Bac_r.clear();
	Bac_ra.clear();
	for (unsigned int i=0; i<Bac_x.size();++i){
		Bac_r.push_back(sqrt((Bac_m[i]/bac_rho + Bac_e_d[i]/bac_rhoe)/PI/bac_h));
		Bac_ra.push_back(sqrt((Bac_m[i]/bac_rho)/PI/bac_h));

	}


	Bac_theta.clear();
	Bac_R.clear();
	
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

	
	for(int i=0;i<100;++i){
			
		/* calculating Rmax*/	
		vector<double>::const_iterator it, RmaxVector;
		long double Rmax; /*ver si no cambia con double */
		RmaxVector = max_element(Bac_R.begin(),Bac_R.end());
		Rmax = *RmaxVector;

		attachment_single_cell(Bac_m,Bac_e_d,Bac_s,Bac_c,Bac_r,Rmax);
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
	return 0;
}



void attachment_single_cell(vector<long double> &Bac_m,vector<long double> &Bac_e_d,vector<int> &Bac_s,vector<int> &Bac_c, vector<long double> &Bac_r,double Rmax)
{
	// -------INPUTS----------
	/* "Rmax" is the radio of the farthest cell from the center*/ /* REMOVE THE OTHERS INPUTS CAREFULLY*/
	
	// -------OUTPUTS---------
	/* "Bac_x", "Bac_y" ... are modified inside the function*/
	
	// -------EXPLANATION----------	
	/* A new cell (random type) is attached to the floc with a random radio */
	
		
	/* generate a random position for the new cell*/
	
		gsl_rng_set (r, rand() + 1);
		// generates a uniform in 0,2PI
		double fi_atta =  gsl_rng_uniform (r)*2*PI;
		long double R_atta = Rmax + 5*Bac_r[0]; 
		

	/* attracting the new cell to the center of the floc*/
	int j=0;
	while (j!=1){
		R_atta = R_atta - Bac_r[0];
		if(R_atta<0){
			Bac_x.push_back(0);
			Bac_y.push_back(0);
			j=1;
			break;
		}
		long double x_atta = R_atta * cos(fi_atta);
		long double y_atta = R_atta * sin(fi_atta);

		for (unsigned int i=0; i<Bac_x.size();++i){					
			long double dist = sqrt(pow((Bac_x[i] - x_atta),2) + pow((Bac_y[i] - y_atta),2));
			if ((dist - 3*Bac_r[i])<0){ /* means new cell is close enough to that cell */
				Bac_x.push_back(x_atta);
				Bac_y.push_back(y_atta);
				j=1;	                  /* breaks the initial for*/
				break;                    /* breaks this for */
			}
		}				 
	}

	/* choosing type of the new cell*/
	
	
	gsl_rng_set (r, rand() + 1);
	// generates a uniform in 0,1
	double ra  = gsl_rng_uniform (r); 

	for(int j=1; j<4; ++j){
		valarray<float> prob1 (prob,j);
		if (ra< prob1.sum()){			
			/* adding the new cells with its properties to the floc*/	
		
			Bac_m.push_back(mass_new_cell[(j-1)]);
			Bac_e_d.push_back(eps_new_cell[(j-1)]);
			Bac_s.push_back(j-1);

			if (j>1){
				vector<int>::const_iterator it, clusterMax;
				clusterMax = max_element(Bac_c.begin(),Bac_c.end());
				Bac_c.push_back(*clusterMax+1);	
			}else{
				Bac_c.push_back(0);	
			}
			break;		
		}			
	} 
	
	//Bac_inf.push_back("att_single_cell");
	Bac_r.push_back(sqrt((Bac_m.back()/bac_rho + Bac_e_d.back()/bac_rhoe)/PI/bac_h));
	Bac_ra.push_back(sqrt((Bac_m.back()/bac_rho)/PI/bac_h));
	
	Bac_theta.push_back(atan2(Bac_y.back(),Bac_x.back()));
	Bac_R.push_back (sqrt(pow(Bac_x.back(),2) + pow(Bac_y.back(),2)));
	
	++number_bacterium;
}
