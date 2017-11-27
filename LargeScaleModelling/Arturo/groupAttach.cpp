// ------------------------------------------------ Master Project ------------------------------------------ 
// ---------------------------------------- Attachment group of cells ----------------------------------------
// --------------------------------------------Arturo Alvarez-Arenas ----------------------------------------
// ------------------------------------------- Newcastle University -----------------------------------------

#include <numeric> // accumulate
#include <iostream>     // std::cout
#include <cstddef>      // std::size_t
#include <valarray>     // std::valarray
#include <ctime>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <numeric>
#include <boost/range/numeric.hpp>
#include <cmath>        // std::abs
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define PI 3.1415926535897932384626433832795

// create a rng in order to create random numbers
gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);

using namespace std;
/* group attach*/
/* we should move all the bacatt_x and bacatt_y in some moment of the code to the (0,0) doing -Lx/2 and Ly/2.
 In this code it is assumed that that stuff have already been done */ 

//void attachment_group_cell(vector<long double> &Bac_x,vector<long double> &Bac_y,vector<long double> &Bac_m,vector<long double> &Bac_e_d,vector<int> &Bac_s,vector<int> &Bac_c,vector<long double> &Bac_r,vector<long double> &Bac_ra,vector<double> &Bac_theta, vector<double> &Bac_R,double Rmax,vector<long double> &Bacatt_x,vector<long double> &Bacatt_y,vector<long double> &Bacatt_m,vector<long double> &Bacatt_e_d, vector<int> &Bacatt_s,vector<int> &Bacatt_c,vector<long double> &Bacatt_r,vector<long double> &Bacatt_ra,vector<double> &Bacatt_theta,vector<double> &Bacatt_R);
void attachment_group_cell(double Rmax,vector<long double> &Bacatt_x ,vector<long double> &Bacatt_y,vector<long double> &Bacatt_m,vector<long double> &Bacatt_e_d,vector<int> &Bacatt_s,vector<int> &Bacatt_c,vector<long double> &Bacatt_r,vector<long double> &Bacatt_ra,vector<double> &Bacatt_theta,vector<double> &Bacatt_R);


// constants


const long double mass_het = 1e-16;
const long double mass_aob = 1e-16;
const long double mass_nob = 1e-16;
const long double mass_eps = 2.6e-17;  // % kg for grid 1e-6 % for having eps density 25 kg/m3
const long double mass_inert = 1.1e-16; //% kg for grid 1e-6 % for having eps density 105 kg/m3

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

long double bac_x[] = {0, -dx, dx, dx, -dx,1.47887e-06};  /* x's position of each cell*/
long double bac_y[] = {0, dy, -dy, dy, -dy,-2.32427e-06};  /* y's position of each cell*/
long double bac_m[] = {mass_het,mass_aob,mass_nob,0,0,mass_het};  /* mass of each cell*/
long double bac_e_d[] = {mass_eps,0,0,mass_eps,mass_inert,mass_eps}; /* mass of the eps and solid debries*/
int bac_s[] = {0,1,2,3,4,0};		  /* type of cell 1-HET 2-AOB 3-NOB 4-EPS 5-DEAD*/	
int bac_c[] = {0,1,2,0,0,0};        /* number of the cluster that the cells belongs to*/



vector<long double> Bac_x (bac_x, bac_x + sizeof(bac_x) / sizeof(bac_x[0]) );
vector<long double> Bac_y (bac_y, bac_y + sizeof(bac_y) / sizeof(bac_y[0]) );
vector<long double> Bac_m (bac_m, bac_m + sizeof(bac_m) / sizeof(bac_m[0]) );
vector<long double> Bac_e_d (bac_e_d, bac_e_d + sizeof(bac_e_d) / sizeof(bac_e_d[0]) );
vector<int> Bac_s (bac_s, bac_s + sizeof(bac_s) / sizeof(bac_s[0]) );
vector<int> Bac_c (bac_c, bac_c + sizeof(bac_c) / sizeof(bac_c[0]) );
vector<long double> Bac_r;
vector<long double> Bac_ra;

vector<double> Bac_theta;
vector<double> Bac_R;


/* los for tienen que estar dentro del main, y haz lo de la pension!!! */


// Inventarse vectores para el group attach.


long double bacatt_x[] = {0+Lx/2, -dx+Lx/2, dx+Lx/2, dx+Lx/2, -dx+Lx/2};  /* x's position of each cell*/
long double bacatt_y[] = {0+Ly/2, dy+Ly/2, -dy+Ly/2, dy+Ly/2, -dy+Ly/2};  /* y's position of each cell*/
long double bacatt_m[] = {mass_het,mass_het,mass_nob,mass_aob,mass_nob};  /* mass of each cell*/
long double bacatt_e_d[] = {mass_eps,mass_eps,0,0,0}; /* mass of the eps and solid debries*/
int bacatt_s[] = {1,1,3,2,3};		  /* type of cell 1-HET 2-AOB 3-NOB 4-EPS 5-DEAD*/	
int bacatt_c[] = {0,0,3,4,3};        /* number of the cluster that the cells belongs to*/


vector<long double> Bacatt_x (bacatt_x, bacatt_x + sizeof(bacatt_x) / sizeof(bacatt_x[0]) );
vector<long double> Bacatt_y (bacatt_y, bacatt_y + sizeof(bacatt_y) / sizeof(bacatt_y[0]) );
vector<long double> Bacatt_m (bacatt_m, bacatt_m + sizeof(bacatt_m) / sizeof(bacatt_m[0]) );
vector<long double> Bacatt_e_d (bacatt_e_d, bacatt_e_d + sizeof(bacatt_e_d) / sizeof(bacatt_e_d[0]) );
vector<int> Bacatt_s (bacatt_s, bacatt_s + sizeof(bacatt_s) / sizeof(bacatt_s[0]) );
vector<int> Bacatt_c (bacatt_c, bacatt_c + sizeof(bacatt_c) / sizeof(bacatt_c[0]) );
vector<long double> Bacatt_r;
vector<long double> Bacatt_ra;

vector<double> Bacatt_theta;
vector<double> Bacatt_R;

//Bacatt_x,Bacatt_y,Bacatt_m,Bacatt_e_d,Bacatt_s,Bacatt_c,Bacatt_r,Bacatt_ra,Bacatt_theta,Bacatt_R
//vector<long double> &Bacatt_x,vector<long double> &Bacatt_y,vector<long double> &Bacatt_m,vector<long double> &Bacatt_e_d, vector<int> &Bacatt_s,vector<int> &Bacatt_c,vector<long double> &Bacatt_r,vector<long double> &Bacatt_ra,vector<double> &Bacatt_theta,vector<double> &Bacatt_R

int main(){
	// initialise the rng
	srand(time(NULL));	
	
	for (unsigned int i=0; i<Bac_x.size();++i){
		Bac_r.push_back(sqrt((Bac_m[i]/bac_rho + Bac_e_d[i]/bac_rhoe)/PI/bac_h)); /* it works just for the first iteration when bac_e_d(AOB&NOB) = 0 */
		Bac_ra.push_back(sqrt((Bac_m[i]/bac_rho)/PI/bac_h));
	}

	for (unsigned int i=0; i<Bac_x.size();++i){				/* changing from cartesian to polar coordinates*/	
		Bac_theta.push_back(atan2(Bac_y[i],Bac_x[i]));
		Bac_R.push_back (sqrt(pow(Bac_x[i],2) + pow(Bac_y[i],2)));
	}
	

	
	for (unsigned int i=0; i<Bacatt_x.size();++i){
		Bacatt_r.push_back(sqrt((Bacatt_m[i]/bac_rho + Bacatt_e_d[i]/bac_rhoe)/PI/bac_h));
		Bacatt_ra.push_back(sqrt((Bacatt_m[i]/bac_rho)/PI/bac_h));
	}



	for (unsigned int i=0; i<Bacatt_x.size();++i){				/* changing from cartesian to polar coordinates*/	
		Bacatt_theta.push_back(atan2(Bacatt_y[i],Bacatt_x[i]));
		Bacatt_R.push_back (sqrt(pow(Bacatt_x[i],2) + pow(Bacatt_y[i],2)));
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
	for (int i=0;i<1;++i){
		cout <<"\n\n i: " << i;
		cout <<"\n\n i: " << i;
			/* calculating Rmax*/	
		vector<double>::const_iterator it, RmaxVector;
		double Rmax;
		RmaxVector = max_element(Bac_R.begin(),Bac_R.end());
		Rmax = *RmaxVector;
		
		attachment_group_cell(Rmax,Bacatt_x ,Bacatt_y,Bacatt_m,Bacatt_e_d,Bacatt_s,Bacatt_c,Bacatt_r,Bacatt_ra,Bacatt_theta,Bacatt_R);
	}
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


void attachment_group_cell(double Rmax,vector<long double> &Bacatt_x ,vector<long double> &Bacatt_y,vector<long double> &Bacatt_m,vector<long double> &Bacatt_e_d,vector<int> &Bacatt_s,vector<int> &Bacatt_c,vector<long double> &Bacatt_r,vector<long double> &Bacatt_ra,vector<double> &Bacatt_theta,vector<double> &Bacatt_R){	
	// -------INPUTS----------
	/* "Rmax" is the radio of the farthest cell from the center*/
	
	// -------OUTPUTS---------
	/* "Bac_x", "Bac_y" ... are modified inside the function*/
	
	// -------EXPLANATION----------	
	/* A new group of cells is attached to the floc. A random angle is selected and the farthest cell close to that radio is choosen.
	   The group of cells is attached to that cell */
	
		
	/* Calculating the center of the new mini floc */
	long double m_tot_att = 0;
	long double m_x_tot_att = 0;
	long double m_y_tot_att = 0;

	for (unsigned int i=0; i<Bacatt_m.size();++i){
		m_tot_att += Bacatt_m[i];                       /* mass total of the mini floc*/
		m_x_tot_att += Bacatt_m[i]*Bacatt_x[i]; 		
		m_y_tot_att += Bacatt_m[i]*Bacatt_y[i];	
	}

	long double xcatt,ycatt;
	xcatt=m_x_tot_att/m_tot_att; /* x coordinate of the minifloc's center */
	ycatt=m_y_tot_att/m_tot_att; /* y coordinate of the minifloc's center */
	
	
	
	long double Bacgroup_R = (sqrt(pow(xcatt,2) + pow(ycatt,2))); /* Radio of the minifloc's center*/

	
	/* Selecting the closest cell with a random angle */
	
	int index = -1;
	long double reffer_R = 0;
	while(index<0){
		

		gsl_rng_set (r, rand() + 1);
		// generates a uniform in [-PI PI]	
		double fi_atta =  gsl_rng_uniform (r) * 2 * PI - PI; /* random angle between [-PI PI]*/
		
		for (unsigned int i=0; i<Bac_theta.size();++i){
			if(abs(Bac_theta[i]-fi_atta)<=(PI/20) && Bac_R[i]>=reffer_R){ /* if the cell is in that angle and the radium is larger than the previous one is selected*/
				reffer_R = Bac_R[i];
				index = i;			
			}
		}
	}


	/* coordinates of the closest cell*/
	long double reffer_theta = Bac_theta[index];
	long double x_reffer = Bac_x[index];
	long double y_reffer = Bac_y[index];
	
	/* Selecting the angle and the initial radio of attachment */
	
	long double Bacgroup_theta = reffer_theta;    /* the angle of the attachment is set as the angle of the initial cell */
	Bacgroup_R += Rmax;							/* the radio of the minifloc's center is recalculated*/	
	
	/* Recalculating new x's and y's position of the floc's center */
	
	long double xcatt_new = Bacgroup_R*cos(Bacgroup_theta);
	long double ycatt_new = Bacgroup_R*sin(Bacgroup_theta);
	
	/* Checking how much the center has been moved and move all the cells*/

	
	long double dxatt = xcatt_new - xcatt;
	long double dyatt = ycatt_new - ycatt;
	
	for (unsigned int i=0; i<Bacatt_x.size();++i){
		Bacatt_x[i] += dxatt;
		Bacatt_y[i] += dyatt;			
	}
	

	
	// putting a generic radio in case the radio is 0
	double dist_move;
	if (Bac_r[index] !=0){
		dist_move = Bac_r[index];
	}else{
		dist_move = 0.538e-06;
	}
	
	for (int j=-1;j<0;--j){

		long double xcatt_trans = Bacgroup_R*cos(Bacgroup_theta);
		long double ycatt_trans = Bacgroup_R*sin(Bacgroup_theta);
		dxatt = xcatt_trans - xcatt_new;
		dyatt = ycatt_trans - ycatt_new;
		
		xcatt_new = xcatt_trans;
		ycatt_new = ycatt_trans;
		for (unsigned int i=0; i<Bacatt_x.size();++i){
			Bacatt_x[i] += dxatt;
			Bacatt_y[i] += dyatt;	
		}		
		
		for (unsigned int i=0; i<Bacatt_x.size();++i){	
			long double dist = sqrt(pow((Bacatt_x[i] - x_reffer),2) + pow((Bacatt_y[i] - y_reffer),2));
			if ((dist - 3*dist_move)<0){ /* means new cell is close enough to that cell */
				j=1;	                  /* breaks the initial while*/
				break;                    /* breaks this for */
			}
			
		}
		
		Bacgroup_R -= dist_move;
	}
	
	/* Adding the new cells to the vectors */
	
	for (unsigned int i=0; i<Bacatt_x.size();++i){
		Bac_x.push_back(Bacatt_x[i]);
		Bac_y.push_back(Bacatt_y[i]);
		Bac_m.push_back(Bacatt_m[i]);
		Bac_e_d.push_back(Bacatt_e_d[i]);
		Bac_s.push_back(Bacatt_s[i]);
		Bac_c.push_back(Bacatt_c[i]);
		Bac_r.push_back(Bacatt_r[i]);
		Bac_ra.push_back(Bacatt_ra[i]);
		Bac_theta.push_back(Bacatt_theta[i]);
		Bac_R.push_back(Bacatt_R[i]);
		//Bac_inf.push_back("att_group_cell");
	}	
}











