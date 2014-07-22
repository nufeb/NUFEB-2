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

#define PI 3.1415926535897932384626433832795


using namespace std;
/* group detach*/


void detachment_group_cell(double Rmax);
void eliminate_cell(int index); /* This function is used in "detachment_group_cell" */

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
vector<long double> Bac_r;
vector<long double> Bac_ra;

vector<double> Bac_theta;
vector<double> Bac_R;







int main(){
	
	
	for (unsigned int i=0; i<Bac_x.size();++i){
		Bac_r.push_back(sqrt((Bac_m[i]/bac_rho + Bac_e_d[i]/bac_rhoe)/PI/bac_h));
		Bac_ra.push_back(sqrt((Bac_m[i]/bac_rho)/PI/bac_h));
	}

	for (unsigned int i=0; i<Bac_x.size();++i){				/* changing from cartesian to polar coordinates*/	
		Bac_theta.push_back(atan2(Bac_y[i],Bac_x[i]));
		Bac_R.push_back (sqrt(pow(Bac_x[i],2) + pow(Bac_y[i],2)));
	}
	
	/* calculating Rmax*/	
	vector<double>::const_iterator it, RmaxVector;
	double Rmax;
	RmaxVector = max_element(Bac_R.begin(),Bac_R.end());
	Rmax = *RmaxVector;
	
	
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

	detachment_group_cell(Rmax);

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

void detachment_group_cell(double Rmax){
	// -------INPUTS----------
	/* "Rmax" is the radio of the farthest cell from the center*/
	
	// -------OUTPUTS---------
	/* "Bac_x", "Bac_y" ... are modified inside the function*/
	
	// -------EXPLANATION----------	
	/* A random angle and radio are selected. The cells out of that radio and the clusters to which it belongs will be removed, */
		
	long double largest_radio = 0; 
	vector<int> ind_theta; /* vector with the index of some cells */
	
	/* A random angle is selected and we check if there is any cell in that angle +-10% and store them in ind_theta*/
	while(ind_theta.size()==0){
			double fi_atta = ((double) rand()/(RAND_MAX))*2*PI - PI;
			for (unsigned int i=0; i<Bac_theta.size();++i){
					if (abs(Bac_theta[i]-fi_atta)<= (PI/10)){
						ind_theta.push_back(i);
						/* For all those cells we check which one is the farthest from the center*/
						if (Bac_R[i]>= largest_radio){
							largest_radio = Bac_R[i];							
						}						
					}				
			}
	}	
	/* Select a random number in [0.5*radio_largest , radio largest] */
	double random_radio = (((double) rand()/(RAND_MAX))*0.5 + 0.5)*largest_radio;
	
	/* All the cells saved before with a radio larger than random_radio are removed */
	
	vector<int> cluster_to_remove;
	const unsigned int ini=ind_theta.size()-1;
	for (int i=ini;i>=0;--i){ /* The loop starts in the last element. So the positions of the other cells are not changed after removing one cell */
		if (Bac_R[ind_theta[i]]>=random_radio){
			cluster_to_remove.push_back(Bac_c[ind_theta[i]]); /* the cluster in which the cell was is saved before removing it*/
			eliminate_cell(ind_theta[i]);
		}	
	}


	/* All the clusters which cells that have been removed are erased */
	unsigned int i=0;
	unsigned int j;
	
	vector<int>::iterator it; 

	while(i<cluster_to_remove.size()){	
		/* find the first element with the same cluster*/
		it = find(Bac_c.begin(),Bac_c.end(),cluster_to_remove[i]); /* "it" is a pointer */
		j = distance(Bac_c.begin(),it);		  /* "j" is the index of the first cell that is in the cluster "cluster_to_remove[i]" */		 
	
		if ( j == Bac_c.size()){ /* if the element it is not found the function "find" returns the size of the vector (Bac_c) */
			++i; /* after removing all the elements of that cluster, we go on to the next cluster*/
		}else{
			eliminate_cell(j);
		}		
	}
}

void eliminate_cell(int index){
	// -------INPUTS----------
	/* "index" is the index of the cell to remove*/
	
	// -------OUTPUTS---------
	/* "Bac_x", "Bac_y" ... are modified inside the function*/
	
	// -------EXPLANATION----------	
	/* this function eliminates all the properties of a cell */
	Bac_x.erase(Bac_x.begin()+index);
	Bac_y.erase(Bac_y.begin()+index);
	Bac_m.erase(Bac_m.begin()+index);
	Bac_e_d.erase(Bac_e_d.begin()+index);
	Bac_c.erase(Bac_c.begin()+index);
	Bac_s.erase(Bac_s.begin()+index);
	Bac_r.erase(Bac_r.begin()+index);
	Bac_ra.erase(Bac_ra.begin()+index);
	Bac_R.erase(Bac_R.begin()+index);
	Bac_theta.erase(Bac_theta.begin()+index);		
}
