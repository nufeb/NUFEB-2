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
#include <math.h>       /* floor */

#define PI 3.1415926535897932384626433832795


using namespace std;
/* group detach*/


void removing_islands();
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

long double bac_x[] = {0 +Lx/2, -dx+Lx/2, dx+Lx/2, dx+Lx/2, -dx+Lx/2, 10*dx+Lx/2};  /* x's position of each cell*/
long double bac_y[] = {0+Ly/2, dy+Ly/2, -dy+Ly/2, dy+Ly/2, -dy+Ly/2, 10*dy+Ly/2};  /* y's position of each cell*/
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
	
	removing_islands();
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


void removing_islands(){
	// -------INPUTS----------
	
	// -------OUTPUTS---------
	/* "Bac_x", "Bac_y" ... are modified inside the function*/
	
	// -------EXPLANATION----------	
	/* All the cells that are not in a defined "main island" will be removed from the floc */
		
	// Define a matrix (dimension Nx*Ny) with 0s if there is no cell, or 1 if there is any cell in the down-left square
	int land[Nx][Ny] = {{ 0 }};
	//cout << "\n antes for (unsigned int i=0; i<Bac_x.size(); ++i) ";
	for (unsigned int i=0; i<Bac_x.size(); ++i){
		int jj,ii; 
		jj = (floor(Bac_x[i]/dx)+1);
		ii = (floor(Bac_y[i]/dy)+1);
		land[ii][jj] = 1;	
	}
	
	// The center of the main island is marked with 3
	//cout << "\n antes for (int i=(Ny/2 - 5); i<=(Ny/2 + 5); ++i) ";
	for (int i=(Ny/2 - 5); i<=(Ny/2 + 5); ++i){
		for (int j=(Nx/2 - 5); j<=(Nx/2 + 5); ++j){
			if (land[i][j]==1){
				land[i][j] = 3;
			}	
		}
	}
	//For each cell, if the square 5x5 that surrounds it intersect the main island, that position is marked as a 3

	int changes = 1; /* necessary to initialize the following loop*/
	while (changes > 0){
		changes = 0;
		for (int i=1;i<(Ny-1);++i){
			for (int j=1;j<(Nx-1);++j){
				if (land[i][j] ==1 && (land[i+1][j]==3 || land[i-1][j]==3 || land[i+2][j]==3 || land[i-2][j]==3 || land[i][j-1]==3 || land[i][j+1]==3 || land[i][j-2]==3 || land[i][j+2]==3 || land[i+1][j+1]==3 || land[i-1][j+1]==3 || land[i+2][j+1]==3 || land[i-2][j+1]==3 || land[i+1][j-1]==3 || land[i-1][j-1]==3 || land[i+2][j-1]==3 || land[i-2][j-1]==3 || land[i+1][j+2]==3 || land[i-1][j+2]==3 || land[i+2][j+2]==3 || land[i-2][j+2]==3 || land[i+1][j-2]==3 || land[i-1][j-2]==3 || land[i+2][j-2]==3 || land[i-2][j-2]==3) ){
				  land[i][j]=3;
				  changes++;
				}
			}
		}
	}
	
	// When there are no changes, all the cells with a 1 in that point is eliminated
	for (int i=Bac_x.size()-1; i>=0; --i){  /* It is done from the end so now the positions are not changed after eliminating cells */
		int ii,jj;
		jj = (floor(Bac_x[i]/dx)+1);
		ii = (floor(Bac_y[i]/dy)+1);
		if (land[ii][jj] == 1){
			eliminate_cell(i);
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
