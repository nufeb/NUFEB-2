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

#define PI 3.1415926535897932384626433832795

using namespace std;


void growth_decay();
void calculate_rates_index(int i);
void decay_one_cell(int a,int ran);
void eliminate_cell(int index);

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

// kinetic parameters
long double mu_max_het = 6;     // 1/d miu_max HET  (Henze ASM2d and others - but consider that here there is no storage
long double mu_max_aob = 0.76;   // 1/d miu_max AOB (Rittmann and McCarty, 2001) at 20 degrees
long double mu_max_nob = 0.81;    // 1/d miu_max NOB (Rittmann and McCarty, 2001) at 20 degrees
long double K_o2_het  = 0.81e-3; //kg/m3
long double K_s_het   = 10e-3;   //kg/m3
long double K_no2_het = 0.3e-3;  //kg/m3 % Alpkvist et al., 2006
long double K_no3_het = 0.3e-3;  //kg/m3
long double K_o2_aob  = 0.5e-3;  // kg/m3 (Rittmann and McCarty, 2001) at 20 degrees
long double K_nh4_aob = 1e-3;    // kg/m3 (Rittmann and McCarty, 2001) at 20 degrees
long double K_o2_nob  = 0.68e-3; // kg/m3 (Rittmann and McCarty, 2001)
long double K_no2_nob = 1.3e-3;  //kg/m3 (Rittmann and McCarty, 2001) at 20 degrees
// yields
long double Y_het = 0.61;       // Ni et al., 2009 (Water Research) gCODx/gCODs
long double Y_eps = 0.18;       // Ni et al., 2009 (Water Research) gCOD_eps/gCODs
long double Y_eps_het = Y_eps/Y_het; 
long double Y_s_het   = 1/Y_het;
long double Y_o2_het  = (1-Y_het-Y_eps)/Y_het;
long double Y_no2_het = (1-Y_het-Y_eps)/1.71/Y_het;
long double Y_no3_het = (1-Y_het-Y_eps)/2.86/Y_het;
long double eta_het   = 0.6;    // reduction factor in anoxic conditions ASM 2d
long double Y_aob     = 0.33;   //(Rittmann and McCarty, 2001) at 20 degrees
long double Y_o2_aob  = (3.42-Y_aob)/Y_aob;
long double Y_nh4_aob = 1/Y_aob;
long double Y_no2_aob = 1/Y_aob;
long double Y_nob     = 0.083;   //(Rittmann and McCarty, 2001) at 20 degrees
long double Y_o2_nob  = (1.15-Y_nob)/Y_nob;
long double Y_no2_nob = 1/Y_nob;
long double Y_no3_nob = 1/Y_nob;

// decay constants
long double Y_I   = 0.4;  // Alpkvist et al., 2006
long double b_het = 0.15; // 1/d Rittman and McCarty, 2001 at 20 degrees
long double b_nob = 0.11; // 1/d Rittmann and McCarty, 2001 at 20 degrees
long double b_aob = 0.11; // 1/d Rittmann and McCarty, 2001 at 20 degrees
long double b_eps = 0.17; // 1/d Ni et al., 2009 (Water Research)


//       1.7) Initializating vectors with cells information

long double bac_x[] = {0, -dx, dx, dx, -dx};  /* x's position of each cell*/
long double bac_y[] = {0, dy, -dy, dy, -dy};  /* y's position of each cell*/
long double bac_m[] = {mass_het,mass_aob,mass_nob,0,0};  /* mass of each cell*/
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



//       1.6) Parameters single attach process

				   /* pHET pAOB pNOB */ 
const float prob[] = {0.6, 0.2, 0.2};  /* probabilities for the new cell in single attach process*/
const long double mass_new_cell[] = {mass_het,mass_aob,mass_nob}; 
const long double eps_new_cell[]={mass_eps,0,0};	

//       1.9) Parameters for growth_decay

/* Initializating vectors */
vector<long double> rates_growth,rates_growth_eps,rates_decay; 

						/*   HET  AOB  NOB	 */
long double total_decay[3] = { 0 , 0 ,  0 };     /* this vector represents all the rates for HET,AOB,NOB as group and not as individual cells */
vector<long double> Total_decay (total_decay, total_decay + sizeof(total_decay) / sizeof(total_decay[0]));
		
					/*   HET  AOB  NOB	 */
long double m_decay[3] = { 0 , 0 , 0 };
vector<long double> M_decay (m_decay, m_decay + sizeof(m_decay) / sizeof(m_decay[0]));

/* Vector of vectors with the indexs of each type of cell */		
vector<vector<double> > Index_cells(4,vector<double>(0,0.0)); /* 1st vector -> HET; 2nd -> AOB; 3rd -> NOB*; 4rd -> EPS */



int main(){
	
	
	for (unsigned int i=0; i<Bac_x.size();++i){
		Bac_r.push_back(sqrt((Bac_m[i]/bac_rho + Bac_e_d[i]/bac_rhoe)/PI/bac_h));
		Bac_ra.push_back(sqrt((Bac_m[i]/bac_rho)/PI/bac_h));
	}

	for (unsigned int i=0; i<Bac_x.size();++i){				/* changing from cartesian to polar coordinates*/	
		Bac_theta.push_back(atan2(Bac_y[i],Bac_x[i]));
		Bac_R.push_back (sqrt(pow(Bac_x[i],2) + pow(Bac_y[i],2)));
	}
	
	
	vector<long double> reset_rates(Bac_x.size(),0.0); 
	rates_growth = reset_rates;
	rates_growth_eps = reset_rates;
	rates_decay = reset_rates;
	vector<vector<double> > reset_index(4,vector<double>(0,0.0));
	Index_cells = reset_index;
	
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
	
	growth_decay();
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

void growth_decay()
{
	// -------INPUTS----------
	
	// -------OUTPUTS---------
	/* "Bac_x", "Bac_y" ... are modified inside the function*/
	
	// -------EXPLANATION----------	
	/* The growth and decay processes are carried out */	

	
	// Calculate rates of each cell, and index of each type of cell
	for (unsigned int i=0; i<Bac_x.size();++i){
		calculate_rates_index(i);								
	}
	
	// Calculate the mass of decay for het,aob,nob
	/* The ODE is solved by the euler method */ /* ODE: d(MassDecay_i)/dt = rate_decay_i  i belongs {HET,AOB,NOB} */
	for (unsigned i=0; i<M_decay.size();++i){
		M_decay[i] += Total_decay[i]*dt;
	}
	
	// Eliminate random cells if the total mass decay exceeds the normal mass		
	for (unsigned i=0; i<M_decay.size();++i){
		while (M_decay[i]>mass_cells[i] && Index_cells[i].size()>0){
			int ran = ((double) rand()/(RAND_MAX))*Index_cells[Bac_s[i]].size(); /* generating an integer random number between 0 and Index_cells[Bac_s[i]].size() */
			int decay_index = Index_cells[i][ran]; /* We have choose a random index of that kind of cells*/
			if (Bac_m[decay_index]!=0){
				decay_one_cell(decay_index,ran); /* That cell will be dead */
			}
		}
	}
	
	// Eliminate EPS granule if it is completely disolved in water
	
	for (int i=Index_cells[3].size()-1;i>=0;--i){ /* Just for the cells Bac_s=3 (EPS) */
		if (Bac_e_d[Index_cells[3][i]] <= 0) /* it means it is completely disolved*/
			eliminate_cell(Index_cells[3][i]);
	}
	
	// Biomass growth
	/* Solving the ODE with Euler method for each cell*/		
	for (unsigned int i=0; i< Bac_m.size();++i){
		Bac_m[i] += rates_growth[i]*dt; 
		Bac_e_d[i] += rates_growth_eps[i]*dt;
	}
	// Recalculating radium for each cell
	for (unsigned int i=0; i<4;++i){
		for (unsigned int j=0; j<Index_cells[i].size();++j){
			int indx = Index_cells[i][j]; 
			if (Bac_s[indx] == 0||3){
				Bac_r[indx] = (sqrt((Bac_m[indx]/bac_rho + Bac_e_d[indx]/bac_rhoe)/PI/bac_h));
			}if(Bac_s[indx] == 1||2){
				Bac_r[indx] = (sqrt((Bac_m[indx]/bac_rho )/PI/bac_h));
			}if(Bac_s[indx] == 4){
				Bac_r[indx] = sqrt(Bac_e_d[indx]/bac_rhod/PI/bac_h);
			}
			Bac_ra[indx] = (sqrt((Bac_m[indx]/bac_rho)/PI/bac_h));
		}
	}
}

void calculate_rates_index(int i){
	// -------INPUTS----------
	/* "i" is the index of the cell. Its rates are calculated*/
	
	// -------OUTPUTS---------
	/* "rates_growth", "rates_growth_eps" and "rates_decay" are modified inside the function
	   "Index_cells" is calculated in this function */
	
	// -------EXPLANATION----------	
	/* Calculate rates of growth, decay, and growth solid debries for the cell "i" in the floc.
	 The concentrations of the soluble components of the up-right point (ii,jj) in the grid is used */
		
	int jj = floor(Bac_x[i]/dx)+1;
	int ii = floor(Bac_y[i]/dy)+1;
	if (Bac_s[i] == 0){
		rates_growth[i] = (bac_m[i]*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_o2[ii][jj]/(K_o2_het+c_o2[ii][jj]) + bac_m[i]*eta_het*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_no2[ii][jj]/(K_no2_het+c_no2[ii][jj])*K_o2_het/(K_o2_het+c_o2[ii][jj]) + bac_m[i]*eta_het*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_no3[ii][jj]/(K_no3_het+c_no3[ii][jj])*K_o2_het/(K_o2_het+c_o2[ii][jj]));
		rates_growth_eps[i] = (rates_growth[i]*Y_eps_het);
		rates_decay[i] = (b_het*Bac_m[i]);
	}
	if (Bac_s[i] == 1){
		rates_growth[i] = (bac_m[i]*mu_max_aob*c_nh4[ii][jj]/(K_nh4_aob+c_nh4[ii][jj])*c_o2[ii][jj]/(K_o2_aob+c_o2[ii][jj]));
		rates_decay[i] = (b_aob*Bac_m[i]);
	}
	if (Bac_s[i] == 2){
		rates_growth[i] = (bac_m[i]*mu_max_nob*c_no2[ii][jj]/(K_no2_nob+c_no2[ii][jj])*c_o2[ii][jj]/(K_o2_nob+c_o2[ii][jj]));
		rates_decay[i] = (b_nob*Bac_m[i]);
	}
	if (Bac_s[i] == 3){
		rates_decay[i] = (b_eps*Bac_e_d[i]);
		Index_cells[3].push_back(i);
	}else if (Bac_s[i]!=4){		
		Total_decay[Bac_s[i]] += rates_decay[i];
		Index_cells[Bac_s[i]].push_back(i);
	}			
}
void decay_one_cell(int a,int ran){
	// -------INPUTS----------
	/* "a" is the index of the cell to die in the decay process */
	/* "ran" is the position of "a" in the vector with indexs of each type*/
	
	// -------OUTPUTS---------
	/* Different parameters are modified inside the function*/
	
	int type;
	type = Bac_s[a];
	
	M_decay[type] -= Bac_m[a];
	m_substrate += Bac_m[a] * (1-Y_I); /* the mass is divided between the substrate released in the bioreactor and the solid debries*/
	if (type==0){
		m_substrate += Bac_e_d[a];
	}
	Bac_e_d[a] = Bac_m[a] * Y_I;
	Bac_m[a] = 0;
	rates_growth[a] = 0;
	Bac_s[a] = 4; /* dead cell */
	Bac_c[a] = 0;
	Index_cells[type].erase(Index_cells[type].begin()+ran); /* It is no longer that kind of cell, now is a dead cell */

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
