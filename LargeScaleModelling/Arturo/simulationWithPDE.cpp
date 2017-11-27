// ------------------------------------------------ Master Project ------------------------------------------ 
// ------------------------------------------------ Model with PDE ------------------------------------------
// --------------------------------------------Arturo Alvarez-Arenas ----------------------------------------
// ------------------------------------------- Newcastle University -----------------------------------------

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
#include <fstream> // create files

//reading files
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <algorithm> // for std::copy
#define PI 3.1415926535897932384626433832795
//PDE
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

using namespace std;

// create a rng in order to create random numbers
gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);

// .........GENERAL FUNCTIONS........

long length(float vector[])
{
	return sizeof(vector)/sizeof(*vector);
}


//   0) ...........INITIALIZATE FUNCTIONS.................


			void attachment_single_cell(vector<long double> &Bac_m,vector<long double> &Bac_e_d,vector<int> &Bac_s,vector<int> &Bac_c, vector<long double> &Bac_r,double Rmax);
			void attachment_group_cell(double Rmax,vector<long double> &Bacatt_x ,vector<long double> &Bacatt_y,vector<long double> &Bacatt_m,vector<long double> &Bacatt_e_d,vector<int> &Bacatt_s,vector<int> &Bacatt_c,vector<long double> &Bacatt_r,vector<long double> &Bacatt_ra,vector<double> &Bacatt_theta,vector<double> &Bacatt_R);
			void detachment_group_cell(double Rmax);
			void eliminate_cell(int index); 
			void removing_islands();
			void growth_decay();
			void calculate_rates_index(int i);
			void decay_one_cell(int a,int ran);
			void exc_eps_division();
			void EPS_excretion(int i);
			void Division(int i);
			void calculate_particulate_reactor();
			void moving_cells();
			
			void calculate_soluble_floc();
			void mass_to_grid();
			void rates_soluble();
			void elliptic(double h, double k, double aa,double bb , std::vector< vector<double> > rMatrix,vector<vector<double> > &c_something,int type_soluble );
			
			void calculate_soluble_reactor();
            void rates_soluble_reactor();
            long double boundary_values(int type_soluble);
//   1) .............SETTING PARAMETERS............

//       1.1) Bioreactor and separator sizes & others

			const short nGrowth = 100; /* each interval time stands for 0.1 days */
			const short nSingleAttach = 1;
			const short nGroupAttach = 2;
			const short nGroupDetach = 2*nGroupAttach+1;
			
			const int Nx = 50;       //grid nodes in x direction
			const int Ny = 50;       // grid nodes in y direction
			const long double dx = 1.1e-6;
			const long double dy = 1.1e-6;
			const long double Lx = Nx*dx;  //% m, space length in x; defining a grid/matrix for placing the bacteria before sending to COMSOL
			const long double Ly = Ny*dy;  //% m, space height in y

			const long double dt = 0.1;
			
			int V_reactor = 1;  /* m^2*/
			
			double alpha = 0.2;    // recycling fraction
			double beta  = 0.01;   // purje (wastege line) fraction

			long double HRT = 0.3;    // days
			long double SDT = 5.1;  // days
			
			long double mBiomass = 0.06;

			
//       1.2) Properties initial bacterium

			const float bac_h    = 1.1e-6; //=dx
			const int bac_rho = 100;
			const int bac_rhod = 105;
			const int bac_rhoe = 25;

			const long double mass_het = 1e-16;
			const long double mass_aob = 1e-16;
			const long double mass_nob = 1e-16;
			const long double mass_eps = 2.6e-17;   // kg for grid 1e-6 % for having eps density 25 kg/m3
			const long double mass_inert = 1.1e-16; // kg for grid 1e-6 % for having eps density 105 kg/m3

			const long double mass_cells[] = {mass_het,mass_aob,mass_nob,mass_eps,mass_inert}; 
			
			long number_bacterium = 5; /* Total number of bacterium, it will be modified during this process*/



//       1.3) Diffusion coefficients in water 
			
			const double D_o2  = 2e-9*3600*24;     /* m^2/d */
			const double D_nh4 = 1.4e-9*3600*24;   /* m^2/d */
			const double D_no2 = 1.2e-9*3600*24;   /* m^2/d */
			const double D_no3 = 1.2e-9*3600*24;   /* m^2/d */
			const double D_s   = 0.5e-9*3600*24;   /* m^2/d */
						
			const double Diffusion_Coefficients[] = {D_o2,D_nh4,D_no2,D_no3,D_s}; /* stored in an array */
			 	
//       1.4) Kinetic parameters 
			// kinetic parameters
			//long double mu_max_het = 6;     // 1/d miu_max HET  (Henze ASM2d and others - but consider that here there is no storage
			long double mu_max_het = 2;     // 1/d miu_max HET  (Henze ASM2d and others - but consider that here there is no storage
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
			long double Y_eps_het = 0.2951;        //  Y_eps/Y_het; 
			long double Y_s_het   = 1.6393;        //   1/Y_het;
			long double Y_o2_het  =  0.3443;             // (1-Y_het-Y_eps)/Y_het;
			long double Y_no2_het =   0.2013;           // (1-Y_het-Y_eps)/1.71/Y_het;
			long double Y_no3_het =  0.1204;             // (1-Y_het-Y_eps)/2.86/Y_het;
			long double eta_het   = 0.6;    // reduction factor in anoxic conditions ASM 2d
			long double Y_aob     = 0.33;   //(Rittmann and McCarty, 2001) at 20 degrees
			long double Y_o2_aob  = 9.3636;		// (3.42-Y_aob)/Y_aob;
			long double Y_nh4_aob = 3.0303;		// 1/Y_aob;
			long double Y_no2_aob = 3.0303;		// 1/Y_aob;
			long double Y_nob     = 0.083;   //(Rittmann and McCarty, 2001) at 20 degrees
			long double Y_o2_nob  = 12.8554;		// (1.15-Y_nob)/Y_nob;
			long double Y_no2_nob = 12.0482;		// 1/Y_nob;
			long double Y_no3_nob = 12.0482;		// 1/Y_nob;
			// decay constants
			long double Y_I   = 0.4;  // Alpkvist et al., 2006
			long double b_het = 0.15; // 1/d Rittman and McCarty, 2001 at 20 degrees
			long double b_nob = 0.11; // 1/d Rittmann and McCarty, 2001 at 20 degrees
			long double b_aob = 0.11; // 1/d Rittmann and McCarty, 2001 at 20 degrees
			long double b_eps = 0.17; // 1/d Ni et al., 2009 (Water Research)
			
//       1.6) Initial concentrations
			
		    /* concentrations of soluble components in the floc */			
			float cf_o2 = 0;
			float cf_nh4 = 0;
			float cf_no2 = 0;
			float cf_no3 = 0;
			float cf_s = 0;
			
			float soluble_floc[] = {cf_o2,cf_nh4,cf_no2,cf_no3,cf_s};  /* stored in an array*/
			
			/* Initial concentrations of soluble components in the inlet*/
			
			const float So2_in = 0.005;
			const float Snh4_in = 0.04;
			const float Sno2_in = 0.001;
			const float Sno3_in = 0.001;
			const float Ss_in = 0.04;
						
			float si_in[] = {So2_in,Snh4_in,Sno2_in,Sno3_in,Ss_in};	/* stored in an array*/
			vector<float> Si_in (si_in,si_in + sizeof(si_in)/sizeof(si_in[0])); 

			const float q_o2 = 0.5;
			const float c_o2_sat = 9e-3;
			
//       1.7) Initializating vectors with cells information

			long double bac_x[] = {0 +Lx/2, -dx+Lx/2, dx+Lx/2, dx+Lx/2, -dx+Lx/2};  /* x's position of each cell*/
			long double bac_y[] = {0+Ly/2, dy+Ly/2, -dy+Ly/2, dy+Ly/2, -dy+Ly/2};  /* y's position of each cell*/
			long double bac_m[] = {mass_het,mass_aob,mass_nob,0,0};  /* mass of each cell*/
			long double bac_e_d[] = {mass_eps,0,0,mass_eps,mass_inert}; /* mass of the eps and solid debries*/
			int bac_s[] = {0,1,2,3,4};		  /* type of cell 1-HET 2-AOB 3-NOB 4-EPS 5-DEAD*/	
			int bac_c[] = {0,1,2,0,0};        /* number of the cluster that the cells belongs to*/
			const char* bac_inf[] = {"INITIAL","INITIAL","INITIAL","INITIAL","INITIAL"};
			
			vector<long double> Bac_x (bac_x, bac_x + sizeof(bac_x) / sizeof(bac_x[0]) );
			vector<long double> Bac_y (bac_y, bac_y + sizeof(bac_y) / sizeof(bac_y[0]) );
			vector<long double> Bac_m (bac_m, bac_m + sizeof(bac_m) / sizeof(bac_m[0]) );
			vector<long double> Bac_e_d (bac_e_d, bac_e_d + sizeof(bac_e_d) / sizeof(bac_e_d[0]) );
			vector<int> Bac_s (bac_s, bac_s + sizeof(bac_s) / sizeof(bac_s[0]) );
			vector<int> Bac_c (bac_c, bac_c + sizeof(bac_c) / sizeof(bac_c[0]) );
			vector<string>  Bac_inf (bac_inf, bac_inf + sizeof(bac_inf) / sizeof(bac_inf[0]) );
			vector<long double> Bac_r;
			vector<long double> Bac_ra;

			vector<double> Bac_theta; /* This angles are in the interval [-pi pi]*/
			vector<double> Bac_R;		
			vector<double> nFlocs;
			
			
//       1.6) Parameters single attach process

						       /* pHET pAOB pNOB */ 
			const float prob[] = {0.5, 0.25, 0.25};  /* probabilities for the new cell in single attach process*/
			const long double mass_new_cell[] = {mass_het,mass_aob,mass_nob}; 
			const long double eps_new_cell[]={mass_eps,0,0};
			
//       1.6) Parameters group attach process

			vector<long double> Bacatt_x ;
			vector<long double> Bacatt_y ;
			vector<long double> Bacatt_m ;
			vector<long double> Bacatt_e_d ;
			vector<int> Bacatt_s ;
			vector<int> Bacatt_c ;
			vector<long double> Bacatt_r;
			vector<long double> Bacatt_ra;

			vector<double> Bacatt_theta;
			vector<double> Bacatt_R;

//       1.9) Parameters for growth_decay

			/* Initializating vectors */
			vector<long double> rates_growth,rates_growth_eps,rates_decay; 

									/*   HET  AOB  NOB	 */
			long double total_decay[3] = { 0 , 0 ,  0 };     /* this vector represents all the rates for HET,AOB,NOB as group and not as individual cells */
			vector<long double> Total_decay (total_decay, total_decay + sizeof(total_decay) / sizeof(total_decay[0]));
					
								/*   HET  AOB  NOB	 */
			long double m_decay[3] = { 0 , 0 , 0 };
			vector<long double> M_decay (m_decay, m_decay + sizeof(m_decay) / sizeof(m_decay[0]));

			long double m_substrate = 0;

			/* Vector of vectors with the indexs of each type of cell */		
			vector<vector<double> > Index_cells(4,vector<double>(0,0.0)); /* 1st vector -> HET; 2nd -> AOB; 3rd -> NOB*; 4rd -> EPS */
			vector<vector<double> > c_s(Nx,vector<double>(Nx,0.08));
			vector<vector<double> > c_o2(Nx,vector<double>(Nx,0.005));
			vector<vector<double> > c_no2(Nx,vector<double>(Nx,0.04));
			vector<vector<double> > c_no3(Nx,vector<double>(Nx,0));
			vector<vector<double> > c_nh4(Nx,vector<double>(Nx,0));

//       1.10) Parameters for division

			const long double bac_mmax[] = {2*mass_het,2*mass_aob,2*mass_nob,100,100}; /* The last two elements are put in order to run the loop in the division processi without diferentiating types of cells */
			const float DivisionFraction = 0.5;

//       1.11) Macro scale 
			vector<vector<double> > mass_reactor(5,vector<double>(0,0.0)); /* concentration of each type in each time */
			vector<vector<double> > soluble_reactor(5,vector<double>(0,0.0)); /* concentration of each type in each time */
		
			long double rR_o2,rR_nh4,rR_no2,rR_no3,rR_s;
//       1.10) Parameters for soluble reactor

			vector<vector<double> > r_s(Nx,vector<double>(Nx,0));
			vector<vector<double> > r_o2(Nx,vector<double>(Nx,0));
			vector<vector<double> > r_no2(Nx,vector<double>(Nx,0));
			vector<vector<double> > r_no3(Nx,vector<double>(Nx,0));
			vector<vector<double> > r_nh4(Nx,vector<double>(Nx,0));

			vector<vector<double> > c_bac_HET(Nx,vector<double>(Nx,0));
			vector<vector<double> > c_bac_AOB(Nx,vector<double>(Nx,0));
			vector<vector<double> > c_bac_NOB(Nx,vector<double>(Nx,0));
			vector<vector<double> > c_bac_EPS(Nx,vector<double>(Nx,0));
			vector<vector<double> > c_bac_DEAD(Nx,vector<double>(Nx,0));

//   2) .................MAIN LOOP.................
int main ()
{


	soluble_reactor[0].push_back(0.005);
	soluble_reactor[1].push_back(0.04);
	soluble_reactor[2].push_back(0.001);
	soluble_reactor[3].push_back(0.001);
	soluble_reactor[4].push_back(0.04);	
			
	srand(time(NULL));
	/* Calculating some Values*/

	long double initialMass=0;
	for(unsigned int i=0; i<Bac_m.size(); ++i){
		initialMass += Bac_m[i] + Bac_e_d[i];
	}
	
	nFlocs.push_back(mBiomass/initialMass);
	for (unsigned int i=0; i<Bac_x.size();++i){
		Bac_r.push_back(sqrt((Bac_m[i]/bac_rho + Bac_e_d[i]/bac_rhoe)/PI/bac_h));
		Bac_ra.push_back(sqrt((Bac_m[i]/bac_rho)/PI/bac_h));
	}	
	
	// reading the vectors for the attachment group of cells
	ifstream is("Bac_x_attach.txt");
    istream_iterator<long double> start(is), end;
    vector<long double> Bacatt_x(start, end);
    
	ifstream is1("Bac_y_attach.txt");
    istream_iterator<long double> start1(is1), end1;
    vector<long double> Bacatt_y(start1, end1);
    
	ifstream is3("Bac_c_attach.txt");
    istream_iterator<int> start3(is3), end3;
    vector<int> Bacatt_c(start3, end3);
    
	ifstream is5("Bac_m_attach.txt");
    istream_iterator<long double> start5(is5), end5;
    vector<long double> Bacatt_m(start5, end5);
    
	ifstream is6("Bac_e_d_attach.txt");
    istream_iterator<long double> start6(is6), end6;
    vector<long double> Bacatt_e_d(start6, end6);
    
	ifstream is7("Bac_s_attach.txt");
    istream_iterator<int> start7(is7), end7;
    vector<int> Bacatt_s(start7, end7);   
          
	for (unsigned int i=0; i<Bacatt_x.size();++i){
		Bacatt_r.push_back(sqrt((Bacatt_m[i]/bac_rho + Bacatt_e_d[i]/bac_rhoe)/PI/bac_h));
		Bacatt_ra.push_back(sqrt((Bacatt_m[i]/bac_rho)/PI/bac_h));
	}

	for (unsigned int i=0; i<Bacatt_x.size();++i){				/* changing from cartesian to polar coordinates*/	
		Bacatt_theta.push_back(atan2(Bacatt_y[i],Bacatt_x[i]));
		Bacatt_R.push_back (sqrt(pow(Bacatt_x[i],2) + pow(Bacatt_y[i],2)));
	}


// --------------------THE SIMULATION STARTS HERE ------------		
	//for (short iterGrowth=1; iterGrowth<=nGrowth; iterGrowth++){
	for (short iterGrowth=1; iterGrowth<=40; iterGrowth++){
// ---------------------------------------------MICRO-SCALE MODEL--------------------------------------------------------
			cout << "\n\n\n\n iterGrowth: " << iterGrowth;	
			cout << "\n\n";							
//       2.1) Soluble concentrations floc 
			
			calculate_soluble_floc();

//       2.2) Saving the actual mass of the floc
			
			//long double mFlocOld = boost::accumulate(Bac_m,0) + boost::accumulate(Bac_e_d,0); /* It will be used when calculating number of flocs */
			long double mFlocOld=0;
			for(unsigned int i=0; i<Bac_m.size(); ++i){
				mFlocOld += Bac_m[i] + Bac_e_d[i];
			}
			/* Moving the cells from the center of the floc to the left-down corner*/																																																																																							
			for (unsigned int i=0; i<Bac_x.size();++i){
				Bac_x[i] -= Lx/2;
				Bac_y[i] -= Ly/2;
			}
			
			Bac_theta.clear();
			Bac_R.clear();
			/* Recalculated in each step Bac_theta and Bac_R*/ 
			for (unsigned int i=0; i<Bac_x.size();++i){				/* changing from cartesian to polar coordinates*/	
				Bac_theta.push_back(atan2(Bac_y[i],Bac_x[i]));
				Bac_R.push_back (sqrt(pow(Bac_x[i],2) + pow(Bac_y[i],2)));
			}
			/* calculating Rmax*/	
			vector<double>::const_iterator it, RmaxVector;
			long double Rmax;
			RmaxVector = max_element(Bac_R.begin(),Bac_R.end());
			Rmax = *RmaxVector;





//       2.3) Single attach process
		    if ( iterGrowth % nSingleAttach == 0 )
			{
				attachment_single_cell(Bac_m,Bac_e_d,Bac_s,Bac_c,Bac_r,Rmax);							
			}	
		
		
//       2.4) Group attach process		
			if ( iterGrowth % nGroupAttach == 0 )
			{
				attachment_group_cell(Rmax,Bacatt_x ,Bacatt_y,Bacatt_m,Bacatt_e_d,Bacatt_s,Bacatt_c,Bacatt_r,Bacatt_ra,Bacatt_theta,Bacatt_R);
			}
			
		
//       2.5) Group dettach process
			if ( iterGrowth % nGroupDetach == 0 )
			{
				detachment_group_cell(Rmax);		
			}			
		
//       2.6) Moving back the cells to the center of the square		
			for (unsigned int i=0; i<Bac_x.size();++i){
				Bac_x[i] += Lx/2;
				Bac_y[i] += Ly/2;
			}
//       2.7) Removing islands process
			if ( iterGrowth % nGroupDetach == 1 )
			{
				removing_islands();			
			}

//       2.8) Saving mass of new floc in order to calculate the nFlocs afterwards
			long double mFlocNew=0;
			for(unsigned int i=0; i<Bac_m.size(); ++i){
				mFlocNew += Bac_m[i] + Bac_e_d[i];
			}

//       2.9) Decay process
			vector<long double> reset_rates(Bac_x.size(),0.0); 
			rates_growth = reset_rates;
			rates_growth_eps = reset_rates;
			rates_decay = reset_rates;
			vector<vector<double> > reset_index(5,vector<double>(0,0.0));
			Index_cells = reset_index;
		
		
			growth_decay();		

//       2.10) Division and EPS excretion process
			exc_eps_division();

		

// ---------------------------------------------MACRO-SCALE MODEL--------------------------------------------------------


//       2.11) Calculating number of Flocs	
			nFlocs.push_back(nFlocs.back()*(1- mFlocOld*dt/(mFlocNew*SDT)));	// mFlocOld = mFlocNew = 0	

//       2.12) Soluble concentrations bioreactor

			calculate_soluble_reactor();
			
//       2.13) Particulate concentrations bioreactor

			calculate_particulate_reactor();			

// ------------------------------------------MOVING CELLS------------------------------------------------------------

//       2.14) Repositioning cells 

			moving_cells();

	/* Putting data in files to be read by R*/
	ofstream myfile;
	myfile.open("c_s.txt");
	for (unsigned int j=0; j<Nx;++j){
		myfile <<"\n";
		for (unsigned int i=0; i<Nx;++i){
			myfile << c_s[j][i];
			myfile << " ";
		}
	}
	myfile.close();	
	myfile.open("c_no2.txt");
	for (unsigned int j=0; j<Nx;++j){
		myfile <<"\n";
		for (unsigned int i=0; i<Nx;++i){
			myfile << c_no2[j][i];
			myfile << " ";
		}
	}
	myfile.close();
	myfile.open("c_no3.txt");
	for (unsigned int j=0; j<Nx;++j){
		myfile <<"\n";
		for (unsigned int i=0; i<Nx;++i){
			myfile << c_no3[j][i];
			myfile << " ";
		}
	}
	myfile.close();
	myfile.open("c_o2.txt");
	for (unsigned int j=0; j<Nx;++j){
		myfile <<"\n";
		for (unsigned int i=0; i<Nx;++i){
			myfile << c_o2[j][i];
			myfile << " ";
		}
	}
	myfile.close();
	myfile.open("c_nh4.txt");
	for (unsigned int j=0; j<Nx;++j){
		myfile <<"\n";
		for (unsigned int i=0; i<Nx;++i){
			myfile << c_nh4[j][i];
			myfile << " ";
		}
	}
	myfile.close();
	myfile.open("soluble_reactorO2.txt");
	for (unsigned int i=0; i<soluble_reactor[0].size();++i){
		myfile << soluble_reactor[0][i];
		myfile << " ";
	}
	myfile.close();

	myfile.open("soluble_reactorNH4.txt");
	for (unsigned int i=0; i<soluble_reactor[1].size();++i){
		myfile << soluble_reactor[1][i];
		myfile << " ";
	}
	myfile.close();
	
		myfile.open("soluble_reactorNO2.txt");
	for (unsigned int i=0; i<soluble_reactor[2].size();++i){
		myfile << soluble_reactor[2][i];
		myfile << " ";
	}
	myfile.close();
	
		myfile.open("soluble_reactorNO3.txt");
	for (unsigned int i=0; i<soluble_reactor[3].size();++i){
		myfile << soluble_reactor[3][i];
		myfile << " ";
	}
	myfile.close();
	
		myfile.open("soluble_reactorS.txt");
	for (unsigned int i=0; i<soluble_reactor[4].size();++i){
		myfile << soluble_reactor[4][i];
		myfile << " ";
	}
	myfile.close();

	myfile.open("mass_reactorHET.txt");
	for (unsigned int i=0; i<mass_reactor[0].size();++i){
		myfile << mass_reactor[0][i];
		myfile << " ";
	}
	myfile.close();

	myfile.open("mass_reactorAOB.txt");
	for (unsigned int i=0; i<mass_reactor[1].size();++i){
		myfile << mass_reactor[1][i];
		myfile << " ";
	}
	myfile.close();
	
		myfile.open("mass_reactorNOB.txt");
	for (unsigned int i=0; i<mass_reactor[2].size();++i){
		myfile << mass_reactor[2][i];
		myfile << " ";
	}
	myfile.close();
	
		myfile.open("mass_reactorEPS.txt");
	for (unsigned int i=0; i<mass_reactor[3].size();++i){
		myfile << mass_reactor[3][i];
		myfile << " ";
	}
	myfile.close();
	
	myfile.open("mass_reactorDEAD.txt");
	for (unsigned int i=0; i<mass_reactor[4].size();++i){
		myfile << mass_reactor[4][i];
		myfile << " ";
	}
	myfile.close();
	
	myfile.open("Bac_x.txt");
	for (unsigned int i=0; i<Bac_x.size();++i){
		myfile << Bac_x[i];
		myfile << " ";
	}
	myfile.close();

	myfile.open("Bac_y.txt");
	for (unsigned int i=0; i<Bac_y.size();++i){
		myfile << Bac_y[i];
		myfile << " ";
	}
	myfile.close();
	
	myfile.open("Bac_r.txt");
	for (unsigned int i=0; i<Bac_r.size();++i){
		myfile << Bac_r[i];
		myfile << " ";
	}
	myfile.close();
	
	myfile.open("Bac_s.txt");
	for (unsigned int i=0; i<Bac_s.size();++i){
		myfile << Bac_s[i];
		myfile << " ";
	}
	myfile.close();

	myfile.open("Bac_m.txt");
	for (unsigned int i=0; i<Bac_m.size();++i){
		myfile << Bac_m[i];
		myfile << " ";
	}
	myfile.close();
	
	myfile.open("Bac_e_d.txt");
	for (unsigned int i=0; i<Bac_e_d.size();++i){
		myfile << Bac_e_d[i];
		myfile << " ";
	}
	myfile.close();

	myfile.open("Bac_ra.txt");
	for (unsigned int i=0; i<Bac_ra.size();++i){
		myfile << Bac_ra[i];
		myfile << " ";
	}
	myfile.close();
	
	myfile.open("Bac_c.txt");
	for (unsigned int i=0; i<Bac_c.size();++i){
		myfile << Bac_c[i];
		myfile << " ";
	}
	myfile.close();	
	
	myfile.open("Bac_inf.txt");
	for (unsigned int i=0; i<Bac_c.size();++i){
		myfile << Bac_inf[i];
		myfile << " ";
	}
	myfile.close();	
	
	myfile.open("nFlocs.txt");
	for (unsigned int i=0; i<nFlocs.size();++i){
		myfile << nFlocs[i];
		myfile << " ";
	}
	myfile.close();
	} /* END OF THE LOOP*/
	
	
/*
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
		cout << "\n Bac_inf: " << Bac_inf[i];
	}*/

	cout << "\n\n number of cells: " << Bac_x.size();
	int dead=0;
	for (unsigned int i=0; i<Bac_x.size();++i){
		if (Bac_s[i]==4){
			dead +=1;
		}
	}
	cout << "\n\n number of dead cells: " << dead;

	cout << "\n\n Nflocs: " << nFlocs.back();
	
	cout << "\n\n concentrations particulate reactor";
	cout << "\n HET : "<< mass_reactor[0].back();
	cout << "\n AOB : "<< mass_reactor[1].back();
	cout << "\n NOB : "<< mass_reactor[2].back();
	cout << "\n EPS : "<< mass_reactor[3].back();
	cout << "\n DEAD : "<< mass_reactor[4].back();


} 
// ---------------------------------------------END OF THE MAIN LOOP--------------------------------------------------------

	

//   4) ..............DEFINED FUNCTIONS.........

// ----------------------------------------------------attachment_single_cell-----------------------------------------------------	


				
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
	
	Bac_inf.push_back("att_single_cell");
	Bac_r.push_back(sqrt((Bac_m.back()/bac_rho + Bac_e_d.back()/bac_rhoe)/PI/bac_h));
	Bac_ra.push_back(sqrt((Bac_m.back()/bac_rho)/PI/bac_h));	
	Bac_theta.push_back(atan2(Bac_y.back(),Bac_x.back()));
	Bac_R.push_back (sqrt(pow(Bac_x.back(),2) + pow(Bac_y.back(),2)));
	
	++number_bacterium;
}

// ----------------------------------------------------attachment_group_cell-----------------------------------------------------					
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
		Bac_inf.push_back("att_group_cell");
	}	
}

// ----------------------------------------------------detachment_group_cell-----------------------------------------------------				


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
			gsl_rng_set (r, rand() + 1);
			// generates a uniform in [-PI PI]	
			double fi_atta = gsl_rng_uniform (r)*2*PI - PI;
			
			for (unsigned int i=0; i<Bac_theta.size();++i){
					// if (abs(Bac_theta[i]-fi_atta)<= (PI/10)){ esta es la que he hecho yo, interpretando que esta en un angulo de +-10%
					if (abs(Bac_theta[i]-fi_atta)<= 0.2*abs(fi_atta)){  // esta es la suya, que segun mi esta mal
						ind_theta.push_back(i);
						/* For all those cells we check which one is the farthest from the center*/
						if (Bac_R[i]>= largest_radio){
							largest_radio = Bac_R[i];							
						}						
					}				
			}
	}	
	/* Select a random number in [0.5*radio_largest , radio largest] */
	

	gsl_rng_set (r, rand() + 1);
	// generates a uniform in [0.5*radio_largest , radio largest]
	double random_radio = (gsl_rng_uniform (r)*0.5 + 0.5)*largest_radio;  
	/* All the cells saved before with a radio larger than random_radio are removed */
	
	vector<int> cluster_to_remove;
	const unsigned int ini=ind_theta.size()-1;
	for (int i=ini;i>=0;--i){ /* The loop starts in the last element. So the positions of the other cells are not changed after removing one cell */
		if (Bac_R[ind_theta[i]]>=random_radio){
			if ((Bac_s[ind_theta[i]]==1)||(Bac_s[ind_theta[i]]==2)){
				cluster_to_remove.push_back(Bac_c[ind_theta[i]]); /* the cluster in which the cell was is saved before removing it*/
			}
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

// ----------------------------------------------------removing islands-----------------------------------------------------				


void removing_islands(){
	// -------INPUTS----------
	
	// -------OUTPUTS---------
	/* "Bac_x", "Bac_y" ... are modified inside the function*/
	
	// -------EXPLANATION----------	
	/* All the cells that are not in a defined "main island" will be removed from the floc */
		
	// Define a matrix (dimension Nx*Ny) with 0s if there is no cell, or 1 if there is any cell in the down-left square
	int land[Nx][Ny] = {{ 0 }};
	for (unsigned int i=0; i<Bac_x.size(); ++i){
		int jj,ii; 
			if(jj>(Nx-1)){
				jj=Nx-1;
			}
			if(jj<0){
				jj=0;
			}
			if(ii>(Nx-1)){
				ii=Nx-1;
			}
			if(ii<0){
				ii=0;
			}
		jj = (floor(Bac_x[i]/dx)+1);
		ii = (floor(Bac_y[i]/dy)+1);
		land[ii][jj] = 1;	
	}
	
	// The center of the main island is marked with 3
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


// ----------------------------------------------------growth decay-----------------------------------------------------				

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
			gsl_rng_set (r, rand() + 1);
			// generates a uniform in 0,1			 
			int ran = gsl_rng_uniform (r)*Index_cells[i].size(); /* generating an integer random number between 0 and Index_cells[Bac_s[i]].size() */
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
	for (unsigned int i=0;i<Bac_x.size();++i){
		if (Bac_s[i] == 0||3){
			Bac_r[i] = (sqrt((Bac_m[i]/bac_rho + Bac_e_d[i]/bac_rhoe)/PI/bac_h));
		}if(Bac_s[i] == 1||2){
			Bac_r[i] = (sqrt((Bac_m[i]/bac_rho )/PI/bac_h));
		}if(Bac_s[i] == 4){
			Bac_r[i] = sqrt(Bac_e_d[i]/bac_rhod/PI/bac_h);
		}
		Bac_ra[i] = (sqrt((Bac_m[i]/bac_rho)/PI/bac_h));
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
	if(jj>Nx-1){
		jj=Nx-1;
	}
	if(jj<0){
		jj=0;
	}
	if(ii>(Nx-1)){
		ii=Nx-1;
	}
	if(ii<0){
		ii=0;
	}

	if (Bac_s[i] == 0){
		rates_growth[i] = (Bac_m[i]*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_o2[ii][jj]/(K_o2_het+c_o2[ii][jj]) + Bac_m[i]*eta_het*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_no2[ii][jj]/(K_no2_het+c_no2[ii][jj])*K_o2_het/(K_o2_het+c_o2[ii][jj]) + Bac_m[i]*eta_het*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_no3[ii][jj]/(K_no3_het+c_no3[ii][jj])*K_o2_het/(K_o2_het+c_o2[ii][jj]));
		rates_growth_eps[i] = (rates_growth[i]*Y_eps_het);
		rates_decay[i] = (b_het*Bac_m[i]);
		Total_decay[Bac_s[i]] += rates_decay[i];
		Index_cells[Bac_s[i]].push_back(i);
	}
	if (Bac_s[i] == 1){
		rates_growth[i] = (Bac_m[i]*mu_max_aob*c_nh4[ii][jj]/(K_nh4_aob+c_nh4[ii][jj])*c_o2[ii][jj]/(K_o2_aob+c_o2[ii][jj]));
		rates_decay[i] = (b_aob*Bac_m[i]);
		Total_decay[Bac_s[i]] += rates_decay[i];
		Index_cells[Bac_s[i]].push_back(i);
	}
	if (Bac_s[i] == 2){
		rates_growth[i] = (Bac_m[i]*mu_max_nob*c_no2[ii][jj]/(K_no2_nob+c_no2[ii][jj])*c_o2[ii][jj]/(K_o2_nob+c_o2[ii][jj]));
		rates_decay[i] = (b_nob*Bac_m[i]);
		Total_decay[Bac_s[i]] += rates_decay[i];
		Index_cells[Bac_s[i]].push_back(i);
	}
	if (Bac_s[i] == 3){
		rates_decay[i] = (b_eps*Bac_e_d[i]);
		Index_cells[3].push_back(i);
	}
	if (Bac_s[i]==4){		
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
	int j=*Index_cells[type].begin()+ran;
	Index_cells[4].push_back(j); /* now it is a dead cell */
	Index_cells[type].erase(Index_cells[type].begin()+ran); /* It is no longer that kind of cell, now is a dead cell */

	Bac_inf[a]="decay_one_cell";
}

// ----------------------------------------------------division eps excretion-----------------------------------------------------		


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
	Bac_inf.push_back("EPS_excretion");
	/* Recalculating some properties previous cell*/
	Bac_e_d[i] *= (1-fraction);
	Bac_r[i] = sqrt((Bac_m[i]/bac_rho + Bac_e_d[i]/bac_rhoe)/PI/bac_h);
	
	/* Updating indexs*/
	Index_cells[3].push_back(Bac_x.size()-1);  
}

void Division(int i){
	// -------INPUTS----------
	/* "index" is the index of the cell that will divide*/
	
	// -------OUTPUTS---------
	/* "Bac_x", "Bac_y" ... are modified inside the function*/
	
	// -------EXPLANATION----------	
	/* the division process is carried out */
	
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
	Bac_inf.push_back("division");
	/* Moving the previous cell*/
	Bac_x[i] -= Bac_r[i]*cos(fi);
	Bac_y[i] -= Bac_r[i]*sin(fi);
	/* Recalculating mass and size*/
	Bac_m[i] *= (1-fraction1);
	Bac_e_d[i] *= (1-fraction1);
	Bac_ra[i] = sqrt((Bac_m[i]/bac_rho)/PI/bac_h);
	Bac_r[i] = sqrt((Bac_m[i]/bac_rho + Bac_e_d[i]/bac_rhoe)/PI/bac_h);
	
	/* Updating indexs*/
	Index_cells[Bac_s[i]].push_back(Bac_x.size()-1);		
}

// ----------------------------------------------------calculate particulate reactor-----------------------------------------------------				

void calculate_particulate_reactor(){
	// -------INPUTS----------
	
	// -------OUTPUTS---------
	/* "mass_reactor" is modified inside the function*/
	
	// -------EXPLANATION----------	
	/* when this function is run concentrations of particulate components in the bioreactor are saved at the end of mass_reactor */
	long double mass;
	for (int i=0; i<5;++i){
		//mass_reactor[i].push_back(boost::accumulate(Index_cells[i],0)*nFlocs.back()/V_reactor);
		mass=0;
		if (i!=3 && 4){
			for (unsigned int j=0;j<Index_cells[i].size();++j){		
				mass += Bac_m[Index_cells[i][j]];
			}
		}else{
			for (unsigned int j=0;j<Index_cells[i].size();++j){			
				mass += Bac_e_d[Index_cells[i][j]];
			}
		}
		mass_reactor[i].push_back(mass*nFlocs.back()/V_reactor);
	}	
}

// ----------------------------------------------------moving cells-----------------------------------------------------

void moving_cells(){
	// -------INPUTS----------
	
	// -------OUTPUTS---------
	/* "Bac_x" and "Bac_y" are modified inside the function */
	
	// -------EXPLANATION----------	
	/* the floc is put in the center of the separator */
	
	/* center of mass for all the cells */
	long double mtot_floc = 0;
	long double mxtot_floc = 0;
	long double mytot_floc = 0;
	for (unsigned int i=0;i<Bac_m.size();++i){
		mtot_floc += Bac_m[i] + Bac_e_d[i];
		mxtot_floc += (Bac_m[i] + Bac_e_d[i])*Bac_x[i];
		mytot_floc += (Bac_m[i] + Bac_e_d[i])*Bac_y[i];
	}
	long double xc_flocOld, yc_flocOld;
	xc_flocOld = mxtot_floc/mtot_floc;
	yc_flocOld = mytot_floc/mtot_floc;
	
	/* calculate the distance to put the center of mass in the center of the square */
	long double dxc,dyc;
	dxc = xc_flocOld - Lx/2;
	dyc = yc_flocOld - Ly/2;
	
	/* moving all the cells that distance */
	for (unsigned int i=0; i<Bac_x.size();++i){
		Bac_x[i] -= dxc;
		Bac_y[i] -= dyc;
	}
}

// ----------------------------------------------------eliminate_cell-----------------------------------------------------				


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
	Bac_inf.erase(Bac_inf.begin()+index);		
}


// ----------------------------------------------------Soluble Floc-----------------------------------------------------				


void calculate_soluble_floc(){
	// -------INPUTS----------
	
	// -------OUTPUTS---------
	/* "c_o2", "c_nh4" ... are modified inside the function*/
	
	// -------EXPLANATION----------	
	/* The concentrations of soluble components in each position of the grid is calculated */
	
	/* put the mass cells in the grid */
	mass_to_grid();
	
	/* calculate rates soluble component */
	rates_soluble();
	
	/* solve PDE */	
	
	elliptic(1,1,Nx,Nx,r_o2,c_o2,0);
	cout << "\n o2 is done";
	elliptic(1,1,Nx,Nx,r_nh4,c_nh4,1);
	cout << "\n nh4 is done";
	elliptic(1,1,Nx,Nx,r_no2,c_no2,2);
	cout << "\n no2 is done";
	elliptic(1,1,Nx,Nx,r_no3,c_no3,3);
	cout << "\n no3 is done";
	elliptic(1,1,Nx,Nx,r_s,c_s,4);
	cout << "\n s is done";
	
	// changing negatives values to 0
	
	for(unsigned int i=0;i<c_o2[0].size();++i){
		for(unsigned int j=0;j<c_o2[0].size();++j){
			if(c_o2[i][j] <0){
				c_o2[i][j] = 0;
			}
		}
	}
	for(unsigned int i=0;i<c_nh4[0].size();++i){
		for(unsigned int j=0;j<c_nh4[0].size();++j){
			if(c_nh4[i][j] <0){
				c_nh4[i][j] = 0;
			}
		}
	}
	for(unsigned int i=0;i<c_no2[0].size();++i){
		for(unsigned int j=0;j<c_no2[0].size();++j){
			if(c_no2[i][j] <0){
				c_no2[i][j] = 0;
			}
		}
	}
	for(unsigned int i=0;i<c_no3[0].size();++i){
		for(unsigned int j=0;j<c_no3[0].size();++j){
			if(c_no3[i][j] <0){
				c_no3[i][j] = 0;
			}
		}
	}
	for(unsigned int i=0;i<c_s[0].size();++i){
		for(unsigned int j=0;j<c_s[0].size();++j){
			if(c_s[i][j] <0){
				c_s[i][j] = 0;
			}
		}
	}
	
}


void mass_to_grid(){
	// -------INPUTS----------
	
	// -------OUTPUTS---------
	/* "c_bac_HET", "c_bac_AOB" ... are modified inside the function*/
	
	// -------EXPLANATION----------	
	/* The mass is put in the nodes of the grid */
	
	// Reseting vectors
	vector<vector<double> > reset_concentrations(Nx,vector<double>(Nx,0));
	c_bac_HET = reset_concentrations;
	c_bac_AOB = reset_concentrations;
	c_bac_NOB = reset_concentrations;
	c_bac_EPS = reset_concentrations;
	c_bac_DEAD = reset_concentrations;
	
	for (unsigned int i=0;i<Bac_x.size();++i){
		int jj = floor(Bac_x[i]/dx) + 1;
		int ii = floor(Bac_y[i]/dy) + 1;
			if(jj>(Nx-1)){
				jj=Nx-1;
			}
			if(jj<0){
				jj=0;
			}
			if(ii>(Nx-1)){
				ii=Nx-1;
			}
			if(ii<0){
				ii=0;
			}
		if (Bac_s[i] ==1){
			c_bac_HET[ii][jj] += Bac_m[i];
		}
		if (Bac_s[i] ==2){
			c_bac_AOB[ii][jj] += Bac_m[i];
		}
		if (Bac_s[i] ==3){
			c_bac_NOB[ii][jj] += Bac_m[i];
		}
		if (Bac_s[i] ==4){
			c_bac_EPS[ii][jj] += Bac_e_d[i];
		}
		if (Bac_s[i] ==5){
			c_bac_DEAD[ii][jj] += Bac_e_d[i];
		}						
	}
}

void rates_soluble(){
	// -------INPUTS----------
	
	// -------OUTPUTS---------
	/* "r_o2", "r_nh4" ... are modified inside the function*/
	
	// -------EXPLANATION----------	
	/* The rate for each point in the grid is calculated */

	for (int ii=0;ii<Nx;++ii){
		for (int jj=0;jj<Nx;++jj){
			r_o2[ii][jj] = -Y_o2_het*(c_bac_HET[ii][jj]*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_o2[ii][jj]/(K_o2_het+c_o2[ii][jj]) + c_bac_HET[ii][jj]*eta_het*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_no2[ii][jj]/(K_no2_het+c_no2[ii][jj])*K_o2_het/(K_o2_het+c_o2[ii][jj]) + c_bac_HET[ii][jj]*eta_het*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_no3[ii][jj]/(K_no3_het+c_no3[ii][jj])*K_o2_het/(K_o2_het+c_o2[ii][jj])) - Y_o2_aob*(c_bac_AOB[ii][jj]*mu_max_aob*c_nh4[ii][jj]/(K_nh4_aob+c_nh4[ii][jj])*c_o2[ii][jj]/(K_o2_aob+c_o2[ii][jj])) - Y_o2_nob*(c_bac_NOB[ii][jj]*mu_max_nob*c_no2[ii][jj]/(K_no2_nob+c_no2[ii][jj])*c_o2[ii][jj]/(K_o2_nob+c_o2[ii][jj]));
			r_nh4[ii][jj] = -Y_nh4_aob*(c_bac_AOB[ii][jj]*mu_max_aob*c_nh4[ii][jj]/(K_nh4_aob+c_nh4[ii][jj])*c_o2[ii][jj]/(K_o2_aob+c_o2[ii][jj]));
			r_no2[ii][jj] = Y_no2_aob * (c_bac_AOB[ii][jj]*mu_max_aob*c_nh4[ii][jj]/(K_nh4_aob+c_nh4[ii][jj])*c_o2[ii][jj]/(K_o2_aob+c_o2[ii][jj])) - Y_no2_nob*(c_bac_NOB[ii][jj]*mu_max_nob*c_no2[ii][jj]/(K_no2_nob+c_no2[ii][jj])*c_o2[ii][jj]/(K_o2_nob+c_o2[ii][jj])) - Y_no2_het*c_bac_HET[ii][jj]*eta_het*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_no2[ii][jj]/(K_no2_het+c_no2[ii][jj])*K_o2_het/(K_o2_het+c_o2[ii][jj]); 
			r_no3[ii][jj] = Y_no3_nob*(c_bac_NOB[ii][jj]*mu_max_nob*c_no2[ii][jj]/(K_no2_nob+c_no2[ii][jj])*c_o2[ii][jj]/(K_o2_nob+c_o2[ii][jj])) - Y_no3_het*c_bac_HET[ii][jj]*eta_het*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_no3[ii][jj]/(K_no3_het+c_no3[ii][jj])*K_o2_het/(K_o2_het+c_o2[ii][jj]);
			r_s[ii][jj] = -Y_s_het*(c_bac_HET[ii][jj]*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_o2[ii][jj]/(K_o2_het+c_o2[ii][jj]) + c_bac_HET[ii][jj]*eta_het*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_no2[ii][jj]/(K_no2_het+c_no2[ii][jj])*K_o2_het/(K_o2_het+c_o2[ii][jj]) + c_bac_HET[ii][jj]*eta_het*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_no3[ii][jj]/(K_no3_het+c_no3[ii][jj])*K_o2_het/(K_o2_het+c_o2[ii][jj]));
		}
	}
}

// ----------------------------------------------------calculate soluble reactor-----------------------------------------------------	
void calculate_soluble_reactor(){
	// -------INPUTS----------
	
	// -------OUTPUTS---------
	/* "soluble_reactor" is modified inside the function*/
	
	// -------EXPLANATION----------	
	/* The concentrations of the soluble components in the reactor is calculated */

	/* put the mass cells in the grid */
	mass_to_grid();
	
	/* calculate rates soluble component */
	rates_soluble_reactor();
	
	/* solve ODE */	
	long double rates_sol[] = {rR_o2,rR_nh4,rR_no2,rR_no3,rR_s};
	for (unsigned int i=0;i<5;++i){
		if (i==0){
			//if (dt*((Si_in[i] - soluble_reactor[i].back())/HRT  - rates_sol[i] + q_o2*(1-(soluble_reactor[i].back()/c_o2_sat))) + soluble_reactor[i].back()<0){
			if (dt*((Si_in[i] - soluble_reactor[i].back())/HRT  - rates_sol[i] ) + soluble_reactor[i].back()<0){
				soluble_reactor[i].push_back(0);
			}else{
				//soluble_reactor[i].push_back(dt*((Si_in[i] - soluble_reactor[i].back())/HRT  - rates_sol[i] + q_o2*(1-(soluble_reactor[i].back()/c_o2_sat))) + soluble_reactor[i].back());
				soluble_reactor[i].push_back(dt*((Si_in[i] - soluble_reactor[i].back())/HRT  - rates_sol[i] ) + soluble_reactor[i].back());
			}	
		}
		if (i==4){
			if (((dt*((Si_in[i] - soluble_reactor[i].back())/HRT  - rates_sol[i]) + soluble_reactor[i].back())+m_substrate*nFlocs.back()/V_reactor)<0){
				soluble_reactor[i].push_back(0);
			}else{
				soluble_reactor[i].push_back((dt*((Si_in[i] - soluble_reactor[i].back())/HRT  - rates_sol[i]) + soluble_reactor[i].back())+m_substrate*nFlocs.back()/V_reactor);
			}
		}	
		if (i==1 || i == 2 || i==3){
			if (dt*((Si_in[i] - soluble_reactor[i].back())/HRT  - rates_sol[i]) + soluble_reactor[i].back()<0){
				soluble_reactor[i].push_back(0);
			}else{	
				soluble_reactor[i].push_back(dt*((Si_in[i] - soluble_reactor[i].back())/HRT  - rates_sol[i]) + soluble_reactor[i].back());
			}
		}
	
	}
}
void rates_soluble_reactor(){
	// -------INPUTS----------
	
	// -------OUTPUTS---------
	/* "soluble_reactor" is modified inside the function*/
	
	// -------EXPLANATION----------	
	/* The concentrations of the soluble components in the floc is calculated.
	   It is calculated adding all the rates for each position */

	rR_o2 = 0;
	rR_nh4 = 0;
	rR_no2 = 0;
	rR_no3 = 0;
	rR_s = 0;
	for (int ii=0;ii<Nx;++ii){
		for (int jj=0;jj<Nx;++jj){
			rR_o2 += -Y_o2_het*(c_bac_HET[ii][jj]*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_o2[ii][jj]/(K_o2_het+c_o2[ii][jj]) + c_bac_HET[ii][jj]*eta_het*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_no2[ii][jj]/(K_no2_het+c_no2[ii][jj])*K_o2_het/(K_o2_het+c_o2[ii][jj]) + c_bac_HET[ii][jj]*eta_het*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_no3[ii][jj]/(K_no3_het+c_no3[ii][jj])*K_o2_het/(K_o2_het+c_o2[ii][jj])) - Y_o2_aob*(c_bac_AOB[ii][jj]*mu_max_aob*c_nh4[ii][jj]/(K_nh4_aob+c_nh4[ii][jj])*c_o2[ii][jj]/(K_o2_aob+c_o2[ii][jj])) - Y_o2_nob*(c_bac_NOB[ii][jj]*mu_max_nob*c_no2[ii][jj]/(K_no2_nob+c_no2[ii][jj])*c_o2[ii][jj]/(K_o2_nob+c_o2[ii][jj]));
			rR_nh4 += -Y_nh4_aob*(c_bac_AOB[ii][jj]*mu_max_aob*c_nh4[ii][jj]/(K_nh4_aob+c_nh4[ii][jj])*c_o2[ii][jj]/(K_o2_aob+c_o2[ii][jj]));
			rR_no2 += Y_no2_aob * (c_bac_AOB[ii][jj]*mu_max_aob*c_nh4[ii][jj]/(K_nh4_aob+c_nh4[ii][jj])*c_o2[ii][jj]/(K_o2_aob+c_o2[ii][jj])) - Y_no2_nob*(c_bac_NOB[ii][jj]*mu_max_nob*c_no2[ii][jj]/(K_no2_nob+c_no2[ii][jj])*c_o2[ii][jj]/(K_o2_nob+c_o2[ii][jj])) - Y_no2_het*c_bac_HET[ii][jj]*eta_het*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_no2[ii][jj]/(K_no2_het+c_no2[ii][jj])*K_o2_het/(K_o2_het+c_o2[ii][jj]); 
			rR_no3 += Y_no3_nob*(c_bac_NOB[ii][jj]*mu_max_nob*c_no2[ii][jj]/(K_no2_nob+c_no2[ii][jj])*c_o2[ii][jj]/(K_o2_nob+c_o2[ii][jj])) - Y_no3_het*c_bac_HET[ii][jj]*eta_het*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_no3[ii][jj]/(K_no3_het+c_no3[ii][jj])*K_o2_het/(K_o2_het+c_o2[ii][jj]);
			rR_s += -Y_s_het*(c_bac_HET[ii][jj]*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_o2[ii][jj]/(K_o2_het+c_o2[ii][jj]) + c_bac_HET[ii][jj]*eta_het*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_no2[ii][jj]/(K_no2_het+c_no2[ii][jj])*K_o2_het/(K_o2_het+c_o2[ii][jj]) + c_bac_HET[ii][jj]*eta_het*mu_max_het*c_s[ii][jj]/(K_s_het+c_s[ii][jj])*c_no3[ii][jj]/(K_no3_het+c_no3[ii][jj])*K_o2_het/(K_o2_het+c_o2[ii][jj]));
		}
	}
	rR_o2 *= nFlocs.back();
	rR_nh4 *= nFlocs.back();
	rR_no2 *= nFlocs.back();
	rR_no3 *= nFlocs.back();
	rR_s *= nFlocs.back();

}
// ----------------------------------------------------Solving PDEs-----------------------------------------------------				

void elliptic(double h, double k, double aa,double bb , std::vector< vector<double> > rMatrix,vector<vector<double> > &c_something,int type_soluble ){
	// INPUTS
	/* c_something (c_s,c_o2...) is needed in order to put the new concentrations in that vector */	
	/* 0 __ _h __ aa	         
	    |__|__|__|  
	   k|__|__|__|  
	    |__|__|__|  
	   bb     	        */

	// AIM
	/* this function tries to find the solutions of the BIVP'

       /  Uxx + Uyy = r(x,y)     U = U(x,y),   0<x<a ,  0<y<b 
      |   U(0,t) = g1(y)   U(a,y) = g2(y)              (boundary conditions) 
       \  U(x,0) = f1(x)   U(x,b) = f2(x)                (initial conditions) */

//   PROCEDURE

/*  The finite difference approximation is
        4*U(i,j) - ( U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1) ) = r(i,j)

  And the values in the border are already known, so we just one to
  calculate the values in the middle. 

    / Known   Known     Known     Known    Known \  
   |  Known   U(2,2)    U(2,3)    U(2,4)   Known  | 
   |  Known   U(3,2)    U(3,3)    U(3,4)   Known  | 
   |  Known   U(4,2)    U(4,3)    U(4,4)   Known  | 
   \  Known   Known     Known     Known    Known /     */
                                                                   
/* In order to work out the extra values we solve the following system 
							m											  <------ b -----> <---r--> 
	   /  4   -1    0   -1    0    0    0     0    0  \   / U(2,2) \    / U(1,2) + U(2,1) + r(2,2) \ 
	  |  -1    4   -1    0   -1    0    0     0    0   | |  U(2,3)  |   | U(1,3)          + r(2,3) | 
	  |   0   -1    4    0    0   -1    0     0    0   | |  U(2,4)  |   | U(1,4) + U(2,5) + r(2,4) | 
	  |  -1    0    0    4   -1    0   -1     0    0   | |  U(3,2)  |   |     U(3,1)      + r(3,1) | 
	n |   0   -1    0   -1    4   -1    0    -1    0   | |  U(3,3)  | = |       0           r(3,3) | 
	  |   0    0   -1    0   -1    4    0     0   -1   | |  U(3,4)  |   |     U(3,5)      + r(3,4) | 
	  |   0    0    0   -1    0    0    4    -1    0   | |  U(4,2)  |   | U(5,2)+  U(4,1) + r(4,2) | 
	  |   0    0    0    0   -1    0   -1     4   -1   | |  U(4,3)  |   | U(5,3)          + r(4,3) | 
	   \  0    0    0    0    0   -1    0    -1    4  /   \ U(4,4) /    \ U(5,4) + U(4,5) + r(4,4)/   */

//                MATRIX M                           VECTOR X              VECTOR B (= b + r)                   

//	# starting procedure for explicit algorithm 

	std::vector<double> I,J;
	// discretation of the x axis (are the columns in the matrix)
	for (double i=0;i<aa; i+=h){
		I.push_back(i);
	}
	// discretation of the y axis (are the rows in the matrix)  
	for (double j=0;j<bb; j+=k){
		J.push_back(j);
	}
	int m,n;
	m = I.size();
	n = J.size();
	
	// first step  , putting the initial conditions into the matrix
	gsl_matrix*U = gsl_matrix_calloc(n, m);
	
	for (int i=0;i<m;++i){
		gsl_matrix_set(U,0,i,boundary_values(type_soluble));
		gsl_matrix_set(U,n-1,i,boundary_values(type_soluble));
	}
	for (int j=0;j<n;++j){
		gsl_matrix_set(U,j,0,boundary_values(type_soluble));
		gsl_matrix_set(U,j,m-1,boundary_values(type_soluble));
	}
	
	// Creating MATRIX M

	
	gsl_matrix*M = gsl_matrix_calloc((n-2)*(m-2),(n-2)*(m-2));	
	
	/* Diagonal values*/
	for (int i=0; i<(n-2)*(m-2);++i){
		gsl_matrix_set(M,i,i,4);
	}
	/* Upper diagonal values */
	for (int i=0; i<(n-2)*(m-2)-1;++i){
			if ((i+1) % (m-2) == 0){
				gsl_matrix_set(M,i,i+1,0);
			}else{
				gsl_matrix_set(M,i,i+1,-1);
			}
	}
	/* Downer diagonal values */
	for (int i=1; i<(n-2)*(m-2);++i){
			if (i % (m-2) == 0){
				gsl_matrix_set(M,i,i-1,0);
			}else{
				gsl_matrix_set(M,i,i-1,-1);
			}
	}
	/* 4th upper diagonal */
	for (int i=0; i<(n-3)*(m-2);++i){
			gsl_matrix_set(M,i,i+(m-2),-1);
	}
	/* 4th downer diagonal */
	for (int i=m-2; i<(n-2)*(m-2);++i){
			gsl_matrix_set(M,i,i-(m-2),-1);
	}

	// Creating VECTOR r	(from the matrix)
	std::vector<double> r;
	for (int i=1;i<n-1;++i){
		for (int j=1; j<m-1;++j){
			r.push_back(rMatrix[i][j]);
		}		
	}

	// Creating VECTOR b
	std::vector<double> v_0(0,m-4);
	std::vector<double> b;

	for (int i=1; i<m-1;++i){
		b.push_back(gsl_matrix_get(U,0,i));
	}
	b[0] += gsl_matrix_get(U,1,0);	
	b[m-3] += gsl_matrix_get(U,1,m-1);
	
	for (int i=3; i<=n-2;++i){		
		b.push_back(gsl_matrix_get(U,i-1,0));
		b.insert(b.end(), v_0.begin(), v_0.end());
		b.push_back(gsl_matrix_get(U,i-1,m-1));
	}
	
	for (int i=1; i<m-1;++i){
		b.push_back(gsl_matrix_get(U,n-1,i));
	}	
	b[b.size()-(m-3)-1] += gsl_matrix_get(U,n-2,0);
	b[b.size()-1] += gsl_matrix_get(U,n-2,m-1);
	
	// Calculating VECTOR B
	std::vector<double> B;
	for(unsigned int i=0; i<r.size();++i){
		B.push_back(r[i] + b[i]);	
	}
	
	// Changing to gsl vector
	int len = B.size();
	gsl_vector*BB = gsl_vector_alloc(len);
	for(int i=0;i<len;++i){
		gsl_vector_set(BB,i,B[i]);
	}
	// Solving the system of linear equations

	gsl_vector*x = gsl_vector_alloc((n-2)*(m-2));	
	gsl_linalg_HH_solve(M,BB,x);
	
	for (int k=0;k<(n-2)*(m-2);++k){
		int row = floor(k/(m-2)) + 1;
		int col = k%(m-2) + 1 ;
		c_something[row][col] = gsl_vector_get(x,k);
	}
	for (unsigned j=0;j<c_something[0].size();++j){
		c_something[0][j] = boundary_values(type_soluble);
		c_something[c_something[0].size()-1][j] = boundary_values(type_soluble);
		c_something[j][0] = boundary_values(type_soluble);
		c_something[j][c_something[0].size()-1] = boundary_values(type_soluble);
		
	}
}
long double boundary_values(int type_soluble){
	// -------INPUTS----------
	/* "type_soluble" is the type of the soluble component that we want to calculate the values in the border of the grid*/
	// -------OUTPUTS---------
	/* the value in the bioreactor of that soluble component*/
	
	// -------EXPLANATION----------	
	/* It uses the last concentration in the bioreactor for that soluble component */

	return soluble_reactor[type_soluble].back();
}


