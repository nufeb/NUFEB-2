// C++ 
// Wastewater treatment
// Arturo Alvarez-Arenas
// Newcastle University


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


using namespace std;


// .........GENERAL FUNCTIONS........

long length(float vector[])
{
	return sizeof(vector)/sizeof(*a)
}


//   0) ...........INITIALIZATE FUNCTIONS.................

		void calculate_soluble_floc(float soluble_floc[], const double Diffusion_coefficients[], double rates_soluble);
		void calculate_soluble_reactor(float Si_in[],soluble_reactor,HRT,rates,nFlocs);
		void attachment_single_cell(vector<long double> &Bac_m,vector<long double> &Bac_e_d,vector<int> &Bac_s,vector<int> &Bac_c, vector<long double> &Bac_r,double Rmax);
		void attachment_group_cell(B);
		void dettachment_group_cell(B);
		void removing_islands(B);
		void decay(B,rates_decay);
		void growth(B,rates_growth);
		void division(B);
		void moving_cells(B);


//   1) .............SETTING PARAMETERS............

//       1.1) Bioreactor and separator sizes & others

			const short nGrowth = 1000;
			const short nSingleAttach = 1;
			const short nGroupAttach = 2;
			const short nGroupDetach = 2*nGroupAttach+1;
			
			
			
//       1.2) Properties initial bacterium

			const long double mass_het = 1e-16;
			const long double mass_aob = 1e-16;
			const long double mass_nob = 1e-16;
			const long double mass_eps = 2.6e-17;   % kg for grid 1e-6 % for having eps density 25 kg/m3
			const long double mass_inert = 1.1e-16; % kg for grid 1e-6 % for having eps density 105 kg/m3
			
			long number_bacterium = 5 /* Total number of bacterium, it will be modified during this process*/

//       1.3) Diffusion coefficients in water 
			
			const double D_o2  = 2e-9*3600*24;     /* m^2/d */
			const double D_nh4 = 1.4e-9*3600*24;   /* m^2/d */
			const double D_no2 = 1.2e-9*3600*24;   /* m^2/d */
			const double D_no3 = 1.2e-9*3600*24;   /* m^2/d */
			const double D_s   = 0.5e-9*3600*24;   /* m^2/d */
						
			const double Diffusion_Coefficients[] = {D_o2,D_nh4,D_no2,D_no3,D_s}; /* stored in an array */
			 	
//       1.4) Kinetic parameters 

//       1.5) Decay constants 

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
						
			float Si_in[] = {So2_in,Snh4_in,Sno2_in,Sno3_in,Ss_in};	/* stored in an array*/


//       1.6) Parameters single attach process

						       /* pHET pAOB pNOB */ 
			const float prob[] = {0.6, 0.2, 0.2};  /* probabilities for the new cell in single attach process*/
			const long double mass_new_cell[] = {mass_het,mass_aob,mass_nob}; 
			const long double eps_new_cell[]={mass_eps,0,0};	
			
//       1.7) Initializating vectors with cells information

			long double bac_x[] = {};  /* x's position of each cell*/
			long double bac_y[] = {};  /* y's position of each cell*/
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

			vector<double> Bac_theta; /* This angles are in the interval [-pi pi]*/
			vector<double> Bac_R;		

//   2) .................MAIN LOOP.................

int main ()
{
	/* Calculating some Values*/

	/* Recalculated in each step Bac_r and Bac_ra*/
	Bac_r.clear();
	Bac_ra.clear();
	for (unsigned int i=0; i<Bac_x.size();++i){
		Bac_r.push_back(sqrt((Bac_m[i]/bac_rho + Bac_e_d[i]/bac_rhoe)/PI/bac_h));
		Bac_ra.push_back(sqrt((Bac_m[i]/bac_rho)/PI/bac_h));
	}	
	
	/* Recalculated in each step Bac_theta and Bac_R*/ 
	Bac_theta.clear();
	Bac_R.clear();
	for (unsigned int i=0; i<Bac_x.size();++i){				/* changing from cartesian to polar coordinates*/	
		Bac_theta.push_back(atan2(Bac_y[i],Bac_x[i]));
		Bac_R.push_back (sqrt(pow(Bac_x[i],2) + pow(Bac_y[i],2)));
	}
		
	for (short iterGrowth=0; iterGrowth<nGrowth; iterGrowth++){

	
// ---------------------------------------------MICRO-SCALE MODEL--------------------------------------------------------

//       2.1) Soluble concentrations floc 

			//double rates_soluble[] ={}  PENSAR COMO NARICES CALCULAR ESTO, SI CON UNA FUNCION O A LAS BRAVAS AQUI MISMO
			

			calculate_soluble_floc (soluble_floc, Diffusion_coefficients, rates_soluble);

//       2.2) Saving the actual mass of the floc
			
			double mFlocOld = boost::accumulate(Bac_m, 0) + boost::accumulate(Bac_e_d, 0);
			
			/* Moving the cells from the center of the floc to the left-down corner*/
			




			
			/* calculating Rmax*/	
			vector<double>::const_iterator it, RmaxVector;
			long double Rmax;
			RmaxVector = max_element(Bac_R.begin(),Bac_R.end());
			Rmax = *RmaxVector;

		

	

//       2.2) Single attach process
		


			if ( iterGrowth % nSingleAttach == 0 )
			{
				attachment_single_cell(Bac_m,Bac_e_d,Bac_s,Bac_c,Bac_r,Rmax);
							
			}
			
//       2.3) Group attach process
	
			if ( iterGrowth % nGroupAttach == 0 )
			{
				attachment_group_cell(B)
			
			}


//       2.4) Group dettach process

			if ( iterGrowth % nGroupDetach == 0 )
			{
				detachment_group_cell(B)
			
			}

//       2.5) Removing islands process

			if ( iterGrowth % nGroupDetach == 1 )
			{
				removing_islands(B)
			
			}

//       2.6) Number of flocs
			
			float mFlocNew
			double mFlocNew = boost::accumulate(Bac_m, 0) + boost::accumulate(Bac_e_d, 0);
			
			// nFlocs(end+1) = nFlocs(end)(1- mFlocOld*dt/(mFlocNew*SDT)) put this formula with c++ syntax
			

//       2.7) Rates growth & decay
	
			//rates_decay = []
			

//       2.8) Decay process

			growth_decay();

//       2.9) Recalculating Bac_r and Bac_ra before division.

			
			
//       2.10) Soluble concentrations bioreactor

			calculate_soluble_reactor(Si_in,soluble_reactor,HRT,rates,nFlocs)


//       2.11) Division and EPS excretion process

			division(B)

//       2.12) Repositioning cells 

			moving_cells(B)

	}

}


	

//   3) ...................PLOTS...................















//   4) ..............DEFINED FUNCTIONS.........

		void calculate_soluble_floc(float soluble_floc[], const double Diffusion_coefficients[], double rates_soluble)
		{
			// define function, as well as the inputs, and change void if neccesary
		}


		void calculate_soluble_reactor(Si_in,soluble_reactor,HRT,rates,nFlocs)
		{
			// define function, as well as the inputs, and change void if neccesary
		}
		
				
		void attachment_single_cell(vector<long double> &Bac_m,vector<long double> &Bac_e_d,vector<int> &Bac_s,vector<int> &Bac_c, vector<long double> &Bac_r,double Rmax)
		{

			
			/* generate a random position for the new cell*/
			
			
				double fi_atta =  ((double) rand() / (RAND_MAX)) * 2 * PI;
				long double R_atta = Rmax + 5*Bac_r[0]; 
				

				/* attracting the new cell to the floc*/
				for(int j=0;j<1;--j){
			
					R_atta = R_atta - Bac_r[0];
					long double x_atta = R_atta * cos(fi_atta);
					long double y_atta = R_atta * sin(fi_atta);


					for (unsigned int i=0; i<Bac_x.size();++i){					
						long double dist = sqrt(pow((Bac_x[i] - x_atta),2) + pow((Bac_y[i] - y_atta),2));
						if ((dist - 3*Bac_r[i])<0){ /* means new cell is close enough to that cell */
							Bac_x.push_back(x_atta);
							Bac_y.push_back(y_atta);
							j=2;	                  /* breaks the initial for*/
							break;                    /* breaks this for */
						}
					}
						 
				}

			/* choosing type of the new cell*/
			
			double r  = ((double) rand() / (RAND_MAX)); /* generating a random number between 0,1*/

			
			for(int j=1; j<4; ++j){
				valarray<float> prob1 (prob,j);
				if (r< prob1.sum()){			
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
			
			
			Bac_r.push_back(sqrt((Bac_m.back()/bac_rho + Bac_e_d.back()/bac_rhoe)/PI/bac_h));
			Bac_ra.push_back(sqrt((Bac_m.back()/bac_rho)/PI/bac_h));
			Bac_theta.push_back (atan(Bac_y.back()/Bac_x.back()));
			Bac_R.push_back (sqrt(pow(Bac_x.back(),2) + pow(Bac_y.back(),2)));
			
			++number_bacterium;
		}
			
			
			
		
		
		void attachment_group_cell(B)
		{
			// define function, as well as the inputs, and change void if neccesary
		}

		void dettachment_group_cell(B)
		{
			// define function, as well as the inputs, and change void if neccesary
		}

		void removing_islands(B)
		{
			// define function, as well as the inputs, and change void if neccesary
		}	
		
		void decay(B,rates_decay)
		{
			// define function, as well as the inputs, and change void if neccesary
		}								

		void growth(B,rates_growth)
		{
			// define function, as well as the inputs, and change void if neccesary
		}	
		
		void division(B)
		{
			// define function, as well as the inputs, and change void if neccesary
		}	
		
		void moving_cells(B)
		{
			// define function, as well as the inputs, and change void if neccesary
		}			
