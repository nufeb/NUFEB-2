// ------------------------------------------------ Master Project ------------------------------------------ 
// -------------------------------------------- PDE solver elliptic ----------------------------------------
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



#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#define PI 3.1415926535897932384626433832795


using namespace std;

void elliptic(double h, double k, double aa,double bb , std::vector< vector<double> > rMatrix,vector<vector<double> > c_something );
int g1(double a);
int g2(double a);
int f1(double a);
int f2(double a);
int r_function(int i, int j);


int main(){
	
	vector<vector<double> > c_something(50,vector<double>(50,0));
	vector<vector<double> > rMatrix(50,vector<double>(50,0));	
	elliptic(1,1,50,50,rMatrix,c_something);
	
	cout << "\n 1 ha acabado" ;
		vector<vector<double> > c_something1(50,vector<double>(50,0));
	vector<vector<double> > rMatrix1(50,vector<double>(50,0));	
	elliptic(1,1,50,50,rMatrix1,c_something1);
	cout << "\n 2 ha acabado" ;
		vector<vector<double> > c_something2(50,vector<double>(50,0));
	vector<vector<double> > rMatrix2(50,vector<double>(50,0));	
	elliptic(1,1,50,50,rMatrix2,c_something2);
	cout << "\n 3 ha acabado" ;
		vector<vector<double> > c_something3(50,vector<double>(50,0));
	vector<vector<double> > rMatrix3(50,vector<double>(50,0));	
	elliptic(1,1,50,50,rMatrix3,c_something3);
	cout << "\n 4 ha acabado" ;
		vector<vector<double> > c_something4(50,vector<double>(50,0));
	vector<vector<double> > rMatrix4(50,vector<double>(50,0));	
	elliptic(1,1,50,50,rMatrix4,c_something4);
    
    //gsl_matrix_calloc  (3, 3);
    //cout << m ;
    
	return (0);
}

void elliptic(double h, double k, double aa,double bb , std::vector< vector<double> > rMatrix,vector<vector<double> > c_something ){
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
	cout << "\n (n-2)*(m-2)" << (n-2)*(m-2);
	cout << "\n m:" << m;
	cout << "\n n:" << n;
	// first step  , putting the initial conditions into the matrix
	gsl_matrix*U = gsl_matrix_calloc(n, m);
	
	for (int i=0;i<m;++i){
		gsl_matrix_set(U,0,i,g1(I[i]));
		gsl_matrix_set(U,n-1,i,g2(I[i]));
	}
	for (int j=0;j<n;++j){
		gsl_matrix_set(U,j,0,f1(J[j]));
		gsl_matrix_set(U,j,m-1,f2(J[j]));
	}
	
	// Creating MATRIX M
/*	
	zero_matrix<int> m1 ((n-2)*(m-2));
	
	matrix<int> M ((n-2)*(m-2),(n-2)*(m-2));
	M = m1;*/
	
	gsl_matrix*M = gsl_matrix_calloc((n-2)*(m-2),(n-2)*(m-2));	
	
	/* Diagonal values*/
	for (int i=0; i<(n-2)*(m-2);++i){
		//M(i,i) = 4;
		gsl_matrix_set(M,i,i,4);
	}
	/* Upper diagonal values */
	for (int i=0; i<(n-2)*(m-2)-1;++i){
			if ((i+1) % (m-2) == 0){
				//M(i,i+1) = 0;
				gsl_matrix_set(M,i,i+1,0);
			}else{
				//M(i,i+1) = -1;
				gsl_matrix_set(M,i,i+1,-1);
			}
	}
	/* Downer diagonal values */
	for (int i=1; i<(n-2)*(m-2);++i){
			if (i % (m-2) == 0){
				//M(i,i-1) = 0;
				gsl_matrix_set(M,i,i-1,0);
			}else{
				//M(i,i-1) = -1;
				gsl_matrix_set(M,i,i-1,-1);
			}
	}
	/* 4th upper diagonal */
	for (int i=0; i<(n-3)*(m-2);++i){
			//M(i,i+(m-2)) = -1;
			gsl_matrix_set(M,i,i+(m-2),-1);
	}
	/* 4th downer diagonal */
	for (int i=m-2; i<(n-2)*(m-2);++i){
			//M(i,i-(m-2)) = -1;
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
	//std::vector<double> x((n-2)*(m-2));
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
		int row = floor(k/(m-2));
		int col = k%(m-2);
		c_something[row][col] = gsl_vector_get(x,k);
	}

	cout << "\n everything is fine";
}
int g1(double a){
	return 0;
}
int g2(double a){
	return 0;
}
int f1(double a){
	return 0;
}
int f2(double a){
	return 0;
}

int r_function(int i, int j){	
	return 3*i + j;	
}
