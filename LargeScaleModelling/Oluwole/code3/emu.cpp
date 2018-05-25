#include <iostream>
#include <iomanip>
#include <fstream>
//#include <armadillo>
using namespace std;
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
//#include <gsl/GaussianCF.h>

// double pred(double newdata)
int main ()
{
int i,j,nout=4;
//int ntime=100;
//int nsim=20; //##number of Monte carlo simulation
//double tau=2.0;
int id2=10;	
int nrow=2950;
int ncol1=4;
int ncol2=9;
FILE *f1,*f2,*f3,*f4;
gsl_matrix *X,*Y,*coef,*input,*H;	
double mij[nrow][ncol1];
double olu;
//double mycor;

//read in X and Y matrix	
 X=gsl_matrix_alloc(nrow,ncol2);
 Y=gsl_matrix_alloc(nrow,ncol1);
 coef=gsl_matrix_alloc(10,5);
 input=gsl_matrix_alloc(1,9);
 H=gsl_matrix_alloc(ncol2,nrow);
 
  f1=fopen("X.txt","r");
  gsl_matrix_fscanf(f1,X);
  fclose(f1);
  
  f2=fopen("Y.txt","r");
  gsl_matrix_fscanf(f2,Y);
  fclose(f2);	
  // gsl_matrix_fprintf(stdout,Y,"%f");
  
  f3=fopen("coef.txt","r");
  gsl_matrix_fscanf(f3,coef);
  fclose(f3);	
 // gsl_matrix_fprintf(stdout,coef,"%f");
  
  f4=fopen("input.txt","r");
  gsl_matrix_fscanf(f4,input);
  fclose(f4);	
 // gsl_matrix_fprintf(stdout,input,"%f");
 
  for (i = 0; i < nrow; i++)
    for (j = 0; j < ncol1; j++)
      {
         mij = gsl_matrix_get (Y, i, j);
       // double aij = gsl_matrix_get (a, i, j);
       // if (mij != aij) k++;
      }
 olu << mij;
//double scale2[] ={0.201935506,0.006313825,0.000010000,0.177466408,0.000010000,6.577416780,7.721354529,6.594165598,4.081750274,0.004609465};
//double theta[]={0.201935506,0.006313825,0.000010000,0.177466408,0.000010000,6.577416780,7.721354529,6.594165598,4.081750274}; //scale2[-*(id2-1)];	
//double scale2b=0.004609465;

//void gsl_matrix_set_all (gsl_matrix * H, double X);
gsl_matrix * newdata=input;
gsl_matrix *  X0=newdata;
int gsl_matrix_transpose_memcpy(gsl_matrix Xinv,gsl_matrix *X);
double gsl_matrix_memcpy(gsl_matrix *H,gsl_matrix X);

//int gsl_matrix_memcpy (gsl_matrix * dest, const gsl_matrix * src)
gsl_matrix * H0=X0;
int m=ncol2; // H->size2;
int n= nrow; //Y->size1;
int n2= 1;  //newdata->size2;
// double A01=mycor(X0,X,theta);

// print out matrix
for (i=0;i<ncol2;i++)
    {
      for (j=0;j<nrow;j++)
	{
	  printf("%f ",gsl_matrix_get(H,i,j));
	}
      printf("\n");
    }
  printf("\n");
  
  for (i=0;i<1;i++)
    {
      for (j=0;j<9;j++)
	{
	  printf("%f ",gsl_matrix_get(H0,i,j));
	}
      printf("\n");
    }
  printf("\n");
  //* release memory */
  gsl_matrix_free(X);
  gsl_matrix_free(Y);
  
// cout << nsim << tau << nout <<endl;
// cout << ntime << id2 << endl; 
//for (i=0;i<id2;i++)
//cout<< theta[i] << endl;
//for (i=0;i<id2;i++)
//cout << scale2[i] << endl;
//cout << H[i][j] << endl;
cout << olu << endl;
 cout << input << endl;
return 0;
}
//



