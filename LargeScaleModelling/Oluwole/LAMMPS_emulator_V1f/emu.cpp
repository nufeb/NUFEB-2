//to compile g++ -O3 -I/usr/include/eigen3/ emu.cpp -o emu
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
//#include <eigen3/Eigen/Dense>
#include <Eigen/Core>
using namespace std;
using namespace Eigen;

//this function computes distance correlation
MatrixXd mycor(const MatrixXd& X1,const MatrixXd& X2,const VectorXd& betai)
{
MatrixXd ones = MatrixXd::Ones(1,X2.rows());
MatrixXd ones2 = MatrixXd::Ones(1,X1.rows());
MatrixXd pdf = betai.asDiagonal();//turn vector to digaonal matrix
MatrixXd R1 = X1*pdf*X1.transpose();
MatrixXd R2 = X2*pdf*X2.transpose();
MatrixXd S1 = R1.diagonal()*ones;
MatrixXd S2 = R2.diagonal()*ones2;
MatrixXd a1 = X1*pdf*X2.transpose(); //t(tx1)%*%pdf%*%tx2
MatrixXd a2 = X2*pdf*X1.transpose();
MatrixXd res = (a1.transpose()+a2-S1.transpose()-S2);
ArrayXXd res2 = res.array();
return(res2.exp());
}
///////////////////////
//main function
int main()
{
	int nrows=4130;
	int ncols=4;
	int ncols1=5; //newline
	int ncols2=9;
	VectorXd betai(9);
betai << 0.003847799,0.001000000,0.001000000,0.077725734,0.001000000,6.862597279,4.714029223,2.095145620,4.58666442;
VectorXd beta2(9);  //newline
beta2 <<  0.16060842, 0.02898964, 0.03102735, 0.57128240, 0.00100000, 12.12281649, 0.00100000, 16.34735376, 2.11941945; //newline
   //read training input dat matrix X
MatrixXd X = MatrixXd::Zero(nrows,ncols2);
ifstream fin1 ("./X.txt");
if (fin1.is_open())
{
for (int row = 0; row < nrows; row++)
for (int col = 0; col < ncols2; col++)
{
float item;
fin1 >> item;
X(row, col) = item;
}
fin1.close();
}
//read in output matrix Y
MatrixXd Y = MatrixXd::Zero(nrows,ncols);
ifstream fin2 ("./Y.txt");
if (fin2.is_open())
{
for (int row = 0; row < nrows; row++)
for (int col = 0; col < ncols; col++)
{
float item;
fin2 >> item;
Y(row, col) = item;
}
fin2.close();
}

//read in new data point for prediction
MatrixXd newdata = MatrixXd::Zero(1,ncols2);
ifstream fin3 ("./input/input.txt");
if (fin3.is_open())
{
for (int row = 0; row < 1; row++)
for (int col = 0; col < ncols2; col++)
{
float item;
fin3 >> item;
newdata(row,col) = item;
}
fin3.close();
}

//read in Ainv matrix
MatrixXd Ainv = MatrixXd::Zero(nrows,nrows);
ifstream fin4 ("./Ainv2.txt");
if (fin4.is_open())
{
for (int row = 0; row < nrows; row++)
for (int col = 0; col < nrows; col++)
{
float item;
fin4 >> item;
Ainv(row, col) = item;
}
fin4.close();
}
//read in betahat matrix
MatrixXd betahat = MatrixXd::Zero(ncols2,ncols);
ifstream fin6 ("./betahat.txt");
if (fin6.is_open())
{
for (int row = 0; row < ncols2; row++)
for (int col = 0; col < ncols; col++)
{
float item;
fin6 >> item;
betahat(row, col) = item;
}
fin6.close();
}

//read in new betahat2 matrix
MatrixXd betahat2 = MatrixXd::Zero(ncols2,ncols1);
ifstream fin8 ("./betahat2.txt");
if (fin8.is_open())
{
for (int row = 0; row < ncols2; row++)
for (int col = 0; col < ncols1; col++)
{
float item;
fin8 >> item;
betahat2(row, col) = item;
}
fin8.close();
}
//read in new outputmatrix Y2 for nutrients
MatrixXd Y2 = MatrixXd::Zero(nrows,ncols1);
ifstream fin9 ("./Y2.txt");
if (fin9.is_open())
{
for (int row = 0; row < nrows; row++)
for (int col = 0; col < ncols1; col++)
{
float item;
fin9 >> item;
Y2(row, col) = item;
}
fin9.close();
}

//read in Ainv matrix for emulator 2
MatrixXd Ainv2 = MatrixXd::Zero(nrows,nrows);
ifstream fin10 ("./Ainv22.txt");
if (fin10.is_open())
{
for (int row = 0; row < nrows; row++)
for (int col = 0; col < nrows; col++)
{
float item;
fin10 >> item;
Ainv2(row, col) = item;
}
fin10.close();
}

//normalization data
MatrixXd s1(1,ncols);
MatrixXd s2(1,ncols);
MatrixXd t1(1,ncols1); //newline
MatrixXd t2(1,ncols1); //newline
s1 << 2.467027e-15,2.339083e+01,8.968750e-10,6.864521e-11; 
s2 << 2.166465e-15,2.079999e+01,4.336633e-10,4.304684e-10; 
t1 << -4.503880e-08,-8.699600e-09,-6.134320e-09,-8.181856e-08,-4.000000e-12; //newline 
t2 << 1.94704e-08,2.02620e-08,8.79828e-09,1.66498e-08,6.00000e-12; 

//emulator 2 for the rate of change of "Biomass concentration","number of particle","Biofilm height","Biofilm surface"
MatrixXd H = X.array().log();//log-tranformed training inputs
MatrixXd H0 = newdata.array().log();//log of test input
int m=H.cols();
int n=Y.rows();
int n2=H0.rows();
MatrixXd A00, A01;
A01.noalias() = mycor(H0,H,betai);//###cross correlation//dont evaluate yet
A00.noalias() = mycor(H0,H0,betai);//##test point correlation
MatrixXd mu_star = H0*betahat+A01.transpose()*Ainv.selfadjointView<Upper>()*(Y-H*betahat);
MatrixXd rate2 = (mu_star.array()*s2.array())+s1.array();

//emulator 1 for the rate of change of nutrient concentrations ("nh3","no2","no3","o2","co2")
 A01.noalias() = mycor(H0,H,beta2);//###cross correlation//dont evaluate yet
 A00.noalias() = mycor(H0,H0,beta2);//##test point correlation
mu_star = H0*betahat2+A01.transpose()*Ainv2.selfadjointView<Upper>()*(Y2-H*betahat2);
MatrixXd rate1 = (mu_star.array()*(t2.array()-t1.array()))+t1.array();
//output results
ofstream a_file ("./output/Rate.txt");
a_file << rate1 << ' '<< rate2;
a_file.close();
//cout << rate1 << rate2 << endl;
}


