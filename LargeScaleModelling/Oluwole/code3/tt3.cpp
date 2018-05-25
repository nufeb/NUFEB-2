//to compile g++ -I/usr/include/eigen3/ tt.cpp -o tt
#include <istream>
#include <string>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Sparse>
//#include <eigen3>
#include <fstream>
#include <eigen3/unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>
using namespace std;
using namespace Eigen;

//this function computes distance correlation
MatrixXf mycor(const MatrixXf& X1,const MatrixXf& X2,const VectorXf& betai)
{
MatrixXf ones = MatrixXf::Ones(1,X2.rows());
MatrixXf ones2 = MatrixXf::Ones(1,X1.rows());
MatrixXf pdf = betai.asDiagonal();//turn vector to digaonal matrix
MatrixXf R1 = X1*pdf*X1.transpose();
MatrixXf R2 = X2*pdf*X2.transpose();
MatrixXf S1 = R1.diagonal()*ones;
MatrixXf S2 = R2.diagonal()*ones2;
MatrixXf a1 = X1*pdf*X2.transpose(); //t(tx1)%*%pdf%*%tx2
MatrixXf a2 = X2*pdf*X1.transpose();
MatrixXf res = (a1.transpose()+a2-S1.transpose()-S2);
ArrayXXf res2 = res.array();
return(res2.exp());
}
///////////////////////
//main function
int main()
{
	int nrows=4130;
	int ncols=4;
	int ncols2=9;
	double tau0=0.001000;
	double tau=0.5;
VectorXf betai(9);
betai << 0.003847799,0.001000000,0.001000000,0.077725734,0.001000000,6.862597279,4.714029223,2.095145620,4.58666442;
//read training input dat matrix X
MatrixXf X = MatrixXf::Zero(nrows,ncols2);
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
MatrixXf Y = MatrixXf::Zero(nrows,ncols);
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
MatrixXf newdata = MatrixXf::Zero(1,ncols2);
ifstream fin3 ("./input/input.txt");
if (fin3.is_open())
{
for (int row = 0; row < 1; row++)
for (int col = 0; col < ncols2; col++)
{
float item;
fin3 >> item;
newdata(row, col) = item;
}
fin3.close();
}

//read in Ainv matrix
MatrixXf Ainv = MatrixXf::Zero(nrows,nrows);
ifstream fin4 ("./Ainv.txt");
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

//
//MatrixXf output=mycor(X1,X2,betai);
//preprocessing data
MatrixXf s1(1,4);
MatrixXf s2(1,4);
s2 << 2.166465e-15,2.079999e+01,4.336633e-10,4.304684e-10; 
s1 << 2.467027e-15,2.339083e+01,8.968750e-10,6.864521e-11; 
//Y=scale(Y,scale=s2,center=s1);
MatrixXf H = X.array().log();
MatrixXf H0 = newdata.array().log();
//H=X
int m=H.cols();
int n=Y.rows();
int n2=H0.rows();
MatrixXf A01 = mycor(H0,H,betai);//###cross correlation
MatrixXf A00 = mycor(H0,H0,betai);//##test point correlation
//MatrixXf A = mycor(H,H,betai);
MatrixXf temp = MatrixXf::Identity(A00.rows(), A00.cols());
MatrixXf A000 = A00+((tau0*tau)*temp);
MatrixXf iOmega = (H.transpose()*Ainv*H).inverse();
MatrixXf betahat= iOmega*(H.transpose()*Ainv*Y);
MatrixXf mu_star = H0*betahat+A01.transpose()*Ainv*(Y-H*betahat);
MatrixXf r1 = H0-(A01.transpose()*Ainv*H);
MatrixXf c_star = A000-(A01.transpose()*Ainv*A01)+(((r1)*(iOmega)*r1.transpose()));
MatrixXf Sigma = ((Y-H*betahat).transpose()*Ainv*(Y-H*betahat))/(n-m);
MatrixXf rate = (mu_star.array()*s2.array())+s1.array();
//output results
ofstream a_file ("./output/Rate2.txt");
a_file << rate;
a_file.close();
cout << rate << endl;
cout << mu_star << endl;

}
