#include <istream>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <iostream>
#include <eigen>

Eigen::MatrixXcd mat(rows, cols);
mat = Eigen::Map<Eigen::VectorXd> (csvData.data(), rows, cols).cast<complex<double> >().transpose

MatrixXf X = MatrixXf::Zero(nrows,ncols);
ifstream fin ("./data.txt");

if (fin.is_open())
{
    for (int row = 0; row < nrows; row++)
        for (int col = 0; col < ncols; col++)
        {
            float item = 0.0;
            fin >> item;
            X(row, col) = item;
        }
    fin.close();
}
cout << "X = " << endl << X << endl;
