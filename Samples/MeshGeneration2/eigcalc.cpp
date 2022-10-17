#include "../../GFD/Discrete/Sparse.hpp"
#include  <vector>
#include <iostream>
#include <Eigen/Dense>
 
using Eigen::MatrixXd;
using namespace gfd;
using namespace std;

int main() {
    Sparse<double> u;
    vector<double> a = {1,2,3,4};
    Buffer<double> vals(a);

    MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    std::cout << m << std::endl;
}