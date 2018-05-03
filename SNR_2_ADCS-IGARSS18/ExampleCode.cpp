#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <numeric>
#include <functional>
#include <stdlib.h>
#include <Eigen/Dense>
#include <thread>
#include <stan/math.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <math.h>
#include <boost/random.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/timer.hpp>

#include "../../ABM_Integrate.h"

using namespace std;
using namespace Eigen;

// Define ODE
VectorXd gravity(VectorXd x, double t){
    VectorXd xDot;
    xDot[0] = x[3];
    xDot[1] = x[4];
    xDot[2] = x[5];
    
    Vector3d acc = acceleration_func(x.head(3));
    
    xDot[3] = acc[0];
    xDot[4] = acc[1];
    xDot[5] = acc[2];
    
    return xDot;
};

int main() {

    VectorXd X0;
    double ti;
    double tf;
    double dt;
    
    
    std::function <VectorXd (VectorXd,double) > odefun = std::bind( &gravity, std::placeholders::_1, std::placeholders::_2);
    ABM_Integrate<double> integrator;
    integrator.set_order(10);
    integrator.set_func(odefun);
    int exit_status = integrator.integ(X0, ti, tf, dt);
   
    return 0;
}

