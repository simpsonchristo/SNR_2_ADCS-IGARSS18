// NewtonRaphson.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <numeric>
#include <functional>
#include <stdlib.h>
#include <Eigen/Dense>
#include <thread>
//#include <stan/math.hpp>


#include "stdafx.h"

using namespace std;
using namespace Eigen;

#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

VectorXd observinglocation();
VectorXd f(MatrixXd X, double t);

#endif /*NEWTONRAPHSON_H*/
