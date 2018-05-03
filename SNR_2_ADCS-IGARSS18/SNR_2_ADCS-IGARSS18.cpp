// SNR_2_ADCS-IGARSS18.cpp : Defines the entry point for the console application.
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
//#include <boost/math/tools/promotion.hpp>
//#include <boost/math/special_functions/legendre.hpp>
#include <math.h>
//#include <boost/random.hpp>
//#include <boost/multiprecision/float128.hpp>
//#include <boost/timer.hpp>
#include <fstream>

#include "stdafx.h"
#include "ABM_Integrate.h"
#include "NewtonRaphson.h"

using namespace Eigen;
using namespace std;

ofstream outputfile;
std::string timefile = "t.csv";
std::string propfile = "Xstar.csv";
std::string estfile = "Xe.csv";


//Define ODE
VectorXd acceleration_func(VectorXd x) {
	VectorXd acc(2);
	//[m/s^2]
	acc[0] = 0;//x
	acc[1] = x[0];//z
	return acc;
};
VectorXd gravity(VectorXd x, double t) {
	//Position
	VectorXd xDot(5);
	xDot[0] = x[2];
	xDot[1] = x[3];

	VectorXd acc = acceleration_func(x.tail(1));

	//Velocity
	xDot[2] = acc[0];
	xDot[3] = acc[1];
	xDot[4] = 0;

	return xDot;
};

int main()
{
	VectorXd X0(5);
	double ti = 0;
	double tf = 4;
	double dt = 1;

	//Assign intial state
	X0[0] = 1.5;//x
	X0[1] = 10.0;//z
	X0[2] = 2.2;//vx
	X0[3] = 0.5;//vz
	X0[4] = -0.3;//g

	std::function <VectorXd (VectorXd, double)> odefun = std::bind( &gravity, std::placeholders::_1, std::placeholders::_2);
	ABM_Integrate<double> integrator;
	integrator.set_order(10);
	integrator.set_func(&odefun);
	integrator.set_init_step_max(1.0);
	int exit_status = integrator.integ(X0, ti, tf, dt);
	MatrixXd Xf = integrator.get_y().transpose();
	VectorXd store_t = integrator.get_t();

	outputfile.open(timefile);
	//write header
	outputfile << "t (s)" << std::endl;
	//write time
	for (int i = 0; i < store_t.rows(); i++)
	{
		outputfile << store_t[i] << ", " << std::endl;
	}
	outputfile.close();

	outputfile.open(propfile);
	//write header
	outputfile << "x" << "," << "y" << "," << "u" << "," << "v" << "," << "g" << std::endl;
	//write state vector
	for (int i = 0; i < Xf.rows(); i++)
	{
		for (int j = 0; j < Xf.cols(); j++)
		{
			outputfile << Xf(i, j) << ", ";
		}
		outputfile << std::endl;
	}
	outputfile.close();
	

	return 0;	
}

