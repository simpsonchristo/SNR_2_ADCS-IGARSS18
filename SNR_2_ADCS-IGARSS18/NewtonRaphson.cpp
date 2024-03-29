#include "NewtonRaphson.h"
//Observing location
VectorXd observinglocation() {
	VectorXd Xs(4);
	Xs[0] = 1.0;//x
	Xs[1] = 1.0;//z
	Xs[2] = 0.0;//vx
	Xs[3] = 0.0;//vz
	return Xs;
}
//range equation
VectorXd f(VectorXd X, double t) {
	VectorXd A(X.rows());
	VectorXd Xs = observinglocation();
	for (int i = 0; i < X.rows(); i++)
	{
		A[i] = pow((pow(X[0] - Xs[0] + (X[2] * t), 2)) +
			(pow(X[1] - Xs[1] + (X[3] * t) - (0.5*X[4] * pow(t, 2)), 2))
			, .5);
	}

	return A;
}
VectorXd df(VectorXd X, double t) {
	VectorXd B(X.rows());
	VectorXd Xs = observinglocation();
	double rho;
	for (int i = 0; i < X.rows(); i++)
	{
		rho =pow((pow(X[0] - Xs[0] + (X[2] * t), 2)) +
			(pow(X[1] - Xs[1] + (X[3] * t) - (0.5*X[4] * pow(t, 2)), 2))
			, .5);

		B[0] = (X[0] - Xs[0] + (X[2] * t)) / rho;
		B[1] = (X[1] - Xs[1] + (X[3] * t) - (0.5*X[4] * pow(t, 2))) / rho;
		B[2] = 0;
		B[3] = 0;
		B[4] = ((-0.5*pow(t,2))*(X[1] - Xs[1] + (X[3] * t) - (0.5*X[4] * pow(t, 2)))) / rho;		
	}
	return B;
}
