#include <Eigen/Dense>
#include <cmath>
#include <math.h>
#include <vector>
#include <functional>
#include <iostream>
//#include <stan\math\rev\core.hpp>


using namespace std;
// using namespace Eigen;

#ifndef ABM_INTEGRATE_H
#define ABM_INTEGRATE_H

template <class Tc>
class ABM_Integrate {
    
    public:
        // Typedefs
        typedef Eigen::Matrix<Tc,Eigen::Dynamic,Eigen::Dynamic> MatrixXp;
        typedef Eigen::Matrix<Tc,Eigen::Dynamic,1> VectorXp;
        typedef Eigen::Matrix<Tc,3,3> Matrix3p;
        typedef Eigen::Matrix<Tc,3,1,0,3,1> Vector3p;
        
        // CONSTRUCTORS
        ABM_Integrate(int integ_order);
        ABM_Integrate();
        
        // Setup
        void set_func(std::function <VectorXp (VectorXp, Tc)> *fun);
        std::function <VectorXp (VectorXp, Tc)> fcn;
        void set_order(int o);
        void set_init_step_max(Tc i);
        
        // Stepper setup
        void set_stepper_dt(Tc h);
        void fill_stepper_y(MatrixXp y);
        void fill_stepper_t(VectorXp t);
        
        // The function to integrate
        VectorXp func(VectorXp y, Tc t);
        
        // Calcs
        int integ(VectorXp &y0, Tc t0, Tc tf, Tc dt);
        int step();
        
        // Intermediate - Integ
        void init_integrator();
        VectorXp f(int phase,int step);
        int mem_index(int m);
        
        // Intermediate - Step
        void latest_stepper_y(VectorXp y_n);
        void init_stepper();
        VectorXp f_step(int step);
        
        // Access
        MatrixXp get_y();
        VectorXp get_t();
        VectorXp get_final_y();
        Tc get_final_t();
        MatrixXp get_error();
        
        // Export
        /*void write_y();
        void write_t();*/
        
    private:
        Eigen::MatrixXi AB;
        Eigen::MatrixXi AM;
        
        MatrixXp store_y;
        VectorXp store_t;
        MatrixXp key_points;
        VectorXp key_time;
        
        MatrixXp mem_dydt;
        std::vector<int> recalc;
        MatrixXp mem_dydt_next;
        std::vector<int> recalc_next;
        
        int order;
        int dim;
        Tc dt;
        int index;
        Tc init_step_max;
        int m;
        int div;
        std::vector<Tc> h;
        int mem_index_shift;
        
        Tc step_dt;

        MatrixXp error;
        
};


#endif /* ABM_INTEGRATE_H */