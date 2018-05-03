#include "ABM_Integrate.h"

// Typedefs
template <class Tc>
using MatrixXp = typename ABM_Integrate<Tc>::MatrixXp;
template <class Tc>
using Matrix3p = typename ABM_Integrate<Tc>::Matrix3p;
template <class Tc>
using VectorXp = typename ABM_Integrate<Tc>::VectorXp;
template <class Tc>
using Vector3p = typename ABM_Integrate<Tc>::Vector3p;

////////////////////////////////////////////////////////////////////////////////
// Constructors
template <class Tc>
ABM_Integrate<Tc>::ABM_Integrate(int integ_order){
    this->set_order(integ_order);
    
    init_step_max = 0.001;
    
    AB = Eigen::MatrixXi::Zero(10,11);
    AB.row(0) << 1,1,0,0,0,0,0,0,0,0,0;
    AB.row(1) << 3,-1,2,0,0,0,0,0,0,0,0;
    AB.row(2) << 23,-16,5,12,0,0,0,0,0,0,0;
    AB.row(3) << 55,-59,37,-9,24,0,0,0,0,0,0;
    AB.row(4) << 1901,-2770,2616,-1274,251,720,0,0,0,0,0;
    AB.row(5) << 4277,-7923,9982,-7298,2877,-475,1440,0,0,0,0;
    AB.row(6) << 198721,-447288,705549,-688256,407139,-134472,19087,60480,0,0,0;
    AB.row(7) << 434241,-1152169,2183877,-2664477,2102243,-1041723,295767,-36799,120960,0,0;
    AB.row(8) << 14097247,-43125206,95476786,-139855262,137968480,-91172642,38833486,-9664106,1070017,3628800,0;
    AB.row(9) << 30277247,-104995189,265932680,-454661776,538363838,-444772162,252618224,-94307320,20884811,-2082753,7257600;
    
    AM = Eigen::MatrixXi::Zero(10,11);
    AM.row(0) << 1,1,0,0,0,0,0,0,0,0,0;
    AM.row(1) << 1,1,2,0,0,0,0,0,0,0,0;
    AM.row(2) << 5,8,-1,12,0,0,0,0,0,0,0;
    AM.row(3) << 9,19,-5,1,24,0,0,0,0,0,0;
    AM.row(4) << 251,646,-264,106,-19,720,0,0,0,0,0;
    AM.row(5) << 475,1427,-798,482,-173,27,1440,0,0,0,0;
    AM.row(6) << 19087,65112,-46461,37504,-20211,6312,-863,60480,0,0,0;
    AM.row(7) << 36799,139849,-121797,123133,-88547,41499,-11351,1375,120960,0,0;
    AM.row(8) << 1070017,4467094,-4604594,5595358,-5033120,3146338,-1291214,312874,-33953,3628800,0;
    AM.row(9) << 2082753,9449717,-11271304,16002320,-17283646,13510082,-7394032,2687864,-583435,57281,7257600;


}
template <class Tc>
ABM_Integrate<Tc>::ABM_Integrate(){
    init_step_max = 0.001;

    
    AB = Eigen::MatrixXi::Zero(10,11);
    AB.row(0) << 1,1,0,0,0,0,0,0,0,0,0;
    AB.row(1) << 3,-1,2,0,0,0,0,0,0,0,0;
    AB.row(2) << 23,-16,5,12,0,0,0,0,0,0,0;
    AB.row(3) << 55,-59,37,-9,24,0,0,0,0,0,0;
    AB.row(4) << 1901,-2770,2616,-1274,251,720,0,0,0,0,0;
    AB.row(5) << 4277,-7923,9982,-7298,2877,-475,1440,0,0,0,0;
    AB.row(6) << 198721,-447288,705549,-688256,407139,-134472,19087,60480,0,0,0;
    AB.row(7) << 434241,-1152169,2183877,-2664477,2102243,-1041723,295767,-36799,120960,0,0;
    AB.row(8) << 14097247,-43125206,95476786,-139855262,137968480,-91172642,38833486,-9664106,1070017,3628800,0;
    AB.row(9) << 30277247,-104995189,265932680,-454661776,538363838,-444772162,252618224,-94307320,20884811,-2082753,7257600;
    
    AM = Eigen::MatrixXi::Zero(10,11);
    AM.row(0) << 1,1,0,0,0,0,0,0,0,0,0;
    AM.row(1) << 1,1,2,0,0,0,0,0,0,0,0;
    AM.row(2) << 5,8,-1,12,0,0,0,0,0,0,0;
    AM.row(3) << 9,19,-5,1,24,0,0,0,0,0,0;
    AM.row(4) << 251,646,-264,106,-19,720,0,0,0,0,0;
    AM.row(5) << 475,1427,-798,482,-173,27,1440,0,0,0,0;
    AM.row(6) << 19087,65112,-46461,37504,-20211,6312,-863,60480,0,0,0;
    AM.row(7) << 36799,139849,-121797,123133,-88547,41499,-11351,1375,120960,0,0;
    AM.row(8) << 1070017,4467094,-4604594,5595358,-5033120,3146338,-1291214,312874,-33953,3628800,0;
    AM.row(9) << 2082753,9449717,-11271304,16002320,-17283646,13510082,-7394032,2687864,-583435,57281,7257600;

}

// Set
template <class Tc>
void ABM_Integrate<Tc>::set_func(std::function <VectorXp (VectorXp, Tc)> *fun){
    fcn = *fun;
}
template <class Tc>
void ABM_Integrate<Tc>::set_order(int integ_order){
    order = integ_order;
}
template <class Tc>
void ABM_Integrate<Tc>::set_init_step_max(Tc i){
    init_step_max = i;
}

// Stepper setup
template <class Tc>
void ABM_Integrate<Tc>::set_stepper_dt(Tc h){
    step_dt = h;
}
template <class Tc>
void ABM_Integrate<Tc>::fill_stepper_y(MatrixXp y){
    store_y.resizeLike(y);
    store_y = y;
}
template <class Tc>
void ABM_Integrate<Tc>::fill_stepper_t(VectorXp t){
    store_t.resizeLike(t);
    store_t = t;
}

// Access
template <class Tc>
Eigen::Matrix<Tc,Eigen::Dynamic,Eigen::Dynamic> ABM_Integrate<Tc>::get_y(){
    return store_y;
}
template <class Tc>
Eigen::Matrix<Tc,Eigen::Dynamic,1> ABM_Integrate<Tc>::get_t(){
    return store_t;
}
template <class Tc>
Eigen::Matrix<Tc,Eigen::Dynamic,1> ABM_Integrate<Tc>::get_final_y(){
    return store_y.col(store_y.cols()-1);
}
template <class Tc>
Tc ABM_Integrate<Tc>::get_final_t(){
    return store_t[store_t.size()-1];
}
template <class Tc>
Eigen::Matrix<Tc,Eigen::Dynamic,Eigen::Dynamic> ABM_Integrate<Tc>::get_error(){
    return error;
}

// The function to integrate
template <class Tc>
Eigen::Matrix<Tc,Eigen::Dynamic,1> ABM_Integrate<Tc>::func(VectorXp y, Tc t){
//    cout << "******** Function Called at Time = " << t << " ********" << endl;
//    cout << "Y = " << y << endl;
    bool breakp = 1;
//    cout << "Function Returned" << endl;
    return fcn(y,t);
}

template <class Tc>
int ABM_Integrate<Tc>::mem_index(int m){
    return fmod(m + mem_index_shift, order+1);
}


////////////////////////////////////////////////////////////////////////////////
////////////        INTEGRATOR        //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <class Tc>
void ABM_Integrate<Tc>::init_integrator(){
    mem_dydt.resize(dim,order+1);
    recalc.resize(order+1);
    std::fill(recalc.begin(), recalc.end(), 1);
    mem_dydt_next.resize(dim,order+1);
    recalc_next.resize(order+1,1);
    key_points.resize(dim,order);
    mem_index_shift = 0;
    error.resize(dim,1);
}

template <class Tc>
Eigen::Matrix<Tc,Eigen::Dynamic,1> ABM_Integrate<Tc>::f(int phase, int step){
    // The purpose of this function is to speed up the integration by holding on
    //     to the last few dydt evaluations. 
    // Primarily, we need an indexing function to lookup what we've stored and 
    //     so we know when to overwrite.
    // Also should separate dt increments and dt_prelim increments.
    // Size of mem_dydt is (dim,order);
    
    // Phase 1, store sequentially, step = index, no looping yet
    // Phase 2, standard integration, but with dt_prelim. Extract dt increments, 
    //     indexer will need to loop
    // Phase 3, use extracted dt increments and looping
    // Phase 4, standard integration procedure with dt, looping
    
    using std::round;
    
    switch (phase) {
        case 1:
            // In case 1, step = index
            if(recalc[step] == 1){
                mem_dydt.col(step) = func(store_y.col(step),store_t(step));
                recalc[step] = 0;
            
                if(fmod(store_t(step) + h[0]/pow(order,2) - store_t(0), h[1]) < init_step_max/pow(order,2)){
                    int idx = round( (store_t(step)-store_t(0) )/h[1] );
                    mem_dydt_next.col(idx) = mem_dydt.col(step);
//                    cout << "Stored time " << store_t(step) << " in mem_dydt_next.col " << idx << endl;
                }
            }
//            cout << "Returning mem_dydt.col(step) = " << mem_dydt.col(step) << endl;
            return mem_dydt.col(step);
            break;
            
            
        case 2:
            // In case 2, step = index
//            cout << "Calling Case 2 with step = " << step << endl;
//            cout << "Time at step = " << store_t(step) << endl;
//            cout << "Mem_index(step) = " << mem_index(step) << endl;
            
            if(recalc[mem_index(step)] == 1){
                mem_dydt.col(mem_index(step)) = func(store_y.col(step),store_t(step));
                recalc[mem_index(step)] = 0;
                recalc[mem_index(step+1)] = 1;
                            
                if(fmod(store_t(step) + h[0]/pow(order,2) - store_t(0), dt/pow(order,div-(m+1)) ) < init_step_max/pow(order,2)){
                    int idx = round( (store_t(step)-store_t(0))/(dt/pow(order,div-(m+1))) );
//                    Tc idx_truth = store_t(step)/h[1];
//                    cout << "Step = " << step << endl;
//                    Tc store_t_at_step = store_t(step);
//                    cout << "store_t at step = " << store_t(step) << endl;
//                    cout << "h[1] = " << h[1] << endl;
                    mem_dydt_next.col(idx) = mem_dydt.col(mem_index(step));
                    recalc_next[idx] = 0;
//                    cout << "Stored at time " << store_t(step) << " in mem_dydt_next.col " << idx << endl;
                }
            }
//            cout << "f Called at index = " << step << " returns mem_dydt at " << mem_index(step) << endl;

            return mem_dydt.col(mem_index(step));
            break;
            
            
        case 3:
            // In case 3, step = # of steps in this h value
//            cout << "Activated f() at Phase 3 with step = " << step << endl;
//            cout << "MemIndex(step) = " << mem_index(step) << endl;
//            cout << "Step " << step << " corresponds with y value " << key_points.col(step) << endl;

            if(recalc[mem_index(step)] == 1){
                mem_dydt.col(mem_index(step)) = func(key_points.col(step),key_time(step));
                recalc[mem_index(step)] = 0;
                recalc[mem_index(step+1)] = 1;
                
                    // Visualize recalc
//                    cout << "Recalc = ";
//                    for(int p=0; p<recalc.size(); p++){
//                        cout << recalc[p];
//                    }
//                    cout << endl;

                
//                cout << "key_time(step) = " << key_time(step) << endl;
//                cout << "fmod inputs = " << key_time(step) + h[0]/pow(order,1) << " and " << h[m+1] << endl;
//                cout << "fmod result = " << fmod(key_time(step) + h[0]/pow(order,1), h[m+1]) << endl;
//                cout << "tolerance = " << 2*init_step_max/order << endl;
                if(fmod(key_time(step) + h[0]/pow(order,2) - store_t(0), dt/pow(order,div-(m+1)) ) < init_step_max/pow(order,2)){
//                    cout << "Attempting to store for extraction at " << mem_index(step) << endl;
                    int idx = round((key_time(step)-store_t(0))/(dt/pow(order,div-(m+1))));
//                    cout << "IDX = " << idx << endl;
//                    cout << "mem_index(step) = " << mem_index(step) << endl;
                    mem_dydt_next.col(idx) = mem_dydt.col(mem_index(step));
                    recalc_next[idx] = 0;
//                    cout << "Stored time " << key_time(step) << " in mem_dydt_next.col " << idx << endl;
                }
            }
//            cout << "f called at step = " << step << " returns mem_dydt at " << mem_index(step) << endl;
            
            return mem_dydt.col(mem_index(step));
            break;
            
            
        case 4:
            // step = index
            if(recalc[mem_index(step)] == 1){
                mem_dydt.col(mem_index(step)) = func(store_y.col(step),store_t(step));
                recalc[mem_index(step)] = 0;
                recalc[mem_index(step+1)] = 1;
            }
            
            return mem_dydt.col(mem_index(step));
            break;
            
        default:
            cout << "Incorrect Phase Input" << endl;
    }
}
//////////////////////////////////////////////////////////////////////////////////
/////////////////////    STAN::MATH::VAR SPECIFIC INSTANCE    ////////////////////
//////////////////////////////////////////////////////////////////////////////////
//template <>
//Eigen::Matrix<stan::math::var,Eigen::Dynamic,1> ABM_Integrate<stan::math::var>::f(int phase, int step){
//    // The purpose of this function is to speed up the integration by holding on
//    //     to the last few dydt evaluations. 
//    // Primarily, we need an indexing function to lookup what we've stored and 
//    //     so we know when to overwrite.
//    // Also should separate dt increments and dt_prelim increments.
//    // Size of mem_dydt is (dim,order);
//    
//    // Phase 1, store sequentially, step = index, no looping yet
//    // Phase 2, standard integration, but with dt_prelim. Extract dt increments, 
//    //     indexer will need to loop
//    // Phase 3, use extracted dt increments and looping
//    // Phase 4, standard integration procedure with dt, looping
//    
//    using std::round;
//    
//    switch (phase) {
//        case 1:
//            // In case 1, step = index
//            if(recalc[step] == 1){
//                mem_dydt.col(step) = func(store_y.col(step),store_t(step));
//                recalc[step] = 0;
//            
//                if(fmod(store_t(step) + h[0]/pow(order,2) - store_t(0), h[1]) < init_step_max/pow(order,2)){
//                    int idx = round( (store_t(step).val()-store_t(0).val() )/h[1].val() );
//                    mem_dydt_next.col(idx) = mem_dydt.col(step);
////                    cout << "Stored time " << store_t(step) << " in mem_dydt_next.col " << idx << endl;
//                }
//            }
////            cout << "Returning mem_dydt.col(step) = " << mem_dydt.col(step) << endl;
//            return mem_dydt.col(step);
//            break;
//            
//            
//        case 2:
//            // In case 2, step = index
////            cout << "Calling Case 2 with step = " << step << endl;
////            cout << "Time at step = " << store_t(step) << endl;
////            cout << "Mem_index(step) = " << mem_index(step) << endl;
//            
//            if(recalc[mem_index(step)] == 1){
//                mem_dydt.col(mem_index(step)) = func(store_y.col(step),store_t(step));
//                recalc[mem_index(step)] = 0;
//                recalc[mem_index(step+1)] = 1;
//                            
//                if(fmod(store_t(step) + h[0]/pow(order,2) - store_t(0), dt/pow(order,div-(m+1)) ) < init_step_max/pow(order,2)){
//                    int idx = round( (store_t(step).val()-store_t(0).val())/(dt.val()/pow(order,div-(m+1))) );
////                    Tc idx_truth = store_t(step)/h[1];
////                    cout << "Step = " << step << endl;
////                    Tc store_t_at_step = store_t(step);
////                    cout << "store_t at step = " << store_t(step) << endl;
////                    cout << "h[1] = " << h[1] << endl;
//                    mem_dydt_next.col(idx) = mem_dydt.col(mem_index(step));
//                    recalc_next[idx] = 0;
////                    cout << "Stored at time " << store_t(step) << " in mem_dydt_next.col " << idx << endl;
//                }
//            }
////            cout << "f Called at index = " << step << " returns mem_dydt at " << mem_index(step) << endl;
//
//            return mem_dydt.col(mem_index(step));
//            break;
//            
//            
//        case 3:
//            // In case 3, step = # of steps in this h value
////            cout << "Activated f() at Phase 3 with step = " << step << endl;
////            cout << "MemIndex(step) = " << mem_index(step) << endl;
////            cout << "Step " << step << " corresponds with y value " << key_points.col(step) << endl;
//
//            if(recalc[mem_index(step)] == 1){
//                mem_dydt.col(mem_index(step)) = func(key_points.col(step),key_time(step));
//                recalc[mem_index(step)] = 0;
//                recalc[mem_index(step+1)] = 1;
//                
//                    // Visualize recalc
////                    cout << "Recalc = ";
////                    for(int p=0; p<recalc.size(); p++){
////                        cout << recalc[p];
////                    }
////                    cout << endl;
//
//                
////                cout << "key_time(step) = " << key_time(step) << endl;
////                cout << "fmod inputs = " << key_time(step) + h[0]/pow(order,1) << " and " << h[m+1] << endl;
////                cout << "fmod result = " << fmod(key_time(step) + h[0]/pow(order,1), h[m+1]) << endl;
////                cout << "tolerance = " << 2*init_step_max/order << endl;
//                if(fmod(key_time(step) + h[0]/pow(order,2) - store_t(0), dt/pow(order,div-(m+1)) ) < init_step_max/pow(order,2)){
////                    cout << "Attempting to store for extraction at " << mem_index(step) << endl;
//                    int idx = round((key_time(step).val()-store_t(0).val())/(dt.val()/pow(order,div-(m+1))));
////                    cout << "IDX = " << idx << endl;
////                    cout << "mem_index(step) = " << mem_index(step) << endl;
//                    mem_dydt_next.col(idx) = mem_dydt.col(mem_index(step));
//                    recalc_next[idx] = 0;
////                    cout << "Stored time " << key_time(step) << " in mem_dydt_next.col " << idx << endl;
//                }
//            }
////            cout << "f called at step = " << step << " returns mem_dydt at " << mem_index(step) << endl;
//            
//            return mem_dydt.col(mem_index(step));
//            break;
//            
//            
//        case 4:
//            // step = index
//            if(recalc[mem_index(step)] == 1){
//                mem_dydt.col(mem_index(step)) = func(store_y.col(step),store_t(step));
//                recalc[mem_index(step)] = 0;
//                recalc[mem_index(step+1)] = 1;
//            }
//            
//            return mem_dydt.col(mem_index(step));
//            break;
//            
//        default:
//            cout << "Incorrect Phase Input" << endl;
//    }
//}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Operations
template <class Tc>
int ABM_Integrate<Tc>::integ(VectorXp &y0, Tc t0, Tc tf, Tc dt_in){
    // Integrates function f from t0 to tf with step size dt starting at y0
    // Stores results in store_y and store_t for access later
    
    // Initialize based on size
    dt = dt_in;
    dim = y0.size();
    store_y.resize(dim,1);
    store_y.col(0) = y0;
    store_t.resize(1);
    store_t(0) = t0;
    index = 0;
    init_integrator();
//    init_step_max = 1;
    
    // Confirmation
//    cout << "Integrating:" << endl;
//    cout << y0 << endl;
//    cout << "from" << endl;
//    cout << t0 << " to " << tf << endl;
//    cout << "With step size " << endl;
//    cout << dt << endl;
    
    // Preliminary -- Build up base to proceed
//    Tc dt_prelim = dt/order;
//    Tc t_endPrelim = dt*(order-1) + t0;
    // How far do we need to subdivide in order to start out accurately?
    div = 1;
    while(dt/(pow(order,div)) > init_step_max){ div++; }
    for(int sub=0; sub<div; sub++){
        h.push_back( dt/(pow(order,div-sub)) );
    }
    h.push_back(dt);
    m = 0;
    

    // Phase 1 - Build from nothing - Increase order as possible
//    cout << "Phase 1" << endl;
    for(int maxO=1; maxO<order; maxO++){
        VectorXp sum = VectorXp::Zero(dim);
        recalc[maxO] = 1;
        // Predict
        for(int i=0; i<maxO; i++){
            sum += AB(maxO-1,i)*f(1,index-i);
        }
//        cout << "Sum = " << sum << endl;;
//        cout << "Stored Previous State = " << store_y.col(index) << endl;
//        cout << "Test Mult and Div = " << dt_prelim*sum/AB(maxO-1,maxO) << endl;
//        cout << "Test Summing = " << store_y.col(index) + dt_prelim*sum/AB(maxO-1,maxO) << endl;
//        cout << "What does store_y look like = " << store_y << endl;
        
        store_y.conservativeResize(Eigen::NoChange, index+2);
//        cout << "Now store_y looks like = " << store_y << endl;
        store_y.col(index+1) = store_y.col(index) + h[0]*sum/AB(maxO-1,maxO);
        store_t.conservativeResize(index+2);
        store_t(index+1) = store_t(index) + h[0];
        error.conservativeResize(Eigen::NoChange, index+2);
        error.col(index+1) = store_y.col(index+1);
        
//        cout << "Prediction = " << store_y.col(index+1) << endl;
        
//        cout << "Now store_y looks like = " << store_y << endl;
//        cout << "Test Storage = " << store_y.col(index+1) << endl;

        // Correct
        sum = VectorXp::Zero(dim);
        for(int i=0; i<maxO; i++){
//            cout << "Corrector debug" << endl;
//            cout << "AM: " << AM(maxO,i) << endl;
//            cout << "Current state: " << store_y.col(index-i+1) << endl;
//            cout << "Current time: " << store_t(index-i+1) << endl;
            sum += AM(maxO-1,i)*f(1,index-i+1);
        }
        store_y.col(index+1) = store_y.col(index) + h[0]*sum/AM(maxO-1,maxO);
        error.col(index+1) -= store_y.col(index+1);
        // This is to update the prediction dydt
        recalc[mem_index(index+1)] = 1;
        
//        cout << "Correction = " << store_y.col(index+1) << endl;
        
        index++;
//        cout << "Most Recent Result at Index = " << index << endl << store_y.col(index) << endl;
        if(abs(store_t(index) - tf) < init_step_max/pow(order,2)) return 0;
    }
// We have now created enough points to integrate at the given order    
    
//    cout << "Recalc Between Phase 1 and 2 = ";
//    for(int p=0; p<recalc.size(); p++){
//        cout << recalc[p];
//    }
//    cout << endl;
    
        // Phase 2 - Build from above - Still small dt
    recalc[order] = 1;
//    store_y.conservativeResize(Eigen::NoChange, order*(order-1)+1);
//    store_t.conservativeResize(order*(order-1)+1);
//    cout << endl << "Phase 2" << endl;
    while( index < order*(order-1) ){
//        recalc_prelim[mem_index(index+1)] = 1;
        VectorXp sum = VectorXp::Zero(dim);
        // Predict
        for(int i=order-1; i>=0; i--){
            sum += AB(order-1,i)*f(2,index-i);
        }
        store_y.conservativeResize(Eigen::NoChange, index+2);
        store_y.col(index+1) = store_y.col(index) + h[0]*sum/AB(order-1,order);
        store_t.conservativeResize(index+2);
        store_t(index+1) = store_t(index) + h[0];
        error.conservativeResize(Eigen::NoChange, index+2);
        error.col(index+1) = store_y.col(index+1);
    
        // Correct
        sum = VectorXp::Zero(dim);
        for(int i=order-1; i>=0; i--){
            sum += AM(order-1,i)*f(2,index-i+1);
        }
        store_y.col(index+1) = store_y.col(index) + h[0]*sum/AM(order-1,order);
        error.col(index+1) -= store_y.col(index+1);
        // This is to update the prediction dydt
        recalc[mem_index(index+1)] = 1;
    
        index++;
//        cout << "Most Recent Result at Index = " << index << endl << store_y.col(index) << endl;
        if(abs(store_t(index) - tf) < init_step_max/pow(order,2)) return 0;
    }
// We have now created enough points to extract out 'order' points.

    
    while(m < div){ // While the step size is less than dt
        // Phase 3 - Pull out states at h[m+1] points
//        cout << endl << "Phase 3 -- Iteration " << m+1 << endl;
//        cout << "Phase 3.1" << endl;
        key_points.resize(dim,order);
        key_time.resize(order);
        key_points = MatrixXp::Zero(dim,order);
        key_time = VectorXp::Zero(order);
        int k=0;
        int j=0;
    //    cout << "Got Past key_point initialization" << endl;
//        cout << "Extracting..." << endl;
        while(k<order){
            // Extract first 'order' states at the given dt intervals
//            cout << "j value of " << j << " gives fmod of " << fmod(store_t(j),h[m+1]) << " at Time = " << store_t(j) << endl;
            if(fmod(store_t(j) + h[0]/pow(order,2) - t0, h[m+1]) < init_step_max/order){
                key_points.col(k) = store_y.col(j);
                key_time(k) = store_t(j);
                k++;
//                cout << "Extracted " << k << " out of " << order << " states." << endl;
//                cout << "Extracted at t = " << store_t(j) << endl;
//                cout << "Y-value = " << store_y.col(j) << endl;
            }
            j++;
            if(j>pow(order,m+2)-order){
//                cout << "Failed to Extract Key Points" << endl;
//                cout << "Extracted " << k << " out of " << order << " Key Points" << endl;
            }
        }
//        cout << "Got past key_point extraction" << endl;
// The points for higher-step-size integration have been removed.
        
        // Prepare for next higher step-size
        mem_dydt = mem_dydt_next;
        mem_index_shift = 0;
        std::fill(recalc.begin(), recalc.end(), 0);
        recalc[recalc.size()-1] = 1;
        recalc[recalc.size()-2] = 1;
        m++;
        
        // Operate on key_points only
        
//        cout << "Phase 3.2" << endl;
        for(int n=0; n<(order); n++){
    //        recalc[mem_index(index)] = 1;
            // Step forward on the extracted values
            VectorXp sum = VectorXp::Zero(dim);
            // Predict
//            cout << endl << "3) PREDICT" << endl;
            for(int i=order-1; i>=0; i--){
//                cout << "Calling AB at " << i << endl;
//                cout << "Calling f at " << n+order-i-1 << endl;
                
                sum += AB(order-1,i)*f(3,(n+order)-i-1);
            }
        
            key_points.conservativeResize(Eigen::NoChange, n+order+1);
            key_points.col(n+order) = key_points.col(n+order-1) + h[m]*sum/AB(order-1,order);
            key_time.conservativeResize(n+order+1);
            key_time(n+order) = key_time(n+order-1) + h[m];
            error.conservativeResize(Eigen::NoChange, index+2);
            error.col(index+1) = key_points.col(n+order);
            
            // Correct
            sum = VectorXp::Zero(dim);
//            cout << endl << "3) CORRECT" << endl;
            for(int i=order-1; i>=0; i--){
//                cout << "Calling AM at " << i << endl;
//                cout << "Calling f at " << n+order-i << endl;
                sum += AM(order-1,i)*f(3,(n+order)-i);
            }
            key_points.col(n+order) = key_points.col(n+order-1) + h[m]*sum/AM(order-1,order);
            error.col(index+1) -= key_points.col(n+order);
//            cout << "Updated prediction in key_points" << endl;
                                
            // Store in store_y for future
            store_y.conservativeResize(Eigen::NoChange, index+2);
            store_y.col(index+1) = key_points.col(n+order);
            store_t.conservativeResize(index+2);
            store_t(index+1) = key_time(n+order);
            // This is to update the prediction dydt
            recalc[mem_index(n+order)] = 1;
            
            index++;
//            cout << "Size of key_points = " << key_points.size() << endl;
//            cout << "Most Recent Result at Index = " << index << endl << "T = " << store_t(index) << endl << "Y = " << store_y.col(index) << endl;
            if(abs(store_t(index) - tf) < init_step_max/pow(order,2)) return 0;
        }
        
//        cout << "Most Recent Timestep = " << store_t(store_t.size()-1) << endl;
//        cout << "Mem_Index(Index) between Phase 3.2 and 3.3 = " << mem_index(index) << endl;
//        cout << "Mem_Index of most recent 'step' = " << mem_index(2*order-1) << endl;
        
        // Transition off of key_points
        mem_index_shift = mem_index(2*order-1) - mem_index(index);
//        cout << "Shifted mem_index(index) = " << mem_index(index) << endl;
        
        // Operate on store_y only
        // Ending condition is that next key_points is full, that is, key_points.cols() = order
        // Can probobly use the same f() storage as Phase 2
//        cout << "Phase 3.3" << endl;
//        cout << "Conditional Checks:" << endl;
//        cout << "Store_t(i) - h(0)/pow(o,2) - t0 = " << (store_t(index) - h[0]/pow(order,2) - t0) << endl;
//        cout << "(dt/pow(order,div-(m+1)))*(order-1) - h(m) = " << ((dt/pow(order,div-(m+1)))*(order-1) - h[m]) << endl;
        
        while( (store_t(index) - h[0]/pow(order,2) - t0) < ((dt/pow(order,div-(m+1)))*(order-1) - h[m]) ){
            VectorXp sum = VectorXp::Zero(dim);
            // Predict
//            cout << "PREDICT" << endl;
            for(int i=order-1; i>=0; i--){
//                cout << "Calling f at " << index-i << endl;
                sum += AB(order-1,i)*f(2,index-i);
            }
            store_y.conservativeResize(Eigen::NoChange, index+2);
            store_y.col(index+1) = store_y.col(index) + h[m]*sum/AB(order-1,order);
            store_t.conservativeResize(index+2);
            store_t(index+1) = store_t(index) + h[m];
            error.conservativeResize(Eigen::NoChange, index+2);
            error.col(index+1) = store_y.col(index+1);
            
//            cout << "Prediction = " << store_y.col(index+1) << endl;
            
            // Correct
//            cout << "CORRECT" << endl;
            sum = VectorXp::Zero(dim);
            for(int i=order-1; i>=0; i--){
//                cout << "Calling f at " << index-i+1 << endl;
                sum += AM(order-1,i)*f(2,index-i+1);
            }
            store_y.col(index+1) = store_y.col(index) + h[m]*sum/AM(order-1,order);
            error.col(index+1) -= store_y.col(index+1);
            // This is to update the prediction dydt
            recalc[mem_index(index+1)] = 1;
            
//            cout << "Correction = " << store_y.col(index+1) << endl;
            
            index++;
//            cout << "Most Recent Result at Index = " << index << endl << "T = " << store_t(index) << endl << "Y = " << store_y.col(index) << endl;
            if(abs(store_t(index) - tf) < init_step_max/pow(order,2)) return 0;
        }
    }
    
    
    // Phase 4 - Run Integration - passed dt
//    cout << "Phase 4" << endl;
    while(store_t(index)<tf){
//        recalc[mem_index(index)] = 1;
        VectorXp sum = VectorXp::Zero(dim);
        // Predict
        for(int i=order-1; i>=0; i--){
            sum += AB(order-1,i)*f(4,index-i);
        }
        store_y.conservativeResize(Eigen::NoChange, index+2);
        store_y.col(index+1) = store_y.col(index) + dt*sum/AB(order-1,order);
        store_t.conservativeResize(index+2);
        store_t(index+1) = store_t(index) + dt;
        error.conservativeResize(Eigen::NoChange, index+2);
        error.col(index+1) = store_y.col(index+1);

        // Correct
        sum = VectorXp::Zero(dim);
        for(int i=order-1; i>=0; i--){
            sum += AM(order-1,i)*f(4,index-i+1);
        }
        store_y.col(index+1) = store_y.col(index) + dt*sum/AM(order-1,order);
        error.col(index+1) -= store_y.col(index+1);
        // This is to update the prediction dydt
        recalc[mem_index(index+1)] = 1;
        
        index++;
//        cout << "Most Recent Result at Index = " << index << endl << "T = " << store_t(index) << endl << "Y = " << store_y.col(index) << endl;
        if(abs(store_t(index) - tf) < init_step_max/pow(order,2)) return 0;
    }
    
    return 1;
}



////////////////////////////////////////////////////////////////////////////////
////////////        STEPPER        /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <class Tc>
void ABM_Integrate<Tc>::latest_stepper_y(VectorXp y_n){
    // Used to update latest prediction with the true value
    store_y.col(index) = y_n;
    recalc[mem_index(index)] = 1;
}

template <class Tc>
void ABM_Integrate<Tc>::init_stepper(){
    dim = store_y.rows();
    mem_dydt.resize(dim,order+1);
    mem_index_shift = 0;
    recalc.resize(order+1);
    std::fill(recalc.begin(), recalc.end(), 1);
    index = store_y.cols()-1;
    error.resize(dim,1);
}

template <class Tc>
int ABM_Integrate<Tc>::step(){
    
    cout << "In Stepper" << endl;
    
    if(store_y.cols() < order){
        cout << "Insufficient fill for specified order" << endl;
        return 0;
    }
    cout << "Current size of store_y = " << store_y.cols() << endl;
    cout << "Current index = " << index << endl;
    
    VectorXp sum = VectorXp::Zero(dim);
    // Predict
    // cout << "Predict \n";
    for(int i=order-1; i>=0; i--){
        // cout << i << endl;
        sum += AB(order-1,i)*f_step(index-i);
    }
    recalc[mem_index(index+1)] = 1;
    
    // Visualize recalc
/*    
*    cout << "Recalc = ";
*    for(int p=0; p<recalc.size(); p++){
*        cout << recalc[p];
*    }
*    cout << endl;
*/
    
    store_y.conservativeResize(Eigen::NoChange, index+2);
    store_y.col(index+1) = store_y.col(index) + step_dt*sum/AB(order-1,order);
    store_t.conservativeResize(index+2);
    store_t(index+1) = store_t(index) + step_dt;
    error.conservativeResize(Eigen::NoChange, index+2);
    error.col(index+1) = store_y.col(index+1);
    
    // Correct
    sum = VectorXp::Zero(dim);
    for(int i=order-1; i>=0; i--){
        sum += AM(order-1,i)*f_step(index-i+1);
    }
    store_y.col(index+1) = store_y.col(index) + step_dt*sum/AM(order-1,order);
    error.col(index+1) -= store_y.col(index+1);
    //This is to update the prediction dydt
    recalc[mem_index(index+1)] = 1;
    
    index++;
    
    cout << "Most Recent Result:" << endl;
    cout << "Step = " << index << endl;
    cout << "T = " << store_t(index) << endl;
    cout << "Y = " << endl << store_y.col(index) << endl << endl;
    
}

template <class Tc>
Eigen::Matrix<Tc,Eigen::Dynamic,1> ABM_Integrate<Tc>::f_step(int step){
    // Step = index
    if(recalc[mem_index(step)] == 1){
        mem_dydt.col(mem_index(step)) = func(store_y.col(step), store_t(step));
        recalc[mem_index(step)] = 0;
//        recalc[mem_index(step+1)] = 1;
        
    }
    
    return mem_dydt.col(mem_index(step));
    
}

////////////////////////////////////////////////////////////////////////////////
/// Instantiations
template class ABM_Integrate<double>;
//template class ABM_Integrate<float>;
template class ABM_Integrate<long double>;
//template class ABM_Integrate<stan::math::var>;