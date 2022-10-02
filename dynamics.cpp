#include "dynamics.h"
//====================================================================
// Boost stuff
//====================================================================

#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>
//using namespace std;
using namespace boost::numeric::odeint;

typedef std::vector< double > state_type;

typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;


class rbdlToBoost {

public:
    rbdlToBoost(Model* model) : model(model) {
        q = VectorNd::Zero(model->dof_count);
        qd = VectorNd::Zero(model->dof_count);
        qdd = VectorNd::Zero(model->dof_count);
        tau = VectorNd::Zero(model->dof_count);

   }

   //3c. Boost uses this 'operator()' function to evaluate the state
   //    derivative of the pendulum.
   void operator() (const state_type &x,
                    state_type &dxdt,
                    const double t){

       //3d. Here we split out q (generalized positions) and qd
       //    (generalized velocities) from the x (state vector)
       //q
       int j = 0;
       for(int i=0; i<model->dof_count; i++){
           q[i] = (double)x[j];
           j++;
       }

       //qd
       for(int i=0; i<model->dof_count; i++){
           qd[i] = (double)x[j];
           j++;
       }

       //3e. Here we set the applied generalized forces to zero
    //    for(int i=0; i<model->dof_count; i++){
    //        tau[i] = 0;
    //    }

       //3f. RBDL's ForwardDynamics function is used to evaluate
       //    qdd (generalized accelerations)
       ForwardDynamics (*model,q,qd,tau,qdd);
       
       //3g. Here qd, and qdd are used to populate dxdt
       //(the state derivative)
       j = 0;
       for(int i = 0; i < model->dof_count; i++){
           dxdt[j] = (double)qd[i];
           j++;
       }
       for(int i = 0; i < model->dof_count; i++){
           dxdt[j] = (double)qdd[i];
           j++;
       }

   }

public:
    Model* model;
    VectorNd q, qd, qdd, tau;
};

struct pushBackStateAndTime
{
    std::vector< state_type >& states;
    std::vector< double >& times;

   pushBackStateAndTime( std::vector< state_type > &states ,
                         std::vector< double > &times )
       : states( states ) , times( times ) { }

   void operator()( const state_type &x , double t )
   {
       states.push_back( x );
       times.push_back( t );
   }
};

void f(const state_type &x, state_type &dxdt, const double t);

rbdlToBoost* rbdlModel;
controlled_stepper_type* controlled_stepper;
double tp = 0;
int step_count =0;
state_type xState;
double kepe0 = 0, pe =0, ke =0;
// void f(const state_type &x, state_type &dxdt, const double t);
Matrix3d dynamics::rpyToMatrix(Vector3d rpy)
{
    Eigen::Matrix3d matrix;
    //!!! order same as in function eulerAngles
    matrix = Eigen::AngleAxisd(rpy[2], Eigen::Vector3d::UnitZ()) *
             Eigen::AngleAxisd(rpy[1], Eigen::Vector3d::UnitX()) *
             Eigen::AngleAxisd(rpy[0], Eigen::Vector3d::UnitZ());

             //!!! order 3,1,3 in eulerAngles difinition
    // std::cout<<matrix.transpose()<<std::endl;
    return matrix;
}

void dynamics::inverseDynamics(const VectorNd &Q,
                               const VectorNd &QDot,
                               const VectorNd &QDDot,
                               VectorNd &Tau,
                               std::vector<SpatialVector> *f_ext)
{
    UpdateKinematics(*model, q, qd, qdd);
    InverseDynamics(*model, Q, QDot, QDDot, Tau, f_ext);
}

void dynamics::nonlinearEffects (const Math::VectorNd &Q,
                                const Math::VectorNd &QDot,
                                Math::VectorNd &Tau,
                                std::vector<Math::SpatialVector> *f_ext)
{
    NonlinearEffects(*model, Q, QDot, Tau, f_ext);
}

void dynamics::forwardDynamics(
    const VectorNd &Q,
    const VectorNd &QDot,
    const VectorNd &Tau,
    VectorNd &QDDot,
    std::vector<SpatialVector> *f_ext)
{
    ForwardDynamics(*model, Q, QDot, Tau, QDDot);
}

bool dynamics::inverseKinematics(Vector3d pos, Matrix3d ori, const VectorNd qin, VectorNd &qres)
{
    // Vector3d local_point = Vector3d::Zero();
    //  std::cout << "inverseKinematics qin=" << q.transpose() << std::endl;
    InverseKinematicsConstraintSet cs;
    cs.AddFullConstraint(end, endpoint, pos, ori);
    cs.AddFullConstraint(end1, endpoint, pos, ori);

    bool result = InverseKinematics(*model, qin, cs, qres);
    // q = qres;
    //  std::cout << "inverseKinematics return " << result << " qout=" << qres.transpose() << std::endl;
    return true;
}

bool dynamics::inverseKinematics(Vector3d pos, Vector3d rpy, const VectorNd qin, VectorNd &qres)
{
    return inverseKinematics(pos, rpyToMatrix(rpy), qin, qres);
}

void dynamics::forwardKinematics(VectorNd qin, Vector3d &endpos, Matrix3d &endori)
{
    endpos = CalcBodyToBaseCoordinates(*model, qin, end, endpoint);
    endori = CalcBodyWorldOrientation(*model, qin, end);
}

void dynamics::forwardKinematics(VectorNd qin, Vector3d &endpos, Vector3d &endrpy)
{
    endpos = CalcBodyToBaseCoordinates(*model, qin, end, endpoint);
    endrpy = CalcBodyWorldOrientation(*model, qin, end).eulerAngles(0, 1, 2);
}

void dynamics::getEndPose(Vector3d &pos, Matrix3d &ori)
{
    pos = CalcBodyToBaseCoordinates(*model, q, end, endpoint);
    ori = CalcBodyWorldOrientation(*model, q, end);

    Vector3d euler = ori.eulerAngles(0, 1, 2);
}

void dynamics::getEndPose(Vector3d &pos, Vector3d &rpy)
{
    getEndPosition(pos);
    getEndOrientation(rpy);
}

void dynamics::getEndPosition(Vector3d &pos)
{
    pos = CalcBodyToBaseCoordinates(*model, q, end, endpoint);
}
void dynamics::getEndOrientation(Vector3d &rpy)
{
    rpy = CalcBodyWorldOrientation(*model, q, end).eulerAngles(0, 1, 2);
}

void dynamics::setEndPose(const Vector3d pos, const Matrix3d ori)
{
    VectorNd qres;
    inverseKinematics(pos, ori, q, qres);
}

void dynamics::setEndPose(const Vector3d pos, const Vector3d rpy)
{
    setEndPose(pos, rpyToMatrix(rpy));
}

void dynamics::setJointPositions(const VectorNd jp)
{
    q = jp;
}

void dynamics::getJointPositions(VectorNd &jp)
{
    jp = q;
}

int dynamics::getDofs()
{
    return model->dof_count;
}

void dynamics::resetDynamics()
{
    for(unsigned int i = 0; i< (model->dof_count); ++i){
        qd[i] = 0;
		qdd[i] = 0;
    }
	
	int j = 0;
    for(unsigned int i = 0; i< (model->dof_count); ++i){
        xState[j++] = q(i);
    }
    for(unsigned int i = 0; i< (model->dof_count); ++i){
        xState[j++] = qd(i);
    }
    tp = 0;
    step_count = 0;
    //std::cout << q.transpose() << std::endl;
	//std::cout << qd.transpose() << std::endl;
}



int dynamics::step(float dt)
{
    
    double t = tp + dt;

    //3h. Here we integrate forward in time between a series of nPts from
    //    t0 to t1
    // pe = Utils::CalcPotentialEnergy(*model, q, true);
    // ke = Utils::CalcKineticEnergy(*model, q, qd, true);
    // if(step_count==0) {kepe0 = (ke+pe);}
    // else{
    //     printf("%f, %f, %f, %f\n",
    //         t, ke,
    //         pe,(ke+pe-kepe0)); 
    // }
    if(step_count==0) {
        int j = 0;
        for(unsigned int i = 0; i< (model->dof_count); ++i){
            xState[j++] = q(i);
        }
        for(unsigned int i = 0; i< (model->dof_count); ++i){
            xState[j++] = qd(i);
        }
    }
    rbdlModel->tau = tau;
    integrate_adaptive( 
        *controlled_stepper ,
        *rbdlModel , xState , tp , t , (t-tp)/10 );
    tp = t;

    int j = 0;
    for(unsigned int k = 0; k < model->dof_count; ++k){

        q[k] = xState[j++];
    }
    for(unsigned int k = 0; k < model->dof_count; ++k){
        qd[k] = xState[j++];
    }
    qdd = rbdlModel->qdd;
    step_count++;
    //std::cout << qd.transpose() << std::endl;
    //UpdateKinematics(*model, q, qd, qdd);
    return 0;
}

void dynamics::printLinksCom()
{
    char lname[16];
    for(int i=0;i<model->dof_count;i++)
    {
        sprintf(lname,"link_%d",i+1);
        unsigned int linkid = model->GetBodyId(lname);
        Vector3d comlocal = model->mBodies[linkid].mCenterOfMass;
        Vector3d p1 = CalcBodyToBaseCoordinates(*model, q, linkid, comlocal , true);
        std::cout << lname << "com_local:" << comlocal.transpose() << "com_world:" << p1.transpose() << std::endl;
    }
}

bool dynamics::loadmodelfile(const char *modelfile)
{
    rbdl_check_api_version(RBDL_API_VERSION);

    model = new Model();

    // 3a. The URDF model is read in here, and turned into a series of
    //     vectors and matricies in model which RBDL uses to evaluate
    //     dynamics quantities

    if (!Addons::URDFReadFromFile(modelfile, model, false, false))
    {
        std::cerr << "Error loading model " << modelfile << std::endl;
        abort();
    }

    int dofs = model->dof_count;
    printf("ndof: %i\n", dofs);
    std::cout << "gravity:" << model->gravity.transpose() << std::endl;
    //model->gravity = Vector3d(0, 0, -9.81);
    //  std::cout << "Degree of freedom overview:" << std::endl;
    //  std::cout << Utils::GetModelDOFOverview(*model);
    std::cout << "Model Hierarchy:" << std::endl;
    std::cout << Utils::GetModelHierarchy(*model);
    q = VectorNd::Zero(dofs);
    qd = VectorNd::Zero(dofs);
    qdd = VectorNd::Zero(dofs);
    UpdateKinematics(*model, q, qd, qdd);


    end = model->GetBodyId("link_6");
    Vector3d p1 = CalcBodyToBaseCoordinates(*model, q, end, {0, 0, 0}, true);
    endpoint = {0,0,0.0}; // CalcBaseToBodyCoordinates(*model, q, end, tiptobase, false);
    

    //3b. Here we instantiate a wrapper class which is needed so that 
    //    Boost can evaluate the state derivative of the model.
    rbdlModel = new rbdlToBoost(model);
    

    //Integration settings
    double absTolVal   = 1e-10;
    double relTolVal   = 1e-6;
    double a_x = 1.0 , a_dxdt = 1.0;
    xState.resize(dofs*2);
    int j = 0;
    for(unsigned int i = 0; i< (model->dof_count); ++i){
        xState[j++] = q(i);
    }
    for(unsigned int i = 0; i< (model->dof_count); ++i){
        xState[j++] = qd(i);
    }


    controlled_stepper = new controlled_stepper_type  
    (
        default_error_checker< double , 
                               range_algebra , 
                               default_operations >
        ( absTolVal , relTolVal , a_x , a_dxdt ) );







    return true;
}

dynamics::dynamics()
{
    model = NULL;
}