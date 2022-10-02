#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <rbdl/rbdl.h>
#include <rbdl/rbdl_utils.h>
#ifndef RBDL_BUILD_ADDON_URDFREADER
#error "Error: RBDL addon URDFReader not enabled."
#endif
#include <rbdl/addons/urdfreader/urdfreader.h>
//using namespace boost::numeric::odeint;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;
// using namespace Eigen;
struct ctrl
{
    VectorNd Kp;
    VectorNd Ki;
    VectorNd Kd;
    VectorNd err;
    VectorNd err_last;
    VectorNd err_i;
    VectorNd err_d;
    VectorNd u;
};

class dynamics
{
public:
    dynamics();
    bool loadmodelfile(const char* modelfile);
    void inverseDynamics(    const VectorNd &Q,
                             const VectorNd &QDot,
                             const VectorNd &QDDot,
                             VectorNd &Tau,
                             std::vector<SpatialVector> *f_ext = NULL);

    void forwardDynamics (  const VectorNd &Q,
                            const VectorNd &QDot,
                            const VectorNd &Tau,
                            VectorNd &QDDot,
                            std::vector<SpatialVector> *f_ext = NULL);
    void nonlinearEffects ( const Math::VectorNd &Q,
                            const Math::VectorNd &QDot,
                            Math::VectorNd &Tau,
                            std::vector<Math::SpatialVector> *f_ext = NULL);
    bool inverseKinematics(Vector3d pos, Matrix3d ori,const VectorNd qin, VectorNd &qres);
    bool inverseKinematics(Vector3d pos, Vector3d rpy,const VectorNd qin, VectorNd &qres);
    void forwardKinematics(VectorNd qin, Vector3d& endpos,Matrix3d& endori);
    void forwardKinematics(VectorNd qin, Vector3d& endpos,Vector3d& endrpy);

    void getEndPose(Vector3d& pos, Matrix3d& ori);
    void getEndPose(Vector3d& pos, Vector3d& rpy);

    void setEndPose(const Vector3d pos, const Matrix3d ori);
    void setEndPose(const Vector3d pos, const Vector3d rpy);
    void getEndPosition(Vector3d& pos);
    void getEndOrientation(Vector3d& rpy);


    void setJointPositions(const VectorNd jp);
    void getJointPositions(VectorNd& jp);

    int getDofs();
    int step(float dt);
	void resetDynamics();
    void printLinksCom();
//added by yhzhu
//private:
    int armtype;
    void addParalelConstrant();
    Matrix3d rpyToMatrix(Vector3d rpy);

    Model* model;
    VectorNd    q;
    VectorNd    qd;
    VectorNd    qdd;
    VectorNd    tau;
    unsigned int    end;
    unsigned int    end1;
    Vector3d endpoint;

};

#endif // dynamics_H
