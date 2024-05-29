#ifndef FLUIDSIM_H
#define FLUIDSIM_H

#include "util.h"
#include "rigidbody.h"

namespace backend {

class FluidSim{
private:
    double l;
    int n1, n2;
    int n_liquid = 0;
    // velocity at the cell's center
    MatrixXd u, v;
    MatrixXd u_new, v_new;
    // level set
    MatrixXd phi;
    // boundary
    MatrixXd phi_solid;
    MatrixXi valid_u, valid_v;
    MatrixXi valid_u_old, valid_v_old;
    
    MatrixXd u_solid, v_solid;
    // for projection
    MatrixXd  weights_u, weights_v;
    MatrixXd  Adiag, Aplusi, Aplusj;
    MatrixXd  d, p;
    MatrixXd  z, s; // z: auxiliary vetor, s: search vector

    // particles
    double particle_radius;
    std::vector<Vector2d> particles;

    // rigid body
    bool add_rigidbody = false;
    RigidBody rigidbody;
    double rigidbody_density;
    double fluid_density;

    MatrixXd phi_rigidbody;
    MatrixXd phi_solid_rigidbody;

    MatrixXd J_x, J_y, J_rot;    

    MatrixXd weights_rigid_u, weights_rigid_v;

    // advance
    // 1. add force
    void applyForce(const double& dt);
    // 2. advect
    void advect(const double& dt);
    void advectParticles(const double& dt);
    Vector2d traceRk2(const Vector2d& position, const double& dt);
    void semiLagrangian(const MatrixXd& field, MatrixXd& field_new, const double& dt, const int& id);
    // 3. project
    void project();
    void computeWeights();        
    void solve(int maxIterations);
    void applyA(const MatrixXd& x, MatrixXd& ans);
    MatrixXd applyA(const MatrixXd& x);
    // some details
    void computePhi();
    void extrapolate();
    void extrapolate(MatrixXd& field, MatrixXd field_new, MatrixXi& valid, MatrixXi& valid_old);
    void constrain();
    void computeN();
    // rigid body
    void updateRigidBodyGrids();// update phi_rigidbody, phi_solid_rigidbody, weights_rigid
    void recomputeRigidBodyMass();
    void recomputeSolidVelocity();


public:
    // init
    FluidSim(const int& n1_init, const int& n2_init, const double& l_init, const MatrixXd& phi_init, const MatrixXd& phi_solid_init);
    void advance(const double& time_step);
    void setVelocity(const Vector2d& vec);

    Vector2d getVelocity(const Vector2d& pos);
    Vector2d getSolidVelocity(const Vector2d& pos);
    double getPhi(const Vector2d& pos);
    double getPhiSolid(const Vector2d& pos);
    double getPhiRigidBody(const Vector2d& pos);
    double getPhiSolidRigidBody(const Vector2d& pos);

    
    // rigid body
    void addRigdBody(
        const std::vector<Vector2d>& vertices_init,
        const std::vector<Face>& faces_init,
        const double& m_init, 
        const double& I_init, 
        const Vector2d& c_init, 
        const double& theta_init, 
        const Vector2d& velocity_init, 
        const double& omega_init,
        const FuncPtr& get_sdf_func_init
    );

    void run(double time_step, int n);// no UI
    void outputVideo(double time_step, int n, int time, int fps, std::string file_name = "output"); 
};

}

#endif