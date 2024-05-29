#ifndef RIGIDBODY_H
#define RIGIDBODY_H

#include "geometry.h"

namespace backend {

class RigidBody {
private:
    std::vector<Vector2d> vertices;// (at theta = 0)
    std::vector<Face> faces;
    std::vector<Vector2d> r_relative;// from mass center to vertices (at theta = 0)
    FuncPtr get_sdf_func;// pointer of the function to compute sdf in S_c
    // generalized
    Vector3d x;
    Vector3d x0;
    Vector3d v;
    Matrix3d M;
    Matrix3d M_inv;

    Vector3d F();

    // for collision
    double e = 0.3f;
    double mu = 0.8f;

    void handleCollision(const MatrixXd& phi_solid, const double& l);    

public:
    RigidBody();
    RigidBody(const std::vector<Vector2d>& vertices_init,
              const std::vector<Face>& faces_init,
              const double& m_init, 
              const double& I_init, 
              const Vector2d& c_init, 
              const double& theta_init, 
              const Vector2d& velocity_init, 
              const double& omega_init,
              const FuncPtr& get_sdf_func_init);
    Matrix2d getRotationMatrix();
    Matrix2d getMass();
    double getI();
    Vector2d getC();
    Vector2d getVelocity();
    Vector2d getVelocity(const Vector2d& pos);
    double getOmega();
    Vector2d getVertexPosition(const int& idx);
    double getSDF(const Vector2d& pos);

    void setVelocity(const Vector2d& velocity);
    void setOmega(const double& omega);

    void setMass(const double& m_u, const double& m_v);

    void advance(const double& time_step, const MatrixXd& phi_solid, const double& l);


};

}


#endif