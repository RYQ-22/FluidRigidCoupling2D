#include "rigidbody.h"

namespace backend {

RigidBody::RigidBody() {}

RigidBody::RigidBody(const std::vector<Vector2d>& vertices_init,
              const std::vector<Face>& faces_init,
              const double& m_init, 
              const double& I_init, 
              const Vector2d& c_init, 
              const double& theta_init, 
              const Vector2d& velocity_init, 
              const double& omega_init,
              const FuncPtr& get_sdf_func_init) :
              vertices(vertices_init),
              faces(faces_init),
              get_sdf_func(get_sdf_func_init) {
        x = Vector3d(c_init(0), c_init(1), theta_init);
        x0 = x;
        v = Vector3d(velocity_init(0), velocity_init(1), omega_init);
        M.setZero(); M(0, 0) = m_init; M(1, 1) = m_init; M(2, 2) = I_init;
        M_inv = M.inverse();            
        for (int i = 0; i < vertices.size(); i++) {
            r_relative.push_back(vertices[i] - Vector2d(x(0), x(1)));
        }        
}

Matrix2d RigidBody::getRotationMatrix() {
    Matrix2d rotationMatrix;
    rotationMatrix(0, 0) = cos(x(2));
    rotationMatrix(1, 1) = cos(x(2));
    rotationMatrix(0, 1) = -sin(x(2));
    rotationMatrix(1, 0) = sin(x(2));
    return rotationMatrix;
}

Matrix2d RigidBody::getMass() {
    Matrix2d mass_matrix;
    mass_matrix.setZero();
    mass_matrix(0, 0) = M(0, 0);
    mass_matrix(1, 1) = M(1, 1);
    return mass_matrix;
}

double RigidBody::getI() {
    return M(2, 2);
}

Vector2d RigidBody::getC() {
    return Vector2d(x(0), x(1));
}

Vector2d RigidBody::getVelocity() {
    return Vector2d(v(0), v(1));
}

Vector2d RigidBody::getVelocity(const Vector2d& pos) {
    Vector2d r_rel = pos - getC();
    return getVelocity() + Vector2d(-r_rel(1), r_rel(0)) * v(2);
}

double RigidBody::getOmega() {
    return v(2);
}

Vector2d RigidBody::getVertexPosition(const int& idx) {
    return Vector2d(x(0), x(1)) + getRotationMatrix() * r_relative[idx];
}

Vector3d RigidBody::F() {
    return Vector3d(0, -0.5 * M(1, 1) * G, 0) - 0.001 * M * v;
}

void RigidBody::setVelocity(const Vector2d& velocity) {
    v(0) = velocity(0);
    v(1) = velocity(1);
    return;
}

void RigidBody::setOmega(const double& omega) {
    v(2) = omega;
    return;
}

void RigidBody::handleCollision(const MatrixXd& phi_solid, const double& l) {
    // handle collision
    Vector2d pos, grad;
    Vector2d vi, vi_n, vi_t;
    Matrix2d K;
    Vector2d r_rel, j, vi_new;
    for (int idx = 0; idx < vertices.size(); idx++) {        
        pos(0) = clamp(getVertexPosition(idx)(0), 0., (phi_solid.rows()-1-1e-6) * l);
        pos(1) = clamp(getVertexPosition(idx)(1), 0., (phi_solid.cols()-1-1e-6) * l);
        if (interpolate_value(pos / l, phi_solid) < 0) {
            interpolate_gradient(grad, pos / l, phi_solid);
            grad.normalize();
            vi = getVelocity(pos);
            if (vi.dot(grad) < 0) {// v_n < 0
                vi_n = vi.dot(grad) * grad;
                vi_t = vi - vi_n;
                K.setZero();
                r_rel = getRotationMatrix() * r_relative[idx];
                K(0, 0) += M_inv(0, 0); K(1, 1) += M_inv(1, 1);
                K(0, 0) += M_inv(2, 2) * r_rel(1) * r_rel(1);
                K(0, 1) += -M_inv(2, 2) * r_rel(0) * r_rel(1);
                K(1, 0) += -M_inv(2, 2) * r_rel(0) * r_rel(1);
                K(1, 1) += M_inv(2, 2) * r_rel(0) * r_rel(0);
                vi_new = -e * vi_n + max(0., 1. - mu * (1. + e) * vi_n.norm() / vi_t.norm()) * vi_t;
                j = K.inverse() * (vi_new - vi);

                v(2) += (r_rel(0) * j(1) - r_rel(1) * j(0)) * M_inv(2, 2);                
                v(0) += j(0) * M_inv(0, 0);
                v(1) += j(1) * M_inv(1, 1);
            }
        }
    }

    // project rigidbody out of solid
    double depth;
    for (int idx = 0; idx < vertices.size(); idx++) {
        pos(0) = clamp(getVertexPosition(idx)(0), 0., (phi_solid.rows()-1-1e-6) * l);
        pos(1) = clamp(getVertexPosition(idx)(1), 0., (phi_solid.cols()-1-1e-6) * l);
        depth = interpolate_value(pos / l, phi_solid) * l;
        if (depth < 0) {
            interpolate_gradient(grad, pos / l, phi_solid);
            grad.normalize();
            x(0) += -depth * grad(0);
            x(1) += -depth * grad(1);
        }
    }

    return;
}

void RigidBody::advance(const double& dt, const MatrixXd& phi_solid, const double& l) {
    // update x and v
    v += M_inv * F() * dt;
    x += v * dt;

    handleCollision(phi_solid, l);

    return;
}

double RigidBody::getSDF(const Vector2d& pos) {
    Vector2d pos_rel = getRotationMatrix().inverse() * (pos - getC());
    return get_sdf_func(pos_rel);
}

void RigidBody::setMass(const double& m_u, const double& m_v) {    
    M(0, 0) = m_u;
    M(1, 1) = m_v;
    return;
}

}
