#include "fluidsim.h"
#include "geometry.h"

int n1 = 120;
int n2 = 120;
double l = 1. / 120.;

double rigidboydSDF(const backend::Vector2d& pos) {
    return backend::computeSDF_cuboid(pos, backend::Vector2d(static_cast<double>(n1) / 6, static_cast<double>(n2) / 6) * l);
}

int main(){ 
    backend::MatrixXd phi_solid = -backend::computeSDF_cuboid(n1, n2, backend::Vector2i(4, 4), backend::Vector2i(n1-8, n2-8));
    backend::MatrixXd phi_fluid = backend::computeSDF_circle(n1, n2, backend::Vector2d((double)n1 / 3, (double)n2 / 3), 24);
    backend::FluidSim sim(n1, n2, l, phi_fluid, phi_solid);
    
    // add rigidbody
    std::vector<backend::Vector2d> vertices;
    std::vector<backend::Face> faces;
    vertices.push_back(backend::Vector2d(n1 * l * 7. / 12., n2 * l * 7. / 12.));// 0
    vertices.push_back(backend::Vector2d(n1 * l * 9. / 12., n2 * l * 7. / 12.));// 1
    vertices.push_back(backend::Vector2d(n1 * l * 9. / 12., n2 * l * 9. / 12.));// 2
    vertices.push_back(backend::Vector2d(n1 * l * 7. / 12., n2 * l * 9. / 12.));// 3
    faces.push_back(backend::Face(0, 1));
    faces.push_back(backend::Face(1, 2));
    faces.push_back(backend::Face(2, 3));
    faces.push_back(backend::Face(3, 0));
    
    
    sim.addRigdBody(
        vertices,
        faces,
        1.,
        2. / 216,
        backend::Vector2d(n1 * l * 2. / 3., n2 * l * 2. / 3.),
        0.2,
        backend::Vector2d(0, 0),
        0.,
        rigidboydSDF
    );
    

    //sim.run(0.002, 1);
    sim.outputVideo(0.005, 10, 5, 20, "coupling7");

    return 0;
}