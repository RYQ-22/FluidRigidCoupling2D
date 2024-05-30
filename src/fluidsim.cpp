#include "fluidsim.h"

namespace backend {

FluidSim::FluidSim(
    const int& n1_init,
    const int& n2_init,    
    const double& l_init,
    const MatrixXd& phi_init,
    const MatrixXd& phi_solid_init) {
        Assert(n1_init+1 == phi_init.rows() && n2_init+1 == phi_init.cols(), "FluidSim::FluidSim", "The shape of phi_init is not correct.");
        Assert(n1_init+1 == phi_solid_init.rows() && n2_init+1 == phi_solid_init.cols(), "FluidSim::FluidSim", "The shape of phi_solid_init is not correct.");
        n1 = n1_init;
        n2 = n2_init;        
        l = l_init;  
        u.resize(n1+1, n2); u_new.resize(n1+1, n2); valid_u.resize(n1+1, n2); valid_u_old.resize(n1+1, n2); weights_u.resize(n1+1, n2); u_solid.resize(n1+1, n2);
        v.resize(n1, n2+1); v_new.resize(n1, n2+1); valid_v.resize(n1, n2+1); valid_v_old.resize(n1, n2+1); weights_v.resize(n1, n2+1); v_solid.resize(n1, n2+1);
        u_solid.setZero(); v_solid.setZero();
        Adiag.resize(n1, n2); Aplusi.resize(n1, n2); Aplusj.resize(n1, n2);
        d.resize(n1, n2); p.resize(n1, n2); z.resize(n1, n2); s.resize(n1, n2);        
        phi_solid = phi_solid_init;
        phi.resize(n1, n2);    
        n_liquid = 0;
        for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) {
            phi(i, j) = interpolate_value(i+0.5, j+0.5, phi_init);
            if (phi(i, j) < 0) {
                n_liquid++;
            }
        }        

    // make the particles large enough so they always appear on the grid
    particle_radius = (l * 1.01 * sqrt(2.0) / 2.0);
    // init particles
    int seed = 0;
	for (int n = 0; n < 2; n++) {
		for (int i = 1; i < n1-1; i++) for (int j = 1; j < n2-1; j++) {
			double a = randhashd(seed++, -1, 1); double b = randhashd(seed++, -1, 1);
			double x_ = static_cast<double>(i) + a, y_ = static_cast<double>(j) + b;
			if (interpolate_value(x_, y_, phi_init) <= -particle_radius) {		
                Assert(x_ > 0 && y_ > 0, "FluidSim::FluidSim", "x_ and y_ should be positive.");		
				if (getPhiSolid(Vector2d(x_ * l, y_ * l)) > 0) {
                    particles.push_back(Vector2d(x_ * l, y_ * l));
                }					
			}
		}
	}
    computePhi();    

    std::cout << "FluidSim is initialized." << std::endl;
}

void FluidSim::advance(const double& time_step){
    // debug
    Assert(time_step > 0, "FluidSim::advance", "Non-positive time step.");
    computeN();
    std::cout 
    << "n_liquid: " << n_liquid 
    << " V_liquid: " << max(u.cwiseAbs().maxCoeff() + v.cwiseAbs().maxCoeff(), 0.) 
    << " p_max: " << p.cwiseAbs().maxCoeff() 
    << " d_max: " << d.cwiseAbs().maxCoeff()
    << " rigidbody v_y: " << rigidbody.getVelocity()(1)
    << std::endl;

    double t = 0;
    bool done = false;
    while (!done) {
        // 5-dx rule
        double dt = 5 * l / (u.cwiseAbs().maxCoeff() + v.cwiseAbs().maxCoeff());
        if (t + dt >= time_step) {
            dt = time_step - t;
            done = true;
        }

        // core steps
        if (add_rigidbody) {
            rigidbody.advance(dt, phi_solid, l);
            updateRigidBodyGrids();
            recomputeRigidBodyMass();
        }
        advectParticles(dt);
        computePhi();
        computeWeights();
        advect(dt);
        applyForce(dt);
        extrapolate();// debug
        if (add_rigidbody) {
            recomputeSolidVelocity();
        }
        constrain();// debug
        project();
        extrapolate();        
        constrain();
        t += dt;
    } 
    return;
}

double FluidSim::getPhi(const Vector2d& pos) {  
    Vector2d position = pos;
    position(0) = clamp(position(0), l, (n1-1) * l);
    position(1) = clamp(position(1), l, (n2-1) * l); 
    return interpolate_value(position / l - Vector2d(0.5, 0.5), phi);
}

double FluidSim::getPhiSolid(const Vector2d& pos) {    
    Assert(pos(0) >= 0 && pos(1) >= 0, "FluidSim::getPhiSolid", "pos(0) and pos(1) should >= 0"); 
    Assert(pos(0) / l < phi_solid.rows()-1 && pos(1) / l < (double)phi_solid.cols()-1, "FluidSim::getPhiSolid", "out of index");  
    return interpolate_value(pos / l, phi_solid);
}

double FluidSim::getPhiRigidBody(const Vector2d& pos) {
    Assert(add_rigidbody, "FluidSim::getPhiRigidBody", "rigidbody not defined.");
    Assert(pos(0) >= 0 && pos(1) >= 0, "FluidSim::getPhiRigidBody", "pos(0) and pos(1) should >= 0"); 
    Assert(pos(0) / l < phi_rigidbody.rows()-1 && pos(1) / l < phi_rigidbody.cols()-1, "FluidSim::getPhiRigidBody", "out of index");    
    return interpolate_value(pos / l, phi_rigidbody);
}

double FluidSim::getPhiSolidRigidBody(const Vector2d& pos) {
    Assert(add_rigidbody, "FluidSim::getPhiSolidRigidBody", "rigidbody not defined.");
    Assert(pos(0) >= 0 && pos(1) >= 0, "FluidSim::getPhiSolidRigidBody", "pos(0) and pos(1) should >= 0"); 
    Assert(pos(0) / l < phi_solid_rigidbody.rows()-1 && pos(1) / l < phi_solid_rigidbody.cols()-1, "FluidSim::getPhiSolidRigidBody", "out of index");  
    return interpolate_value(pos / l, phi_rigidbody);
}

Vector2d FluidSim::getVelocity(const Vector2d& pos) {
    Vector2d lattice_position = pos / l;
    // project_to_interior
    lattice_position(0) = clamp(lattice_position(0), 1., n1-1.);
    lattice_position(1) = clamp(lattice_position(1), 1., n2-1.);    
    
    double u_ave = interpolate_value(lattice_position(0), lattice_position(1) - 0.5, u);
    double v_ave = interpolate_value(lattice_position(0) - 0.5, lattice_position(1), v);
    
    return Vector2d(u_ave, v_ave);
}

Vector2d FluidSim::getSolidVelocity(const Vector2d& pos) {
    Vector2d lattice_position = pos / l;
    // project_to_interior
    lattice_position(0) = clamp(lattice_position(0), 1., n1-1.);
    lattice_position(1) = clamp(lattice_position(1), 1., n2-1.);    
    
    double u_ave = interpolate_value(lattice_position(0), lattice_position(1) - 0.5, u_solid);
    double v_ave = interpolate_value(lattice_position(0) - 0.5, lattice_position(1), v_solid);
    
    return Vector2d(u_ave, v_ave);
}

Vector2d FluidSim::traceRk2(const Vector2d& position, const double& dt){
    Vector2d input = position;
	Vector2d velocity = getVelocity(input);
	velocity = getVelocity(input + 0.5 * dt * velocity);
	input += dt * velocity;
	return input;
}

MatrixXd FluidSim::applyA(const MatrixXd& x){
    MatrixXd ans(x.rows(), x.cols());
    ans.setZero();
    for (int i = 0; i < x.rows(); i++) for (int j = 0; j < x.cols(); j++) {
		ans(i, j) = Adiag(i, j) * x(i, j);
        if(i < x.rows()-1)
		    ans(i, j) += Aplusi(i, j) * x(i + 1, j);
        if(j < x.cols()-1)    
		    ans(i, j) += Aplusj(i, j) * x(i, j + 1);        
        if(i > 0)
		    ans(i, j) += Aplusi(i - 1, j) * x(i - 1, j);
        if(j > 0)
		    ans(i, j) += Aplusj(i, j - 1) * x(i, j - 1);       
	}
    if (add_rigidbody) {
        ans += fluid_density * (J_x.reshaped()).dot(x.reshaped()) * J_x / rigidbody.getMass()(0, 0);
        ans += fluid_density * (J_y.reshaped()).dot(x.reshaped()) * J_y / rigidbody.getMass()(1, 1);
        ans += fluid_density * (J_rot.reshaped()).dot(x.reshaped()) * J_rot / rigidbody.getI();
    }
    return ans;
}

void FluidSim::applyA(const MatrixXd& x, MatrixXd& ans){
    ans.resize(x.rows(), x.cols());
    ans.setZero();
    for (int i = 0; i < x.rows(); i++) for (int j = 0; j < x.cols(); j++) {
		ans(i, j) = Adiag(i, j) * x(i, j);
        if(i < x.rows()-1)
		    ans(i, j) += Aplusi(i, j) * x(i + 1, j);
        if(j < x.cols()-1)    
		    ans(i, j) += Aplusj(i, j) * x(i, j + 1);        
        if(i > 0)
		    ans(i, j) += Aplusi(i - 1, j) * x(i - 1, j);
        if(j > 0)
		    ans(i, j) += Aplusj(i, j - 1) * x(i, j - 1);       
	}
    if (add_rigidbody) {
        ans += fluid_density * (J_x.reshaped()).dot(x.reshaped()) * J_x / rigidbody.getMass()(0, 0);
        ans += fluid_density * (J_y.reshaped()).dot(x.reshaped()) * J_y / rigidbody.getMass()(1, 1);
        ans += fluid_density * (J_rot.reshaped()).dot(x.reshaped()) * J_rot / rigidbody.getI();            
    }
    return;
}


void FluidSim::extrapolate(){
    valid_u.setZero();
    valid_v.setZero();    
    // reset valid_u
    for (int i = 1; i < n1; i++) for (int j = 0; j < n2; j++) {
        if (weights_u(i, j) > 0 && (phi(i, j) < 0.0 || phi(i - 1, j) < 0.0)) {
			valid_u(i, j) = 1;
		}		
    }
    // reset valid_v
    for (int i = 0; i < n1; i++) for (int j = 1; j < n2; j++) {
        if (weights_v(i, j) > 0 && (phi(i, j) < 0 ||  phi(i, j - 1) < 0)) {
			valid_v(i, j) = 1;
		}		
    }   
    // set u = 0 if valid_u = 0
    for (int idx = 0; idx < u.size(); idx++){
        if(valid_u(idx) == 0){
            u(idx) = 0;
        }
    }
    // set v = 0 if valid_v = 0
    for (int idx = 0; idx < v.size(); idx++){
        if(valid_v(idx) == 0){
            v(idx) = 0;
        }
    }    

	// Apply several iterations in all directions
	extrapolate(u, u_new, valid_u, valid_u_old);
	extrapolate(v, v_new, valid_v, valid_v_old);
    return;
}

void FluidSim::extrapolate(MatrixXd& field, MatrixXd field_new, MatrixXi& valid, MatrixXi& valid_old){
    field_new = field;
	for (int num = 0; num < 20; num++) {
		double sum = 0;
		int count = 0;
		valid_old = valid;
		for (int i = 1; i < field.rows() - 1; i++) for (int j = 1; j < field.cols() - 1; j++) {
			sum = 0;
			count = 0;
			if (!valid_old(i, j)) {

				if (valid_old(i + 1, j)) {
					sum += field(i + 1, j);
					count++;
				}
				if (valid_old(i - 1, j)) {
					sum += field(i - 1, j);
					count++;
				}
				if (valid_old(i, j + 1)) {
					sum += field(i, j + 1);
					count++;
				}
				if (valid_old(i, j - 1)) {
					sum += field(i, j - 1);
					count++;
				}				

				if (count > 0) {
					field_new(i, j) = sum / count;
					valid(i, j) = 1;
				}
			}
		}
		field = field_new;
	}
    return;
}

void FluidSim::advectParticles(const double& dt){
    for (int i=0; i<particles.size(); i++) {
        particles[i] = traceRk2(particles[i], dt);
        particles[i](0) = clamp(particles[i](0), l, (n1-1) * l);
        particles[i](1) = clamp(particles[i](1), l, (n2-1) * l);
        // solid boundary        
        double phi_val = getPhiSolid(particles[i]);
        if (add_rigidbody) {
            phi_val = getPhiSolidRigidBody(particles[i]);  
        }
        if(phi_val < 0){
            Vector2d grad;
            interpolate_gradient(grad, particles[i] / l, phi_solid);
            if (add_rigidbody) {
                interpolate_gradient(grad, particles[i] / l, phi_solid_rigidbody);
            }
            grad.normalize();
            particles[i] -= phi_val * grad;
        }
    }
    return;
}

void FluidSim::computePhi(){
    phi.setConstant(3*l);
    for(int p=0; p<particles.size(); p++){
        Vector2i cell_ind = (particles[p] / l - Vector2d(0.5, 0.5)).cast<int>();        
        for (int j = max(0, cell_ind[1] - 1); j <= min(cell_ind[1] + 1, n2 - 1); j++) 
        for (int i = max(0, cell_ind[0] - 1); i <= min(cell_ind[0] + 1, n1 - 1); i++) {
            Vector2d sample_pos((i + 0.5) * l, (j + 0.5) * l);                    
            double test_val = (sample_pos - particles[p]).norm() - particle_radius;
            if (test_val < phi(i, j)) {
                phi(i, j) = test_val;
            }
        }            		
    }
    
    // extend phi into solid
    for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) {
		if (phi(i, j) < 0.5 * l && getPhiSolid(Vector2d(i+0.5, j+0.5) * l) < 0) {
			phi(i, j) = -0.5 * l;
		}
	}

    return;
}

void FluidSim::advect(const double& dt){
    semiLagrangian(u, u_new, dt, 0);
    semiLagrangian(v, v_new, dt, 1);    
    u = u_new;
    v = v_new;    
    return;
}

void FluidSim::semiLagrangian(const MatrixXd& field, MatrixXd& field_new, const double& dt, const int& id) {
    assert(id == 0 || id == 1);
    Vector2d offset;
    if(id == 0){
        offset = Vector2d(0, 0.5);
    }
    else if(id == 1){
        offset = Vector2d(0.5, 0);
    }   
    for (int i = 1; i < field.rows()-1; i++) for (int j = 1; j < field.cols()-1; j++) {
		Vector2d pos = Vector2d(i * l, j * l) + offset * l;        
        field_new(i, j) = getVelocity(traceRk2(pos, -dt))(id);
	}
    return;
}

void FluidSim::applyForce(const double& dt){
    v.array() += -G * dt;
    return;
}

void FluidSim::computeWeights(){
    // compute weights_u
    for(int i=0; i<n1+1; i++) for(int j=0; j<n2; j++) {
        weights_u(i, j) = 1 - computeFraction(phi_solid(i, j), phi_solid(i, j+1));
        if (add_rigidbody) {
            weights_u(i, j) -= weights_rigid_u(i, j);
        }
        weights_u(i, j) = clamp(weights_u(i, j), 0., 1.);
    }
    // compute weights_v
    for(int i=0; i<n1; i++) for(int j=0; j<n2+1; j++) {
        weights_v(i, j) = 1 - computeFraction(phi_solid(i, j), phi_solid(i+1, j));
        if (add_rigidbody) {
            weights_v(i, j) -= weights_rigid_v(i, j);
        }
        weights_v(i, j) = clamp(weights_v(i, j), 0., 1.);
    }   
    return;
}


void FluidSim::project(){
    // compute J_x, J_y, J_rot
    if (add_rigidbody) {
        J_x.setZero(); J_y.setZero(); J_rot.setZero();        
        for(int i = 0; i < n1; i++) for(int j = 0; j < n2; j++) {
            if (phi(i, j) < 0) {
                J_x(i, j) = weights_rigid_u(i+1, j) - weights_rigid_u(i, j);
                J_y(i, j) = weights_rigid_v(i, j+1) - weights_rigid_v(i, j);
                Vector2d pos((i+0.5) * l, (j+0.5) * l);
                pos = pos - rigidbody.getC();
                J_rot(i, j) = pos(0) * J_y(i, j) - pos(1) * J_x(i, j);
            }            
        }         
    }

    // compute A and d
    Adiag.setZero(); Aplusi.setZero(); Aplusj.setZero(); d.setZero();
    double theta = 0;
    for(int i = 0; i < n1; i++) for(int j = 0; j < n2; j++) {
        if(phi(i, j) < 0){
            // right
            if(weights_u(i+1, j) > 0){
                if(phi(i+1, j) < 0){
                    Adiag(i, j) += weights_u(i+1, j);
                    Aplusi(i, j) += -weights_u(i+1, j);
                }
                else{
                    theta = computeFraction(phi(i, j), phi(i+1, j));
                    if(theta < 0.01) theta = 0.01;                    
                    Adiag(i, j) += weights_u(i+1, j) / theta;
                }
                d(i, j) += -weights_u(i+1, j) * u(i+1, j);
            }
            // left            
            if(weights_u(i, j) > 0){
                if(phi(i-1, j) < 0){
                    Adiag(i, j) += weights_u(i, j);
                }
                else{
                    theta = computeFraction(phi(i, j), phi(i-1, j));
                    if(theta < 0.01) theta = 0.01;                
                    Adiag(i, j) += weights_u(i, j) / theta;
                }
                d(i, j) += weights_u(i, j) * u(i, j);
            }
            // top
            if(weights_v(i, j+1) > 0){
                if(phi(i, j+1) < 0){
                    Adiag(i, j) += weights_v(i, j+1);
                    Aplusj(i, j) += -weights_v(i, j+1);
                }
                else{
                    theta = computeFraction(phi(i, j), phi(i, j+1));
                    if(theta < 0.01) theta = 0.01;                    
                    Adiag(i, j) += weights_v(i, j+1) / theta;
                }
                d(i, j) += -weights_v(i, j+1) * v(i, j+1);
            }
            // bottom            
            if(weights_v(i, j) > 0){
                if(phi(i, j-1) < 0){
                    Adiag(i, j) += weights_v(i, j);
                }
                else{
                    theta = computeFraction(phi(i, j), phi(i, j-1));
                    if(theta < 0.01) theta = 0.01;                    
                    Adiag(i, j) += weights_v(i, j) / theta;
                }
                d(i, j) += weights_v(i, j) * v(i, j);
            }
            
        }        
    }

    // add rigidbody terms
    if (add_rigidbody) {
        // rhs
        d += -J_x * rigidbody.getVelocity()(0);
        d += -J_y * rigidbody.getVelocity()(1);
        d += -J_rot * rigidbody.getOmega();
    }
    
    // ensure A is positive definite
    for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) {
        if(Adiag(i, j)==0 && Aplusi(i, j)==0 && Aplusj(i, j)==0 && (i==0 || Aplusi(i-1, j)==0) && (j==0 || Aplusj(i, j-1)==0)){
            Adiag(i, j) = 1;
            d(i, j) = 0;
        }		
	}
    // TODO: remove the 1D null space
    // solve A * p = d
    solve(1000);

    // update u
    for (int i = 1; i < n1; i++) for (int j = 0; j < n2; j++) {
        if (weights_u(i, j) > 0 && (phi(i, j) < 0 || phi(i - 1, j) < 0)) {
			theta = 1;
            if (phi(i, j) > 0 || phi(i - 1, j) > 0) theta = computeFraction(phi(i, j), phi(i - 1, j));
            if (theta < 0.01) theta = 0.01;
            //theta = 1;// debug
            u(i, j) += -((p(i, j) - p(i - 1, j))) / theta;
		}		
    }
    // update v
    for (int i = 0; i < n1; i++) for (int j = 1; j < n2; j++) {
        if (weights_v(i, j) > 0 && (phi(i, j) < 0 ||  phi(i, j - 1) < 0)) {
			theta = 1;
            if (phi(i, j) > 0 || phi(i, j - 1) > 0) theta = computeFraction(phi(i, j), phi(i, j - 1));
            if (theta < 0.01) theta = 0.01;
            //theta = 1;// debug
            v(i, j) += -((p(i, j) - p(i, j - 1))) / theta;
		}		
    }
    // update rigidbody
    if (add_rigidbody) {        
        Vector2d velocity_temp = rigidbody.getVelocity();
        double omega_temp = rigidbody.getOmega();
        for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) {
            if (phi(i, j) < 0) {
                velocity_temp(0) += fluid_density * J_x(i, j) * p(i, j) / rigidbody.getMass()(0, 0);
                velocity_temp(1) += fluid_density * J_y(i, j) * p(i, j) / rigidbody.getMass()(1, 1);
                omega_temp += fluid_density * J_rot(i, j) * p(i, j) / rigidbody.getI();
            }
        }
        rigidbody.setVelocity(velocity_temp);
        rigidbody.setOmega(omega_temp);        
    }

    return;
}

void FluidSim::solve(int maxIterations){
    // set p = 0
	//p.setConstant(0);
	if (d.cwiseAbs().maxCoeff() == 0) {		
		p.setConstant(0);
        return;
	}
    d = d - applyA(p);
	// TODO
	// applyPrecon();
    // ######
    z = d;
	s = z;
	double sigma = (z.reshaped()).dot(d.reshaped());
	double sigma_new;
	double alpha;
    while (d.cwiseAbs().maxCoeff() > 1e-6) {
        Assert(maxIterations >= 0, "FluidSim::Solve", "CG failed.");
        applyA(s, z);
		alpha = sigma / (z.reshaped()).dot(s.reshaped());
		p += s * alpha;
		d -= z * alpha;		
        // TODO
		//applyPrecon();
        // ######		
		z = d;
		sigma_new = (z.reshaped()).dot(d.reshaped());
		alpha = sigma_new / sigma;
		sigma = sigma_new;        		
        s = z + alpha * s;
        maxIterations--;
    }	
	return;
}

void FluidSim::constrain(){
    u_new = u;
    v_new = v;    

    Vector2d pos, vec, vec_solid, normal_solid;
    // constrain u
    for(int i=1; i<n1; i++) for(int j=1; j<n2-1; j++) {
        if(weights_u(i, j) == 0){
            pos = Vector2d(i, j+0.5) * l;
            vec = getVelocity(pos);
            vec_solid = getSolidVelocity(pos);
            if (add_rigidbody) {
                interpolate_gradient(normal_solid, pos / l, phi_solid_rigidbody);
            }
            else {
                interpolate_gradient(normal_solid, pos / l, phi_solid);
            }
            normal_solid.normalize();            
            vec -= vec.dot(normal_solid) * normal_solid;
            vec += vec_solid.dot(normal_solid) * normal_solid;
            u_new(i, j) = vec(0);
        }
    }
    // constrain v
    for(int i=1; i<n1-1; i++) for(int j=1; j<n2; j++) {
        if(weights_v(i, j) == 0){
            pos = Vector2d(i+0.5, j) * l;
            vec = getVelocity(pos);
            vec_solid = getSolidVelocity(pos);
            if (add_rigidbody) {
                interpolate_gradient(normal_solid, pos / l, phi_solid_rigidbody);
            }
            else {
                interpolate_gradient(normal_solid, pos / l, phi_solid);
            }
            normal_solid.normalize();
            vec -= vec.dot(normal_solid) * normal_solid;
            vec += vec_solid.dot(normal_solid) * normal_solid;
            v_new(i, j) = vec(1);
        }
    }    
    // boundary
    for(int i=0; i<n1+1; i++) for(int j=0; j<n2; j++) if(i==0 || i==n1) 
        u_new(i, j) = 0;
    for(int i=0; i<n1; i++) for(int j=0; j<n2+1; j++) if(j==0 || j==n2) 
        v_new(i, j) = 0;

    u = u_new;
    v = v_new;    
    return;
}

void FluidSim::computeN(){
    n_liquid = 0;
    for(int i=0; i<n1; i++) for(int j=0; j<n2; j++) {
        if(phi(i, j) <= 0) n_liquid++; 
    }
}

void FluidSim::setVelocity(const Vector2d& vec){
    u.setConstant(vec(0));
    v.setConstant(vec(1));
    return;
}

// rigid body

void FluidSim::addRigdBody(
        const std::vector<Vector2d>& vertices_init,
        const std::vector<Face>& faces_init,
        const double& m_init, 
        const double& I_init, 
        const Vector2d& c_init, 
        const double& theta_init, 
        const Vector2d& velocity_init, 
        const double& omega_init,
        const FuncPtr& get_sdf_func_init
    ) {
        Assert(!add_rigidbody, "FluidSim::addRigidBody", "Only one rigid body is supported.");
        std::cout << "Rigidbody is added." << std::endl;
        add_rigidbody = true;

        // init 
        rigidbody = RigidBody(vertices_init, faces_init, m_init, I_init, c_init, theta_init, velocity_init, omega_init, get_sdf_func_init);
        J_x.resize(n1, n2); J_y.resize(n1, n2); J_rot.resize(n1, n2);
        phi_rigidbody.resize(n1+1, n2+1); phi_solid_rigidbody.resize(n1+1, n2+1);
        weights_rigid_u.resize(n1+1, n2); weights_rigid_v.resize(n1, n2+1);

        updateRigidBodyGrids();
        // compute rigidbody_density
        rigidbody_density = (rigidbody.getMass()(0, 0) + rigidbody.getMass()(1, 1)) / (weights_rigid_u.array().abs().sum() + weights_rigid_v.array().abs().sum());
        fluid_density = rigidbody_density * 0.7;

        return;
}

void FluidSim::updateRigidBodyGrids() {
    // update phi_rigidbody
    for (int i = 0; i < n1+1; i++) for (int j = 0; j < n2+1; j++) {
        Vector2d pos(i * l, j * l);
        phi_rigidbody(i, j) = rigidbody.getSDF(pos);
    }
    // update phi_solid_rigidbody (no intersection)
    for (int i = 0; i < n1+1; i++) for (int j = 0; j < n2+1; j++) {
        if (abs(phi_rigidbody(i, j)) > abs(phi_solid(i, j))) {
            phi_solid_rigidbody(i, j) = phi_solid(i, j);
        }
        else {
            phi_solid_rigidbody(i, j) = phi_rigidbody(i, j);
        }
    }    
    // compute weights_rigid_u
    for(int i=0; i<n1+1; i++) for(int j=0; j<n2; j++) {
        weights_rigid_u(i, j) = computeFraction(phi_rigidbody(i, j), phi_rigidbody(i, j+1));
        weights_rigid_u(i, j) = clamp(weights_rigid_u(i, j), 0., 1.);
    }
    // compute weights_rigid_v
    for(int i=0; i<n1; i++) for(int j=0; j<n2+1; j++) {
        weights_rigid_v(i, j) = computeFraction(phi_rigidbody(i, j), phi_rigidbody(i+1, j));
        weights_rigid_v(i, j) = clamp(weights_rigid_v(i, j), 0., 1.);
    }    
    return;
}

void FluidSim::recomputeRigidBodyMass() {
    // recompute effective mass
    rigidbody.setMass(rigidbody_density * weights_rigid_u.array().abs().sum(), rigidbody_density * weights_rigid_v.array().abs().sum());
}

void FluidSim::recomputeSolidVelocity() {
    // compute u_solid
    u_solid.setZero();
    for(int i=0; i < n1+1; i++) for(int j = 0; j < n2; j++) {
        Vector2d pos(i * l, (j+0.5) * l);        
        pos(0) = clamp(pos(0), 0., (n1-1e-6) * l);
        pos(1) = clamp(pos(1), 0., (n2-1-1e-6) * l);                    
        if (getPhiSolid(pos) < 0) {
            u_solid(i, j) = 0.;
        }
        else if (getPhiRigidBody(pos) < 0) {
            u_solid(i, j) = rigidbody.getVelocity(pos)(0);
        }        
    }
    
    // compute v_solid
    v_solid.setZero();
    for(int i=0; i < n1; i++) for(int j = 0; j < n2+1; j++) {
        Vector2d pos((i+0.5) * l, j * l);
        pos(0) = clamp(pos(0), 0., (n1-1e-6) * l);
        pos(1) = clamp(pos(1), 0., (n2-1-1e-6) * l); 
        if (getPhiSolid(pos) < 0) {
            v_solid(i, j) = 0.;            
        }
        else if (getPhiRigidBody(pos) < 0) {
            v_solid(i, j) = rigidbody.getVelocity(pos)(1);
        }
    }
    return;    
}

void FluidSim::run(double time_step, int n) {
    while (1) {				
		for (int i = 0; i < n; i++) {
			advance(time_step);
		}
	}	
    return;
}

void FluidSim::outputVideo(double time_step, int n, int time, int fps, std::string file_name) {
    int resolution = 40;
    cv::Mat image = cv::Mat::zeros(resolution* n1, resolution * n2, CV_8UC3);
    cv::VideoWriter writer;
    writer.open(file_name + ".avi", cv::VideoWriter::fourcc('M', 'J', 'P', 'G'), fps, cv::Size(resolution * n1, resolution * n2));
    Vector2d pos;
    int frame_number = time * fps;
    for (int idx = 0; idx < frame_number; idx++) {
        
        for (int i = 0; i < resolution * n1; i++) for (int j = 0; j < resolution * n2; j++) {
            pos = Vector2d(i+0.5, j+0.5) / resolution * l;
            if (getPhiRigidBody(pos) <= 0) {
                image.at<cv::Vec3b>(resolution * n2 - j, i) = cv::Vec3b(200, 200, 200);// BGR
            }
            else if (getPhi(pos) <= 0 && getPhiSolid(pos) > 0) {
                image.at<cv::Vec3b>(resolution * n2 - j, i) = cv::Vec3b(200, 0, 0);// BGR
            }
            else {
                image.at<cv::Vec3b>(resolution * n2 - j, i) = cv::Vec3b(0, 0, 0);
            }
        }
        

        /*
        for (int i = 0; i < 40 * n1; i++) for (int j = 0; j < 40 * n2; j++) {            
            image.at<cv::Vec3b>(40 * n2 - j, 40 * n1 -i) = cv::Vec3b(0, 0, 0);            
        }
        for (int idx = 0; idx < particles.size(); idx++) {

        }
        */
        std::cout << "frame: " << idx << std::endl;
        writer.write(image);
        for (int k = 0; k < n; k++) {
            advance(time_step);
        }
    }
    //writer.release();
    return;
}


}