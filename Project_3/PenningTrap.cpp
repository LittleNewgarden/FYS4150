#include <iostream>
#include <armadillo>

#include "PenningTrap.hpp"



// function defining the Penning trap system
PenningTrap::PenningTrap(double B0, double V0, double d, double f, double Wv)
{
    B0_= B0;
    V0_ = V0;
    d_ = d;
    f_ = f;
    Wv_ = Wv;

}


// Add a particle to the trap
void PenningTrap::add_particle(Particle p)
{
    All_Particles.push_back(p);
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r, bool time_dependance, double t)
{
    // different potential if there is time dependance
    if (time_dependance == true)
    {
        double V = V0_*(1 + f_*std::cos(Wv_*t));
        arma::vec E = V/(d_*d_)*arma::vec{r[0],r[1],-2*r[2]};
        return E;
    }   
    else
    
    {
        arma::vec E = V0_/(d_*d_)*arma::vec{r[0],r[1],-2*r[2]};
        return E;
    }
}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j)
{
    // extracting properties of particle i and j
    arma::vec r = All_Particles[i].r_; 
    arma::vec rj = All_Particles[j].r_; 
    int q = All_Particles[i].q_; 
    int qj = All_Particles[j].q_;

    arma::vec r_rj = r-rj;
    double norm = arma::norm(r_rj);
    
    arma::vec E = Ke_*qj*(r_rj)/(norm*norm*norm);

    arma::vec F = q*E;

    return F;
}

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i, bool time_dependance ,  double t)
{
    // extracting properties of particle i
    arma::vec r = All_Particles[i].r_; 
    arma::vec v = All_Particles[i].v_; 

    double norm = arma::norm(r);

    // check if particle i is inside the Penning trap
    if (norm>d_)
    {
        arma::vec F =  arma::vec{0,0,0};
        return F;
    }

    else
    {
        int q =  All_Particles[i].q_; 
        arma::vec F = q*external_E_field(r, time_dependance, t) + q*arma::vec{v[1]*B0_, -v[0]*B0_, 0};
        return F;
    }
    

}

// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i)
{
    // summing over the force from all particles
    arma::vec F = arma::vec{0, 0, 0};
    for (int j=0; j<  (int)  All_Particles.size(); j++)
    {
        if (j!=i)
        {
            F += force_particle(i, j);
        }
    }
    return F;
}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i, bool interactions, bool time_dependance, double t)
{
    // different force depending on wether Coloumb forces between particles are considerd
    if (interactions==true)
    {
        arma::vec F = total_force_external(i, time_dependance, t) + total_force_particles(i);
        return F;
    }

    else 
    {
        arma::vec F = total_force_external(i, time_dependance, t);
        return F;
    }

}

 // Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt, bool interactions, bool time_dependance , double t)
{
        arma::mat r = arma::mat(3,All_Particles.size());
        arma::mat v = arma::mat(3,All_Particles.size());

        // save the possitions and velocities of all particles before we evlove the system in time
        for (int i = 0; i<  (int)  All_Particles.size(); i++)
        {
            r.col(i) = All_Particles[i].r_;
            v.col(i) = All_Particles[i].v_;   
        }

        arma::mat k1_r = arma::mat(3,All_Particles.size());     
        arma::mat k2_r = arma::mat(3,All_Particles.size());
        arma::mat k3_r = arma::mat(3,All_Particles.size());
        arma::mat k4_r = arma::mat(3,All_Particles.size());

        arma::mat k1_v = arma::mat(3,All_Particles.size());     
        arma::mat k2_v = arma::mat(3,All_Particles.size());
        arma::mat k3_v = arma::mat(3,All_Particles.size());
        arma::mat k4_v = arma::mat(3,All_Particles.size());

        // finding k1 for the two coupled diff equations
        for (int i = 0; i<  (int)  All_Particles.size(); i++)
        {
            k1_r.col(i) = All_Particles[i].v_*dt;
            k1_v.col(i) = total_force(i, interactions, time_dependance, t)/All_Particles[i].m_*dt;
        }

        // uppdating the possition and velocity using k1
        for (int i = 0; i<  (int)  All_Particles.size(); i++)
        {
            All_Particles[i].r_ = r.col(i) + 0.5*k1_r.col(i);
            All_Particles[i].v_ = v.col(i) + 0.5*k1_v.col(i);
        }

        // finding k2 using updated possitions        
        for (int i = 0; i<  (int)  All_Particles.size(); i++)
        {
            k2_r.col(i) = All_Particles[i].v_*dt;
            k2_v.col(i) = total_force(i, interactions, time_dependance, t + dt/2)/All_Particles[i].m_*dt;
        }

        // uppdating the possition and velocity starting form the initial r and v vectors using k2
        for (int i = 0; i<  (int)  All_Particles.size(); i++)
        {
            All_Particles[i].r_ = r.col(i) + 0.5*k2_r.col(i);
            All_Particles[i].v_ = v.col(i) + 0.5*k2_v.col(i);
        }

        // finding k3 using updated possitions    
        for (int i = 0; i<  (int)  All_Particles.size(); i++)
        {
            k3_r.col(i) = All_Particles[i].v_*dt;
            k3_v.col(i) = total_force(i, interactions, time_dependance, t + dt/2)/All_Particles[i].m_*dt;

        }

        // uppdating the possition and velocity starting form the initial r and v vectors using k3
        for (int i = 0; i<  (int)  All_Particles.size(); i++)
        {
            All_Particles[i].r_ = r.col(i) + k3_r.col(i);
            All_Particles[i].v_ = v.col(i) + k3_v.col(i);
        }

        // finding k4 using updated possitions            
        for (int i = 0; i<  (int)  All_Particles.size(); i++)
        {
            k4_r.col(i) = All_Particles[i].v_*dt;
            k4_v.col(i) = total_force(i, interactions, time_dependance, t + dt)/All_Particles[i].m_*dt;
        }

        // finding r and v after the time step dt
        for (int i = 0; i<  (int)  All_Particles.size(); i++)
        {
            All_Particles[i].r_ = r.col(i) + (k1_r.col(i) + 2*k2_r.col(i) + 2*k3_r.col(i) + k4_r.col(i))/6;
            All_Particles[i].v_ = v.col(i) + (k1_v.col(i) + 2*k2_v.col(i) + 2*k3_v.col(i) + k4_v.col(i))/6;
        }
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt, bool interactions, bool time_dependance , double t)
{
    arma::mat a = arma::mat(3, (int)  All_Particles.size());
    arma::vec Old;
    double m;

    // finding the total force on all particles using the r and v vectors from before the time step
    for (int i=0; i<  (int)  All_Particles.size(); i++)
    {
        m = All_Particles[i].m_;
        a.col(i) = total_force(i, interactions, time_dependance, t)/m;
    }

    // update the v and r vectors of the particles
    for (int i=0; i<  (int)  All_Particles.size(); i++)
    {
        All_Particles[i].r_ +=  All_Particles[i].v_ *dt;
        All_Particles[i].v_ += a*dt;
    }
}


// counting number of particles inside the trap
int PenningTrap::count()
{
    double norm;
    int number_of_particles_inside = 0;

    for (int i = 0; i <  (int)  All_Particles.size(); i++)
    {
        norm = arma::norm(All_Particles[i].r_);

        if (norm<d_)
        {
            number_of_particles_inside += 1;
        }
    }
    return number_of_particles_inside;
}