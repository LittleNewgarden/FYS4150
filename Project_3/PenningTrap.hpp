#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <iostream>
#include <armadillo>

#include "Particle.hpp"



class PenningTrap
{
    public:
        // initializing variables 
        double B0_; double V0_; double d_; double f_; double Wv_;

        // vector containign particles
        std::vector<Particle> All_Particles;        

        //constant
        double Ke_ = 1.38935333e5;

        // function defining the Penning trap system
        PenningTrap(double B0, double V0, double d, double f=0, double Wv=0);

        // Add a particle to the trap
        void add_particle(Particle p);

        // Add all particles to the trap
        void add_all_particles(std::vector<Particle> p );

        // External electric field at point r=(x,y,z)
        arma::vec external_E_field(arma::vec r, bool time_dependance, double t);  

        // External magnetic field at point r=(x,y,z)
        arma::vec external_B_field(arma::vec r);  

        // Force on particle_i from particle_j
        arma::vec force_particle(int i, int j);

        // The total force on particle_i from the external fields
        arma::vec total_force_external(int i, bool time_dependance, double t);

        // The total force on particle_i from the other particles
        arma::vec total_force_particles(int i);

        // The total force on particle_i from both external fields and other particles
        arma::vec total_force(int i, bool interactions, bool time_dependance, double t);
        
        // Evolve the system one time step (dt) using Runge-Kutta 4th order
        void evolve_RK4(double dt, bool interactions, bool time_dependance=false, double t=0);

        // Evolve the system one time step (dt) using Forward Euler
        void evolve_forward_Euler(double dt, bool interactions, bool time_dependance=false, double t=0);

        // counting number of particles inside the trap
        int count();



};





















#endif