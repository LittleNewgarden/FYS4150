#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <iostream>
#include <armadillo>



class Particle
{
    public:
        // initializing variables
        int q_; double m_; arma::vec r_; arma::vec v_;

        // function that defines particle properties
        Particle(int q, double m, arma::vec r, arma::vec v);
};


#endif