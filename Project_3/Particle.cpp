#include "Particle.hpp"



// function that defines particle properties
Particle::Particle(int q, double m, arma::vec r, arma::vec v)
{
    q_ = q;
    m_ = m;
    r_ = r;
    v_ = v;
}
