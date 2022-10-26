#include <iostream>
#include <armadillo>

#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace std;
using namespace arma;




// function designed to return number of particles remaining after itterating over time with RK4 or forward euler, made for task9 
int N_particles_end(int n, PenningTrap A, double dt, bool interactions = true, bool time = false)
{
    // choosing if we are using time dependant potential
    if (time==true)
    {
        double t=0;

        // integration loop
        for (int k = 0; k < n; k++)
        {  
            A.evolve_RK4(dt, interactions,time,t);
            t+=dt;
        }

        return A.count();
    }

    else{
        for (int k = 0; k < n; k++)
        {  
            A.evolve_RK4(dt, interactions);
        }  

        return A.count();
    }
}
        





int main()
{

    // Penning trap paramiters
    double B0 = 9.65e1; // u/(micro s)/e
    double V0 = 2.41e6; // u(micro m)^2/(micro s)^2/e
    double d = 500;  //micro m

    // particle paramiters
    int q = 1; 
    double m = 40.078;

    // used to determine what the evolve_save function will do
    string method1 = "RK4";
    string method2 = "Euler";
    bool interactions_false = false;
    bool interactions_true = true;
    bool r_write_true = true;
    bool r_write_false = false;
    bool v_write_false = false;
    bool v_write_true = true;


    // creating two Penning trap instances in order to generate two datasets to average over.
    PenningTrap task9a = PenningTrap(B0,V0,d);
    PenningTrap task9b = PenningTrap(B0,V0,d);


    // adding 100 particles to the Penning trap objects task9a and task9b
    for (int i = 0; i < 100;i++)
    {
        vec ra = vec(3).randn() * 0.1 * task9a.d_;  // random initial position
        vec va = vec(3).randn() * 0.1 * task9a.d_;  // random initial velocity
        
        task9a.add_particle(Particle(q, m, ra, va));


        vec rb = vec(3).randn() * 0.1 * task9b.d_;  // random initial position
        vec vb = vec(3).randn() * 0.1 * task9b.d_;  // random initial velocity
        
        task9b.add_particle(Particle(q, m, rb, vb));
    }    
    
    
    // list of amplitudes
    vec f = vec{0.1,0.4,0.7};

    // minimum stepsize
    int steps = (int) ((2.5-0.2)/0.02) + 1;
    
    // list range of frequencies
    vec Wv = linspace(0.2,2.5,steps);
    

    mat n_particles_500ms = mat(f.n_rows, Wv.n_rows);

    // integration paramiters
    double t = 500; //micro s
    int n = 40000;
    double dt = t/n;       
    string name;

    // itterating over f and Wv to find number of particles left after 500 micro s for each combination without interactions
    for (int i = 0; i <  (int)  f.n_rows; i++)
    {
        task9a.f_ = f[i]; 
        for (int j = 0; j<  (int)  Wv.n_rows; j++)
        {
            task9a.Wv_ = Wv[j];

            n_particles_500ms.row(i)[j] = N_particles_end( n, task9a, dt, interactions_false, true);
        } 
    }

    name = "100_particles_500ms_interactions=0.txt";
    n_particles_500ms.save(name, raw_ascii);



    // zooming in to on area of Wv

    int w_steps = 50;

    vec Wv_zoomed = linspace(2.12,2.22,w_steps);

    vec n_particles_500ms_zoomed_a = vec(Wv_zoomed.n_rows);
    vec n_particles_500ms_zoomed_b = vec(Wv_zoomed.n_rows);

    vec n_particles_500ms_zoomed_interactions_a = vec(Wv_zoomed.n_rows);
    vec n_particles_500ms_zoomed_interactions_b = vec(Wv_zoomed.n_rows);

    task9a.f_ = f[0];
    task9b.f_ = f[0];

    // itterating over Wv to find number of particles left after 500 micro s with and without interactions
    //  for two different instances of Penning trap
    for (int j = 0; j < (int) Wv_zoomed.n_rows; j++)
    {
        task9a.Wv_ = Wv_zoomed[j];
        task9b.Wv_ = Wv_zoomed[j];

        n_particles_500ms_zoomed_a[j] = N_particles_end( n, task9a, dt, interactions_false, true);
        n_particles_500ms_zoomed_b[j] = N_particles_end( n, task9b, dt, interactions_false, true);

        n_particles_500ms_zoomed_interactions_a[j] = N_particles_end( n, task9a, dt, interactions_true, true);
        n_particles_500ms_zoomed_interactions_b[j] = N_particles_end( n, task9b, dt, interactions_true, true);

        cout << j+1 << "\n\n";
    } 

    name = "100_particles_500ms_interactions=0_zoomed_inn_n=" + to_string(w_steps) + "_a.txt";
    n_particles_500ms_zoomed_a.save(name, raw_ascii);

    name = "100_particles_500ms_interactions=0_zoomed_inn_n=" + to_string(w_steps) + "_b.txt";
    n_particles_500ms_zoomed_b.save(name, raw_ascii);

    name = "100_particles_500ms_interactions=1_zoomed_inn_n=" + to_string(w_steps) + "_a.txt";
    n_particles_500ms_zoomed_interactions_a.save(name, raw_ascii);

    name = "100_particles_500ms_interactions=1_zoomed_inn_n=" + to_string(w_steps) + "_b.txt";
    n_particles_500ms_zoomed_interactions_b.save(name, raw_ascii);

    
    return 0;
}