#include <iostream>
#include <armadillo>

#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace std;
using namespace arma;


// function designed to write out the r and v vectors after itterating over time with RK4 or forward euler, made for task8 
void evolve_save_r_v(int n, PenningTrap A, double dt, string met = "RK4", bool interactions = true, bool r_write = true, bool v_write = false){
    int N = A.All_Particles.size();
    cube r = cube(n+1,3,N);
    cube v = cube(n+1,3,N);
    double t;

    // choice of integration method
    if (met == "RK4")
    {
        //  integration loop
        for (int i = 0; i < n; i++)
        {  
            // loop over particles
            for (int j = 0; j < N; j++)
            { 
                r.slice(j).row(i) = (A.All_Particles[j].r_).t();
                v.slice(j).row(i) = (A.All_Particles[j].v_).t();    
            }

            t+=dt;
            A.evolve_RK4(dt, interactions);
        }

        // loop over particles to add last step
        for (int j = 0; j < N; j++)
        { 
            r.slice(j).row(n) = (A.All_Particles[j].r_).t();
            v.slice(j).row(n) = (A.All_Particles[j].v_).t();  
        }


        string name_r;
        string name_v;

        // wirte files for r and v vectors
        if (r_write==true)
        {
            for (int j = 0; j < (int) A.All_Particles.size(); j++)
            {
                name_r = "r_n_iterations=" + std::to_string(n) + "_method=" + met + "_number_particles=" + std::to_string(N)
                + "_particle=" +  std::to_string(j) + "_interactions=" + std::to_string(interactions) + ".txt";
                r.slice(j).save(name_r,raw_ascii);
            }
        }
        
        if (v_write ==true)
        {
            for (int j = 0; j < (int) A.All_Particles.size(); j++)
            {
                name_v = "v_n_iterations=" + std::to_string(n) + "_method=" + met + "_number_particles=" + std::to_string(N)
                + "_particle=" +  std::to_string(j) + "_interactions=" + std::to_string(interactions) + ".txt";
                v.slice(j).save(name_v,raw_ascii);
            }     
        }
    }

    if (met == "Euler")
    {
        
        //  integration loop
        for (int i = 0; i < n; i++)
        {  
            // loop over particles
            for (int j = 0; j < N; j++)
            { 
                r.slice(j).row(i) = (A.All_Particles[j].r_).t();
                v.slice(j).row(i) = (A.All_Particles[j].v_).t();    
            }

            t+=dt;
            A.evolve_forward_Euler(dt, interactions);
        }

        // loop over particles to add last step
        for (int j = 0; j < N; j++)
        { 
            r.slice(j).row(n) = (A.All_Particles[j].r_).t();
            v.slice(j).row(n) = (A.All_Particles[j].v_).t();    
        }
    
        string name_r;
        string name_v;

        // wirte files for r and v vectors
        if (r_write==true)
        {
            for (int j = 0; j < (int) A.All_Particles.size(); j++)
            {
                name_r = "r_n_iterations=" + std::to_string(n) + "_method=" + met + "_number_particles=" + std::to_string(N)
                + "_particle=" +  std::to_string(j) + "_interactions=" + std::to_string(interactions) + ".txt";
                r.slice(j).save(name_r,raw_ascii);
            }
        }

        if (v_write ==true)
        {
            for (int j = 0; j < (int) A.All_Particles.size(); j++)
            {
                name_v = "v_n_iterations=" + std::to_string(n) + "_method=" + met + "_number_particles=" + std::to_string(N)
                + "_particle=" +  std::to_string(j) + "_interactions=" + std::to_string(interactions) + ".txt";
                v.slice(j).save(name_v,raw_ascii);
            }     
        }
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
    vec r = vec{20, 0.0, 20};
    vec v = vec{0.0,25,0.0};
    vec r1 = vec{25, 25, 0};
    vec v1 = vec{0,40,5};

    // defining particle 1
    Particle Ca =  Particle(q, m, r, v);

    // defining particle 2
    Particle Ca2 =  Particle(q, m, r1, v1);
    
    // Penning trap object used for testing
    PenningTrap test1 = PenningTrap(B0,V0,d);
   
    // adding particle 1 to test1
    test1.add_particle(Ca);
 
    
    double t = 50; // micro s
    double dt; 

    // used to determine what the evolve_save function will do
    string method1 = "RK4";
    string method2 = "Euler";
    bool interactions_false = false;
    bool interactions_true = true;
    bool r_write_true = true;
    bool r_write_false = false;
    bool v_write_false = false;
    bool v_write_true = true;



    // single particle n=[4000, 8000, 16000, 32000], 50 micro s, euler and RK4


    vec n_l = vec{4000, 8000, 16000, 32000};

    // saving r at each time step using different time steps and different integration methods
    for (int i=0; i < (int) n_l.size(); i++)
    {
        dt = t/n_l[i]; 
        evolve_save_r_v(n_l[i],test1, dt, method1, interactions_true, r_write_true, v_write_false);
        evolve_save_r_v(n_l[i],test1, dt, method2, interactions_true, r_write_true, v_write_false);
    }


    // two particles, n= 4000 with and without interactions, RK4


    dt = t/n_l[0];
    
    // adding secound particle
    test1.add_particle(Ca2);

    // saving r and v
    evolve_save_r_v(n_l[0],test1, dt, method1, interactions_true, r_write_true, v_write_true);

    evolve_save_r_v(n_l[0],test1, dt, method1, interactions_false, r_write_true, v_write_true);



   
    return 0;
}