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
    // n = 4000 is used for the wz plot


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





    // TASK 9 

    PenningTrap task9 = PenningTrap(B0,V0,d);

    // adding 100 particles to the Penning trap object task9
    for (int i = 0; i < 100;i++)
    {
        vec r = vec(3).randn() * 0.1 * task9.d_;  // random initial position
        vec v = vec(3).randn() * 0.1 * task9.d_;  // random initial velocity
        
        task9.add_particle(Particle(q, m, r, v));
    }    
    
    
    // list of amplitudes
    vec f = vec{0.1,0.4,0.7};

    // minimum stepsize
    int steps = (int) ((2.5-0.2)/0.02) + 1;
    
    // list range of frequencies
    vec Wv = linspace(0.2,2.5,steps);
    

    mat n_particles_500ms = mat(f.n_rows, Wv.n_rows);

    // integration paramiters
    t = 500; //micro s
    int n = 40000;
    dt = t/n;       
    string name;

    // itterating over f and Wv to find number of particles left after 500 micro s for each combination without interactions
    for (int i = 0; i <  (int)  f.n_rows; i++)
    {
        task9.f_ = f[i]; 
        for (int j = 0; j<  (int)  Wv.n_rows; j++)
        {
            task9.Wv_ = Wv[j];

            n_particles_500ms.row(i)[j] = N_particles_end( n, task9, dt, interactions_false, true);
        } 
    }

    name = "100_particles_500ms_no_interactions.txt";
    n_particles_500ms.save(name, raw_ascii);



    // zooming in to on area of Wv

    int w_steps = 50;

    vec Wv_zoomed = linspace(2.12,2.22,w_steps);
    vec n_particles_500ms_zoomed = vec(Wv_zoomed.n_rows);


    vec n_particles_500ms_zoomed_interactions = vec(Wv_zoomed.n_rows);

    // itterating over f and Wv to find number of particles left after 500 micro s for each combination without and without interactions
    for (int i = 0; i < (int) 1; i++)
    {
        task9.f_ = f[i];
        for (int j = 0; j < (int) Wv_zoomed.n_rows; j++)
        {
            task9.Wv_ = Wv_zoomed[j];
        
            n_particles_500ms_zoomed[j] = N_particles_end( n, task9, dt, interactions_false, true);
            n_particles_500ms_zoomed_interactions[j] = N_particles_end( n, task9, dt, interactions_true, true);
            cout << j+1 << "\n\n";
        } 
    }

    name = "100_particles_500ms_no_interactions_zoomed_inn_n=" + to_string(w_steps) + ".txt";
    n_particles_500ms_zoomed.save(name, raw_ascii);

    name = "100_particles_500ms_no_interactions_zoomed_inn_n=" + to_string(w_steps) + "_interactions.txt";
    n_particles_500ms_zoomed_interactions.save(name, raw_ascii);


    
    return 0;
}