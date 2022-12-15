#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <assert.h>

#include "include/all_functions.hpp"


using namespace std;
using namespace arma;


int main(int argc, const char* argv[])
{

    if (argc != 14)
    {
        cout << "14 arguments are requierd, you gave : " + to_string(argc) << "\n" ;
        exit(0) ;
    }

    // innput arguments 
    double h = atof(argv[1]) ;
    double dt = atof(argv[2]) ;
    double T = atof(argv[3]) ;
    double xc = atof(argv[4]) ;
    double yc = atof(argv[5]) ;
    double sigma_x = atof(argv[6]) ;
    double sigma_y = atof(argv[7]) ;
    double px = atof(argv[8]) ;
    double py = atof(argv[9]) ;
    double v0 = atof(argv[10]) ;
    string name = argv[11] ;
    string slit = argv[12] ;
    string save_u = argv[13] ;

    // number of points along x and y in simulation box
    int M = 1/h + 1;
    
    // declear V
    mat V ;

    // creating the wanted potential 
    if (slit == "single")
    {
        single_slit_V(V, M, v0) ;
    }
    if (slit == "double")
    {
        double_slit_V(V, M, v0) ;
    }
    if (slit == "triple")
    {
        triple_slit_V(V, M, v0) ;
    }
    if (slit == "wall")
    {
        wall_V(V, M, v0) ;
    }
    if (slit == "no_slit")
    {
        V = mat(M-2, M-2).fill(0) ;
    }


    // initialising the wave packet
    cx_mat u0_boundary;
    initialise_u0(u0_boundary, M, xc, yc, sigma_x, sigma_y, px, py);

    // defining the inner wave packet
    arma::cx_mat u0 = u0_boundary( arma::span(1,M-2), arma::span(1,M-2));


    // creating the matrices A and B
    sp_cx_mat A;
    sp_cx_mat B;

    AB_sparce(M, h, dt, V, A, B);

    // number of time steps 
    int Tn = T/dt ;
    
    // creating the necesary cubes used to compute and save P and U
    cube P = cube(M, M, Tn+1) ;
    cx_cube U = cx_cube(M, M, Tn+1) ;
    cx_cube U_temp = cx_cube(M-2, M-2, Tn+1) ;
    cx_cube P_temp = cx_cube(M-2, M-2, Tn+1) ;

    // finding the probability distribution p 
    cx_mat u0_conj = conj( u0 ) ;
    cx_mat p0 = u0_conj % u0 ;

    // initial values
    P_temp.slice(0) = p0 ;
    U_temp.slice(0) = u0 ;


    // used to find the runtime
    auto t1 = std::chrono::high_resolution_clock::now();

    // iterating in time 
    for (int i = 0; i < Tn; i++)
    {           
        // vectorising matirx u_ij
        cx_vec ui_vec = vectorise( U_temp.slice(i) ) ;
        
        // finding b = B*u
        cx_mat b = B*ui_vec ;
        
        // solving A * u^(n+1) = b
        cx_mat solved = spsolve(A, b) ;

        // calculating p 
        cx_mat solved_conj = conj( solved ) ;
        cx_mat p = solved_conj % solved ;
        
        // saving p and u in cubes 
        for (int j = 0; j < (int) solved.size(); j++)
        {
            U_temp.slice(i+1)[j] = solved[j] ;
            P_temp.slice(i+1)[j] = p[j] ;
        }
    }

    // end timer
    auto t2 = std::chrono::high_resolution_clock::now();
    double duration_seconds = std::chrono::duration<double>(t2 - t1).count();
    
    // print time 
    cout << "MCMC time :" << duration_seconds << "\n";

    // impliment boundary conditions for each time 
    P(span(1, M-2), span(1, M-2), span::all) = real(P_temp) ;
    U(span(1, M-2), span(1, M-2), span::all) = U_temp ;

    // addition to names of real and imaginary when saving U
    string name1 = "_real.txt" ;
    string name2 = "_imag.txt" ;
    
    // when saving must use raw_ascii since we did not manage to get the numpy armadillo binary function to work

    // save P
    P.save(name + ".txt", raw_ascii) ;

    // save U
    if (save_u == "true")
    {
        cube Usave_r = (real(U)) ;
        Usave_r.save(name + name1, raw_ascii) ;
    
        cube Usave_i = (imag(U)) ;
        Usave_i.save(name + name2, raw_ascii) ;
    }
    
    return 0;
}





// task 7 :

// 0.005 0.000025 0.008 0.25 0.5 0.05 0.05 200 0 0 task7_1 no_slit false 
// 0.005 0.000025 0.008 0.25 0.5 0.05 0.1 200 0 10000000000 task7_2 double false

// task 8 :

// 0.005 0.000025 0.002 0.25 0.5 0.05 0.2 200 0 10000000000 task8 double true

// task 9 :

// 0.005 0.000025 0.002 0.25 0.5 0.05 0.2 200 0 10000000000 task9_single single false
// 0.005 0.000025 0.002 0.25 0.5 0.05 0.2 200 0 10000000000 task9_double double false
// 0.005 0.000025 0.002 0.25 0.5 0.05 0.2 200 0 10000000000 task9_triple triple false
