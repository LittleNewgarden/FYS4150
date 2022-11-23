#include <iostream>
#include <vector>
// #include <string>
// #include <math.h>
// #include <iomanip>
// #include <fstream>
#include <armadillo>
// #include <chrono>

using namespace std;
using namespace arma;

#include "utils/mcmc.hpp"




int main()
{
    // variables that needs to be given when the program is run   
    int L ;            //size of the lattice
    double percent ;   // percentage of spin down in the inital lattice
    double T_start ;   // start teperature
    double T_end ;     // end teperature
    int Tn ;           // number of temperatures, If only oen temperature: T_start = T_end, Tn = 1
    int N ;            // number of Monte Carlo cycles
    int full ;         // if full == 1: the value of energy and magnetization will be saved for each MC cycle
    int average ;      // if average == 1: the average value of energy and magnetization for each temperature will be saved
    string name ;      // the name of the file, _average.txt or _all.txt will be added to the end of the names depnding on full and average
    
    cout<<"Latice side lenght :"<<flush ;
    cin >> L ;
    cout<<"percentage of spin -1 :"<<flush ;
    cin >> percent ;
    cout<<"Start temperature :"<<flush ;
    cin >> T_start ;
    cout<<"End temoerature :"<<flush ;
    cin >> T_end ;
    cout<<"Number of temperatures :"<<flush ;
    cin >> Tn ;
    cout<<"number of MC sycles :"<<flush ;
    cin >> N ;
    cout<<"Write all E and M values to file? :"<<flush ;
    cin >> full ;
    cout<<"Write average values of E and M to file? :"<<flush ;
    cin >> average ;
    cout<<"Name of file :"<<flush;
    cin >> name ;

    int NN ;           // number of spins in lattice
    NN = L*L ; 

    arma_rng::set_seed_random() ;                // getting a seed

    mat S0 = mat(L, L, fill::randu) ;            // initializing the lattice and magnitization
    int M0 = initial_S_M0(S0, percent) ;
    vec T_list = linspace(T_start, T_end, Tn) ;  // vector of T values


    mat EM_out ;  //declearing the matrixes uesd to save wanted values
    mat out ; 


    if (average==0 && full==0)  // varning that no values will be saved
    {
        cout << "No file will be written \n\n";
    }

    if (average==1)      // defining the matix used to save the average values in addition to Cv and chi
    {
        out = mat(Tn, 6) ;
    }

    if (full==1)         // defining the matrix used to save the valeus of energy and magnetization ofter each MC cycle for each temperature
    {
        EM_out = mat(N, Tn*2) ;
    }
 
    // finding the time before the parallelized or non-parallelized code
    auto t1 = std::chrono::high_resolution_clock::now();  

    // parallized code

    #ifdef _OPENMP    
    {
        int nThreads = omp_get_max_threads() ;  // geting max number of threads

        if (Tn<nThreads)     // if the number of threads is greater than the number of temperatures, dont use all threads
        {
            nThreads = Tn ;
        }

        omp_set_num_threads(nThreads); // use nTreads number of threads

        #pragma omp parallel for
        for (int i = 0; i < Tn; i++)   // parallelized for loop
        {
            arma_rng::set_seed_random() ;  
            mat S = S0 ;                           // initial state is the same for each temperature
            vector<vector<double *>> Sref ;        // declearing the reference matrix 
            int M_0 = M0 ;                         // magnetization is the same for each temperature
            double E_0 = initial_Sref_E0(S, Sref); // initiallization of Sref so that it points to S, and get initial energy

            mat EM = MCMC(N, NN, T_list[i], S, Sref, E_0, M_0);   // using MCMC to find the energy and magnetization of lattice
            

            //  CALCULATING THE WANTED QUANTITIES

            if (average==1)
            {
                int burn_in = 20000 ;  
                mat EM_burn = EM(span(burn_in, EM.n_rows-1), span()) ;   // removing data from burn-in period

                mat EM2 = EM_burn%EM_burn ;  // energy and magnitization squared

                mat EM_abs_M = EM_burn ; 
                EM_abs_M.col(1) = abs(EM_burn.col(1)) ;  // absolute value of magnetization

            
                mat em_mean = (mean(EM_abs_M, 0))/(NN) ;  // mean energy and magnetization over number of spins
                mat em2_mean = (mean(EM2, 0))/(NN*NN) ;   // mean energy and magnetization over number of spins
    
                out.row(i)[0] = em_mean[0] ;              // saving in matrix
                out.row(i)[1] = em2_mean[0] ;
                out.row(i)[2] = em_mean[1] ;
                out.row(i)[3] = em2_mean[1] ;             
                out.row(i)[4] = NN/T_list[i]/T_list[i]*( em2_mean[0] - em_mean[0]*em_mean[0] ); 
                out.row(i)[5] = NN/T_list[i]*( em2_mean[1] - em_mean[1]*em_mean[1] );  
            }

            if (full==1)  // saving all energy and magnetization values in matrix
            {
                EM_out.col(i*2) = EM.col(0) ;
                EM_out.col(i*2+1) = EM.col(1) ;
            }
        }
    }
    #else
    {
        for (int i = 0; i < Tn; i++)
        {
            mat S = S0 ;                           // initial state is the same for each temperature
            vector<vector<double *>> Sref ;        // declearing the reference matrix 
            int M_0 = M0 ;                         // magnetization is the same for each temperature
            double E_0 = initial_Sref_E0(S, Sref); // initiallization of Sref so that it points to S, and get initial energy

            mat EM = MCMC(N, NN, T_list[i], S, Sref, E_0, M_0);   // using MCMC to find the energy and magnetization of lattice


            //  CALCULATING THE WANTED QUANTITIES

            if (average==1)
            {
                int burn_in = 20000 ;  
                mat EM_burn = EM(span(burn_in, EM.n_rows-1), span()) ;   // removing data from burn-in period

                mat EM2 = EM_burn%EM_burn ;  // energy and magnitization squared

                mat EM_abs_M = EM_burn ; 
                EM_abs_M.col(1) = abs(EM_burn.col(1)) ;  // absolute value of magnetization

            
                mat em_mean = (mean(EM_abs_M, 0))/(NN) ;  // mean energy and magnetization over number of spins
                mat em2_mean = (mean(EM2, 0))/(NN*NN) ;   // mean energy and magnetization over number of spins
    
                out.row(i)[0] = em_mean[0] ;              // saving in matrix
                out.row(i)[1] = em2_mean[0] ;
                out.row(i)[2] = em_mean[1] ;
                out.row(i)[3] = em2_mean[1] ;             
                out.row(i)[4] = NN/T_list[i]/T_list[i]*( em2_mean[0] - em_mean[0]*em_mean[0] ); 
                out.row(i)[5] = NN/T_list[i]*( em2_mean[1] - em_mean[1]*em_mean[1] );  
            }

            if (full==1)  // saving all energy and magnetization values in matrix
            {
                EM_out.col(i*2) = EM.col(0) ;
                EM_out.col(i*2+1) = EM.col(1) ;
            }
        }
    }
    #endif


    auto t2 = std::chrono::high_resolution_clock::now() ;                       // finding the time after 
    double duration_seconds = std::chrono::duration<double>(t2 - t1).count() ;  // calculating the difference
        
    cout << "MCMC time :" << duration_seconds << "\n" ;    // printin the time used for the parallelized or non-parallelized code
    

    if (average==1)  // saving the chosen values
    {    
        out.save(name + "_average.txt", raw_ascii) ;
    }

    if (full==1)
    {    
        EM_out.save(name + "_all.txt", raw_ascii) ; 
    }


    return 0 ;
}

