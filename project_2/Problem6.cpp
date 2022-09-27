#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <chrono>

using namespace std;
using namespace arma;

void writefile(vec x, vec y, string name);
mat sym_tri_diag_A(int N, double d, double a);


double off_diag(const mat& A, int& k, int& l);

void Jacobi_rotation_algorithem(mat A, double eps, vec& eigenvalues, mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged);

void analytical_sym_tri_N6(int N, double d, double a, mat& v, vec& lam);

int main(){

    double epsilon = pow(10,-8); // convergens limit


    vector<int> N_list = vector<int>{11,101}; 
    for (int i = 0; i<2; i++){
        int N = N_list[i];     // defining size of matrix
        // int N = 101;



        double h = (1.)/(N-1);  // defining step length
        double d = 2./(h*h);    // defining diagonals
        double a = -1./(h*h);

        // creating tri diagonal symetric matrix
        // not solving for edges, therfore N-2
        mat A = sym_tri_diag_A(N-2, d, a);

        vec eigenvalues = vec(N-2);
        mat eigenvectors = mat(N-2,N-2);

        double maxiter = pow(10,7);  // maximum number of allowed iterations of the jacobi rotation 
        int iterations;  // used to count itterations of jacobi rotations
        bool converged; // used to check if solution converged

        // function finding eigen values and eigen vectors of matrix A using jacobi rotations
        Jacobi_rotation_algorithem(A, epsilon, eigenvalues, 
                eigenvectors, maxiter, iterations, converged);


        // finfing the idices of the eigenvalues in rising order
        uvec index = sort_index(eigenvalues); 

        // we want only the three smallest eigen vectors
        // and now have size N since we want to add endpoints
        mat eigen_vec_full = mat(N,3).fill(0.);


        // inptting the vanted eigen vectors into eigen_vec_full
        for (int i = 0; i< (int) eigen_vec_full.n_rows-2; i++){
            for (int j = 0; j< (int) eigen_vec_full.n_cols; j++){
                eigen_vec_full(i+1,j) = eigenvectors(i,index(j));

            }
        }

        // saving the eigen vectors to a file
        string name = "prob6_" +  std::to_string(N) + ".txt";
        eigen_vec_full.save(name,raw_ascii);


        // now finding eigen values and eigen vectors for the same matrix using the analytical method
        mat analytical_eigenvectors = mat(N-2,N-2);
        vec analytical_eigenvalues = vec(N-2);

        analytical_sym_tri_N6(N-2, d, a, analytical_eigenvectors, analytical_eigenvalues);


        // finfing the idices of the eigenvalues in rising order
        index = sort_index(analytical_eigenvalues);

        // we want only the three smallest eigen vectors
        // and now have size N since we want to add endpoints
        mat analytical_eigen_vec_full = mat(N,3).fill(0.);

        // inptting the vanted eigen vectors into eigen_vec_full
        for (int i = 0; i< (int) eigen_vec_full.n_rows-2; i++){
            for (int j = 0; j< (int) eigen_vec_full.n_cols; j++){
                analytical_eigen_vec_full(i+1,j) = analytical_eigenvectors(i,index(j));
            }
        }

        // saving the eigen vectors to a file
        name = "prob6_analytical_" +  std::to_string(N) + ".txt";
        analytical_eigen_vec_full.save(name,raw_ascii);
    }
    return 0;
}













    // function finding eigen values and eigen vectors of matrix A
void Jacobi_rotation_algorithem(mat A, double eps, vec& eigenvalues, mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged){
    
    // indices of highest of diagonal value
    int k; 
    int l;

    iterations = 0; // number of iterations start at 0 

    mat R = eye(A.n_rows,A.n_cols); // R stsrts of as identity matrix

    double max = off_diag(A, k, l); //findign highest value and indices of A
    
    double t, c, s; // defining tan, cos and s

    while (abs(max)>eps){ // while loop going untill highest of diagonal value is less than chosen convergence limit

        double tau = (A(l,l) - A(k,k))/(2*A(k,l));

        if (tau>=0){  // chosing the smallest rotation
            t = 1./(tau + sqrt(1.+pow(tau,2)));
        }

        if (tau<0){
            t = -1./(-tau + sqrt(1.+pow(tau,2)));
        }

        c = 1/sqrt(1.+pow(t,2));
        s = c*t;


        // updating values in A
        //saving A(k,k) since it is needed after its been updated
        double a_kk = A(k,k) ;

        A(k,k) = A(k,k)*pow(c,2.) - 2*A(k,l)*c*s + A(l,l)*pow(s,2.);
        A(l,l) = A(l,l)*pow(c,2.) + 2*A(k,l)*c*s + a_kk*pow(s,2.);

   
        A(k,l) = 0.;
        A(l,k) = 0.;

        for (int i = 0; i< (int) A.n_cols; i++){
            if (i!=k && i!=l){

                double a_ik = A(i,k);

                A(i,k)= A(i,k)*c-A(i,l)*s;
                A(k,i) = A(i,k);

                A(i,l) = A(i,l)*c + a_ik*s;
                A(l,i) = A(i,l);
                
            }

            double r_ik = R(i,k);

            R(i,k) = R(i,k)*c - R(i,l)*s;
            R(i,l) = R(i,l)*c + r_ik*s;

        }

        // one itteration compleat
        iterations = iterations + 1 ;

        if (iterations>=maxiter){ // stop loop if it runs to long before convergence
            break ;
        }
    
        // findign new max value after A has been updated
        max = off_diag(A, k, l);

    }


    // testing if the algorithem converged
    if (iterations>=maxiter){
        converged = false;
    } else{
        converged = true;
    }

    // saves the eigen values and eigen vectors
    for (int i = 0; i< (int) A.n_cols; i++){
        eigenvalues(i) = A(i,i);
        eigenvectors.col(i) = normalise(R.col(i));
    }

}





// function for finding the largest off diagonal value and indices
double off_diag(const mat& A, int& k, int& l){

    double max = 0.; // inital value of max

    for (int i = 0; i < (int) A.n_rows-1; i++){    // nested loops, iterating only over the upper triangular part
        for (int j = i+1; j < (int) A.n_cols; j++){        // not including the main diagonal

            if (abs(A(i,j)) > abs(max)){ // testing if the element (i,j) is greater than max
                max = A(i,j);  // saving max, k and l
                k = j;
                l = i;
            }
        }
    }

    return max;
}






// function for creating tri diagonal symetric matrix of sice N
mat sym_tri_diag_A(int N, double d, double a){

    mat A = mat(N,N).fill(0.);
    
    for (int i = 0; i< N-1; i++){ // filling the diagonals with values a and d

        A(i,i) = d;
        A(i,i+1) = a;
        A(i+1,i) = a;

    }

    A(N-1,N-1) = d; //filling last value of main diagonal

    return A;

}







// function for analytically finding eigen vectors and eigen values 
void analytical_sym_tri_N6(int N, double d, double a, mat& v, vec& lam){

    for (int i=1; i<N+1; i++){  //finding eigen values and eigen values analytically

        lam(i-1) = d + 2*a*cos(i*datum::pi/(N+1));
        
        for (int j=1; j<N+1; j++){
        
            v(i-1,j-1) = sin(i*j*datum::pi/(N+1)) ;
        
        }
    }
    v = normalise(v.t());
}