#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <armadillo>



arma::mat sym_tri_diag_A(int N,double d, double a);



int main(){

    // defining size of matrix
    int N = 6;

    double h = 1./(N-1);    // defining step length

    double d = 2./(h*h);    // defining diagonals
    double a = -1./(h*h);

    arma::mat A = sym_tri_diag_A(N,d,a);    // creating tri diagonal symetric matrix


    arma::vec eigval;   
    arma::mat eigvec;

    arma::eig_sym(eigval, eigvec, A);  //finfing eigen value and eigen vector using armadillo


    std::cout << "armadillo:\n\n"; //printing eigen values and vectors found using armadillo

    std::cout << eigvec ; 
    
    std::cout << "\n\n";
    
    std::cout << eigval  ;


    arma::vec lam = arma::vec(N).fill(0.);  
    arma::mat v = arma::mat(N,N).fill(0.);

    for (int i=1; i<N+1; i++){      //finding eigen values and eigen values analytically
        
        lam(i-1) = d + 2*a*cos(i*arma::datum::pi/(N+1));        
        
        for (int j=1; j<N+1; j++){
        
            v(i-1,j-1) = sin(i*j*arma::datum::pi/(N+1)) ;
        
        }
    }
    v = arma::normalise(v.t()); 


    std::cout << "\n\n";

    std::cout << "analytically:\n\n";  // printing eigen values and eigen vectors found analytically

    std::cout << (v);

    std::cout << "\n\n";
    
    std::cout << lam;


    return 0;
}














// function for creating tri diagonal symetric matrix of sice N
arma::mat sym_tri_diag_A(int N, double d, double a){

    arma::mat A = arma::mat(N,N).fill(0.);
    
    for (int i = 0; i< N-1; i++){ // filling the diagonals with values a and d

        A(i,i) = d;
        A(i,i+1) = a;
        A(i+1,i) = a;

    }

    A(N-1,N-1) = d; //filling last value of main diagonal

    return A;

}