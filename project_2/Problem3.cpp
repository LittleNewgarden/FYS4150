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



double off_diag(const mat& A, int& k, int& l); 


int main(){
    // defining size of matrix
    int N = 4;

    mat A = eye(N,N);  // creating identity matrix
    vec d = vec{0.5,-0.7,-0.7,0.5}; // the values of the other diagonal

    for (int i = 1; i<=N; i++){ // compleating the A matrix
        A(N-i,i-1) = d(i-1);  

    }
    cout << A << "\n\n";  // printing the A matrix

    int k; 
    int l;


    double  max_val = off_diag(A, k, l);  // finding the maximum off diagonal value and indices

    cout << max_val << ", " << k << ", " << l << "\n\n";

    return 0;
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