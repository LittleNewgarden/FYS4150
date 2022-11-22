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

// function used to initiallise the state matirx and return the initial magnetization based on the percentage of spin down
int initial_S_M0(mat& S, double percent);


// function used to initiallise the reference matirx and return the initial energy
double initial_Sref_E0(mat& S, vector<vector<double *>>& Sref);


// Initialising the refence matrix using the state matrix
vector<vector<double *>> ref_mat(mat& S);


// One MC cycle finding the new energy, new magnetization and the new state S
void MC(int NN, double T, mat& S,vector<vector<double *>>& Sref, double& E0, int& M0);


// outer loop of the MCMC method. saves the values after each MC cycle
mat MCMC(int N, int NN, double T, mat& S, vector<vector<double *>>& Sref, double& E0, int& M0);