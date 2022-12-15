#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <assert.h>


// function used to find index k in vectorised matrix u_ij. 
// We dont use this since armadillo does it automaticly by caling u_ij(k), but creating the function was one of the tasks, so here it is
int indx(int i, int j, arma::mat A) ;


// Function that creates the A and B matricies 
void AB_sparce(int M, double h, double dt, arma::mat V, arma::sp_cx_mat& A, arma::sp_cx_mat& B) ;


// initialises the wave packet u
void initialise_u0(arma::cx_mat& u0, int M, double xc, double yc, double sigma_x, double sigma_y, double px, double py) ;


// creates a solid wall potential
void wall_V(arma::mat& V, int M, double v0) ;

// creates a single slit potential
void single_slit_V(arma::mat& V, int M, double v0) ;

// creates a double slit potential
void double_slit_V(arma::mat& V, int M, double v0) ;

// creates a triple slit potential
void triple_slit_V(arma::mat& V, int M, double v0) ;

