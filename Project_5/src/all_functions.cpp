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
int indx(int i, int j, arma::mat A)
{
    int k_indx;
    k_indx = j*A.n_rows + i;
    
    return k_indx;
}



// Function that creates the A and B matricies 
void AB_sparce(int M, double h, double dt, arma::mat V, arma::sp_cx_mat& A, arma::sp_cx_mat& B)
{
    // define constants idt and r
    arma::cx_double i_dt = dt*1j;

    arma::cx_double r;
    r = dt*1j/(2*h*h);
    


    //Pick out the indices i,j from V and transform to k, is commented out for previousl mentioned reasons

    // int q = 0 ;
    // arma::vec k = arma::vec((M-2)*(M-2),arma::fill::zeros); 
    // for (int i = 0; i < (int) V.n_rows; i++){
    //     for (int j=0; j < (int) V.n_cols; j++){
    //         k(q) = indx(j, i, V);
    //         q+=1 ;
    //     }
    // }



    // fill the a and b
    arma::cx_vec ak = arma::cx_vec((M-2)*(M-2)).fill(1. + 4.*r); 
    arma::cx_vec bk = arma::cx_vec((M-2)*(M-2)).fill(1. - 4.*r); 

    arma::vec v = arma::vectorise  (V) ;
    for (int i = 0; i< (M-2)*(M-2); i++)
    {   
        arma::cx_double c = i_dt/2.*v(i) ;  
        ak[i] +=  c;
        bk[i] += -c;
    }


    //Making A and B matrices
    A = arma::sp_cx_mat((M-2)*(M-2), (M-2)*(M-2));
    B = arma::sp_cx_mat((M-2)*(M-2), (M-2)*(M-2));

    // //Filling the main diagonals
    A.diag(0) = ak;
    B.diag(0) = bk;

    //Filling the third super-and subdiagonals
    int diag = M-2 ;

    A.diag(diag).fill(-r);
    A.diag(-diag).fill(-r);
    // std::cout << A << "5\n\n";
    B.diag(diag).fill(r);
    B.diag(-diag).fill(r);

    // //First super-and subdiagonal of A and B that skips every third
    for (int i = 0; i<(M-2)*(M-2)-1; i++){

        if((i+1)%(M-2) == 0){
            continue;
        }

        else{
            A(i,i+1) = -r;
            A(i+1, i) = -r;

            B(i,i+1) = r;
            B(i+1, i) = r;
        }
    }
}



// initialises the wave packet u
void initialise_u0(arma::cx_mat& u0, int M, double xc, double yc, double sigma_x, double sigma_y, double px, double py)
{
    arma::vec x = arma::linspace(0, 1, M);
    arma::vec y = arma::linspace(0, 1, M);

    u0 = arma::cx_mat(M, M).fill(0);

    // filling the initial wave packet u_ij
    for (int i=0; i<M; i++) 
    {
        for (int j=0; j<M; j++)
        {
            double xi_c = x[i] - xc;
            double yi_c = y[j] - yc;
            double top_x = xi_c*xi_c;
            double bottom_x = 2*sigma_x*sigma_x;
            double top_y = yi_c*yi_c;
            double bottom_y = 2*sigma_y*sigma_y;

            std::complex<double> exponent = - top_x/bottom_x - top_y/bottom_y + 1j*px*( xi_c ) + 1j*py*( yi_c );
            u0.col(i)[j] = std::exp( exponent  );
        }
    }

    // defining boundarys
    for (int i=0; i<M; i++)
    {
        u0.row(i)[0] = 0;
        u0.row(i)[M-1] = 0;
        u0.col(i)[0] = 0;
        u0.col(i)[M-1] = 0;
    }
    
    // normalising the wave packet
    arma::cx_mat u0_conj = arma::conj(u0);
    arma::cx_mat u02 = u0_conj % u0;
    
    std::complex<double> sum = arma::accu(u02);

    u0 = u0/std::sqrt(sum);
}



// creates a solid wall potential
void wall_V(arma::mat& V, int M, double v0)
{
    V = arma::mat(M-2, M-2).fill(0);

    // defining x and y axis together as they are the same
    arma::vec xy = arma::linspace(0,1,M)(arma::span(1,M-2)) ;

    arma::vec xi = arma::vec(M-2) ;

    //  finding center index using indx_min
    for (int i = 0; i < M-2; i++)
    {
        xi[i] = std::abs(xy[i] - 0.5) ;
    }    

    int x_i = arma::index_min(xi) ;

    // defining wall thickness
    int wall_thickness = 0.02/(1./M) ;

    // creating wall
    for (int i=0; i< wall_thickness+1; i++)
    {
        V.col(x_i + i).fill(v0) ;
    }
}




// creates a single slit potential
void single_slit_V(arma::mat& V, int M, double v0) // , W_thicknes, s_lenght, W_x
{
    V = arma::mat(M-2, M-2).fill(0);
    
    // defining x and y axis together as they are the same
    arma::vec xy = arma::linspace(0,1,M)(arma::span(1,M-2)) ;

    arma::vec xi = arma::vec(M-2) ;
    arma::vec yi = arma::vec(M-2) ;

    arma::vec yi_edge_max = arma::vec(M-2) ;

    // finding the indix where used to define slit, min/max may be named wrong
    for (int i = 0; i < M-2; i++)
    {
        xi[i] = std::abs(xy[i] - 0.5) ;
        yi[i] = std::abs(xy[i] - 0.5) ;

        yi_edge_max[i] = std::abs(xy[i] - 0.5 - 0.025) ;
 
    }    

    int x_i = arma::index_min(xi) ;
    int y_i = arma::index_min(yi) ;

    int y_e_max = arma::index_min(yi_edge_max) ;


    // defining wall thickness   
    int wall_thickness = 0.02/(1./M) ;

    // creating wall with slit length 11
    for (int i=0; i< wall_thickness+1; i++)
    {
        V.col(x_i + i).fill(v0) ;

        for (int j = y_i + (y_i - y_e_max); j < y_e_max + 1; j++)
        {
            V.col(x_i+i)[j] = 0 ;
        }
    }
}



// creates a double slit potential
void double_slit_V(arma::mat& V, int M, double v0)
{
    V = arma::mat(M-2, M-2).fill(0);

    // defining x and y axis together as they are the same
    arma::vec xy = arma::linspace(0,1,M)(arma::span(1,M-2)) ;
    // std::cout << xy << " " << xy.n_rows << " \n\n";
    arma::vec xi = arma::vec(M-2) ;
    arma::vec yi = arma::vec(M-2) ;

    arma::vec yi_center_edge_max = arma::vec(M-2) ;
    arma::vec yi_slit_edge_max = arma::vec(M-2) ;

    // finding the indicies where used to define the seperation wall and slits, min/max may be named wrong
    for (int i = 0; i < M-2; i++)
    {
        xi[i] = std::abs(xy[i] - 0.5) ;
        yi[i] = std::abs(xy[i] - 0.5) ;

        yi_center_edge_max[i] = std::abs(xy[i] - 0.5 - 0.025) ;
        yi_slit_edge_max[i] = std::abs(xy[i] - 0.5 - 0.025 - 0.05) ;
 
    }    

    int x_i = arma::index_min(xi) ;
    int y_i = arma::index_min(yi) ;

    int y_c_e_max = arma::index_min(yi_center_edge_max) ;
    int y_s_e_max = arma::index_min(yi_slit_edge_max) ;

    // defining wall thickness   
    int wall_thickness = 0.02/(1./M) ;

    // creating wall with slits, center wall length 11, slits length 10
    for (int i=0; i< wall_thickness+1; i++)
    {
        V.col(x_i + i).fill(v0) ;

        for (int j = y_c_e_max+1; j < y_s_e_max+1; j++)
        {
            V.col(x_i+i)[j] = 0 ;
            V.col(x_i+i)[y_i + (y_i - j)] = 0 ;
        }
    }
}



// creates a triple slit potential
void triple_slit_V(arma::mat& V, int M, double v0)
{
    V = arma::mat(M-2, M-2).fill(0);

    // defining x and y axis together as they are the same
    arma::vec xy = arma::linspace(0,1,M)(arma::span(1,M-2)) ;

    arma::vec xi = arma::vec(M-2) ;
    arma::vec yi = arma::vec(M-2) ;

    arma::vec yi_center_edge_min = arma::vec(M-2) ;
    arma::vec yi_slit_edge_min = arma::vec(M-2) ;
    arma::vec yi_slit_edge_max = arma::vec(M-2) ;

    // finding the indicies where used to define the seperation walls and slits, min/max may be named wrong
    for (int i = 0; i < M-2; i++)
    {
        xi[i] = std::abs(xy[i] - 0.5) ;
        yi[i] = std::abs(xy[i] - 0.5) ;

        yi_center_edge_min[i] = std::abs(xy[i] - 0.5 + 0.025) ;
        yi_slit_edge_max[i] = std::abs(xy[i] - 0.5 + 0.025 + 0.05) ;
        yi_slit_edge_min[i] = std::abs(xy[i] - 0.5 + 0.025 + 0.05 + 0.05) ;

    }    

    int x_i = arma::index_min(xi) ;
    int y_i = arma::index_min(yi) ;

    int y_c_e_min = arma::index_min(yi_center_edge_min) ;
    int y_s_e_max = arma::index_min(yi_slit_edge_max) ;
    int y_s_e_min = arma::index_min(yi_slit_edge_min) ;


    // defining wall thickness   
    int wall_thickness = 0.02/(1./M) ;

    // finding the lenght of the slits/walls
    int len_slits = (y_i + (y_i - y_c_e_min)) - y_c_e_min ;
    
    // creating wall with slits 
    for (int i=0; i< wall_thickness+1; i++)
    {
        V.col(x_i + i).fill(v0) ;

        // center slit is lenght 11
        for (int j = 0; j < len_slits + 1; j++)
        {
            V.col(x_i+i)[y_c_e_min + j] = 0 ;
        }
        // outer slits/wall are length 10
        for (int j = 0; j < len_slits; j++)
        {
            V.col(x_i+i)[y_s_e_min + j] = 0 ;
            V.col(x_i+i)[   y_i + (y_i - y_s_e_max) + 1 + j    ] = 0 ;
        }
    }
}


