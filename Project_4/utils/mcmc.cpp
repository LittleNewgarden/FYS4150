#include <iostream>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;


// Initialising the refence matrix using the state matrix
vector<vector<double *>> ref_mat(mat& S)
{
    int height = S.n_rows;
    int width = S.n_cols;

    vector<vector<double *>> Sref(height+2,vector<double*>(width+2, 0)) ;

    for (int i=0; i< height; i++)
    {
        for (int j = 0; j<width;j++)
        {
            Sref[i+1][j+1] = &S.row(i)[j];
        }
    }

    for (int i = 0; i < height; i++)
    {
        Sref[i+1][0] = &S.row(i)[width-1];
        Sref[i+1][width+1] = &S.row(i)[0];
    }

    for (int j=0; j<width;j++)
    {
        Sref[0][j+1] = &S.row(height-1)[j];
        Sref[height+1][j+1] = &S.row(0)[j];
    }

    return Sref;
}


// One MC cycle finding the new energy, new magnetization and the new state S
void MC(int NN, double T, mat& S,vector<vector<double *>>& Sref, double& E0, int& M0)
{   
    int height = S.n_rows;
    int width = S.n_cols;

    for (int q=0; q<NN;q++)
    {
        // CHOOSE RANDOM INDEX TO FLIP

        int indx_i =  randi(distr_param(0,height-1)) ;
        int indx_j =  randi(distr_param(0,width-1)) ;


        int ref_i = indx_i + 1 ;  // indecis of reference matrix
        int ref_j = indx_j + 1 ;


        // CALCULATE P/P if spin is fliped

        double dE = 2*(*Sref[ref_i][ref_j])*( *Sref[ref_i-1][ref_j] + *Sref[ref_i+1][ref_j] + *Sref[ref_i][ref_j-1] + *Sref[ref_i][ref_j+1] );

        double beta = 1/(T);

        double P_P = exp(-beta*dE);


        // ACCEPT REGECT

        double A = 1;

        if (P_P < A)
        {
            A = P_P;
        }

        double r = randu();
        int dM;

        if (r <= A)
        {
            S.row(indx_i)[indx_j] = S.row(indx_i)[indx_j]*(-1); 
            dM = S.row(indx_i)[indx_j]*(2);
        }
        else
        {
            dE = 0;
            dM = 0;
        }

        E0 = E0 + dE;
        M0 = M0 + dM;
    }
}


// outer loop of the MCMC method. saves the values after each MC cycle
mat MCMC(int N, int NN, double T, mat& S, vector<vector<double *>>& Sref, double& E0, int& M0)
{
    vec E = vec(N);
    vec M = vec(N);
    E[0] = E0;
    M[0] = M0;
    for (int i = 0; i< N; i++)
    {
        MC(NN, T, S, Sref, E0, M0);
        E[i+1] = E0;
        M[i+1] = M0;

    }
    mat EM = mat(N,2);
    EM.col(0) = E;
    EM.col(1) = M;
    
    return EM;
}



// function used to initiallise the state matirx and return the initial magnetization based on the percentage of spin down
int initial_S_M0(mat& S, double percent)
{
    int M_0=0 ; 
    for (int i=0; i < (int) S.size(); i++) 
    {
        if ( S[i]>=percent)
        {
            S[i] = 1 ;
            M_0 += 1 ;
        }
        else
        {
            S[i] = -1 ;
            M_0 += -1 ;
        }
    }
    return M_0 ;
}


// function used to initiallise the reference matirx and return the initial energy
double initial_Sref_E0(mat& S, vector<vector<double *>>& Sref)
{
    Sref = ref_mat(S) ; // initialising the reference matirx

    double E_0 = 0 ;
    
    int height = S.n_rows ;
    int width = S.n_cols ;
    for (int i=1; i<height+1; i++ )  // cfinding initial energy
    {
        for (int j=1; j<width+1; j++)
        {
            E_0  += *Sref[i][j]*( *Sref[i+1][j] + *Sref[i][j+1] ) ;
        }
    }
    E_0 = -E_0 ;
    return E_0 ;
}

