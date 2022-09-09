#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <chrono>

// declearing functions
std::vector<double> thomas_algorithm_special_case( std::vector<double> g);

std::vector<double> thomas_algorithm(std::vector<double> a, std::vector<double> b
                                  , std::vector<double> c, std::vector<double> g);

int writefile(std::vector<double> x, std::vector<double> res, std::string name);

int main(){
    //three dimantional vector containing the averaged running time of the general and special algorithems. shape: (2,6,10)
    std::vector< std::vector< std::vector<double> > > Times(2, std::vector<std::vector<double>> (6, std::vector<double>(10)));

    // vector containing number of steps up to 10^6
    std::vector<int> n{10,100,1000,10000,100000,1000000};


    //itterating over the different n
    for (int i = 0; i<(int) n.size(); i++){

        // vectors a, b and c making up the diagonals in addition to x and g
        //b has the same lenght as g. n steps in compleat solution gives n+1 points. 
          //removing the boundaries leaves us with n-1 points. a and c has one less point
        std::vector<double> a (n[i]-2,-1);
        std::vector<double> b (n[i]-1,2);
        std::vector<double> c (n[i]-2,-1);
        
        std::vector<double> x(n[i]+1 );
        std::vector<double> g(n[i]-1);

        // stepsize
        double h = 1./n[i];
        // filling the x axis
        for (int j = 0; j<=n[i]; j++){
            x[j] = h*j;
        }
        
        //known values
        int v_0 = 0;
        int v_n = 0;

        // finding g
        g[0] = 100*std::exp(-10*x[1])*h*h + v_0;
        for (int j = 1; j<n[i]-1; j++){
            g[j] = 100*std::exp(-10*x[j+1])*h*h;
        }

        g[n[i]-2] = 100*std::exp(-10*x[n[i]-1])*h*h + v_n;

        // for each different n rund the general and special algorithems 10 times
        for (int p = 0; p<10; p++){




            // Start measuring time
            auto t1 = std::chrono::high_resolution_clock::now();


            // general algorithem
            thomas_algorithm(a, b, c, g); 


            // Stop measuring time
            auto t2 = std::chrono::high_resolution_clock::now();

            // Calculate the elapsed time
            // We use chrono::duration<double>::count(), which by default returns duration in seconds
            double duration_seconds_general = std::chrono::duration<double>(t2 - t1).count();






            // Start measuring time
            auto t12 = std::chrono::high_resolution_clock::now();


            // special algorithem
            thomas_algorithm_special_case(g); 


            // Stop measuring time
            auto t22 = std::chrono::high_resolution_clock::now();

            // Calculate the elapsed time
            // We use chrono::duration<double>::count(), which by default returns duration in seconds
            double duration_seconds_special = std::chrono::duration<double>(t22 - t12).count();


            // saveing the time used by the two algorithems in Times
            Times[0][i][p] = duration_seconds_general;
            Times[1][i][p] = duration_seconds_special;


        }
    }


    double sum1;
    double sum2;
    std::vector<double> avrg_general(6);
    std::vector<double> avrg_special(6);

    // itterating over different n
    for (int i = 0; i<6 ; i++){
        // itterating over the 10 tests
        for (int j= 0; j<10; j++){
            // summing all the 10 test for the two algorithems
            sum1 = sum1 + Times[0][i][j];
            sum2 = sum2 + Times[1][i][j];

        //findign the average time used by dividing by 10
        avrg_general[i] = sum1/10;
        avrg_special[i] = sum2/10;

        }
    }

    //writing the avraged times to a file
    std::string name = "general_vs_special.txt";
    writefile(avrg_general, avrg_special, name);


    return 0;
}


//thomas algorithem
std::vector<double> thomas_algorithm(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> g){
    
    // creating a copys of the vectors whoes values will be changed
    std::vector<double> G(g.size());
    std::vector<double> B(b.size());
    std::vector<double> V(g.size());

    int n = g.size();
    // inital values of B and G are not changed
    B.at(0) = b.at(0);
    G.at(0) =  g.at(0);
    
    // forward elimintion
    for (int i=1; i<n; i++){

        double a_B = a.at(i-1)/B.at(i-1);//a[i-1] because the first value in row-2 is at index 0

        G.at(i) = g.at(i) - a_B*G.at(i-1);
        B.at(i) = b.at(i) - a_B*c.at(i-1);

    }
    
    // bacward elimination
    // first find the solution at last index
    V.at(n-1) = G.at(n-1)/B.at(n-1);
    for (int i = 2; i<=n ; i++){

        V.at(n-i) = (G.at(n-i) - V.at(n-i+1)*c.at(n-i))/B.at(n-i);
    }

    return V;

}







// thomas algorithm for special case
std::vector<double> thomas_algorithm_special_case( std::vector<double> g){
    
    // creating a copys of the vectors whoes values will be changed
    std::vector<double> G(g.size());
    std::vector<double> B(g.size());
    std::vector<double> V(g.size());

    int n = g.size();
    // inital values of B and G are not changed
    B.at(0) = 2.;
    G.at(0) =  g.at(0);
    
    // forward elimintion
    for (int i=1; i<n; i++){

        G.at(i) = g.at(i) + G.at(i-1)/B.at(i-1);
        B.at(i) = (i+2)/((double) (i+1) );       

    }

    // bacward elimination
    // first find the solution at last index
    V.at(n-1) = G.at(n-1)/B.at(n-1);
    for (int i = 2; i<=n ; i++){

        V.at(n-i) = (G.at(n-i) + V.at(n-i+1))/B.at(n-i);
    }

    return V;

}







// function used to write a file with two colums
int writefile(std::vector<double> x, std::vector<double> y, std::string name) {

    std::ofstream ofile;
    ofile.open(name);

    // parameters used for writing to file neatly
    int width = 12;
    int prec = 6;

    //writing one x and y value per line
    for (int i = 0; i < x.size(); i++){
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x[i] << "    " 
        << std::setw(width) << std::setprecision(prec) << std::scientific  << y[i] << "\n" ;

    }

    ofile.close();

    return 0;

}