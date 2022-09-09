#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <iomanip>
#include <fstream>



// declearing functions
std::vector<double> thomas_algorithm(std::vector<double> a, std::vector<double> b
                                  , std::vector<double> c, std::vector<double> g);

int writefile(std::vector<double> x, std::vector<double> res, std::string name);

int main(){
    // vector of n's that will be needed for task 7 and 8
    // n is her the number of steps in the compleat solution
    std::vector<int> n{10,100,1000,10000,100000,1000000,10000000};
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

        //finding the solution 
        std::vector <double> V = thomas_algorithm(a, b, c, g); 

        //name of each file
        std::string file_name = "P7_" + std::to_string(n[i]) + ".txt";

        // adding the first and last value
        V.push_back(0.);
        V.insert(V.begin(), 0.);

        //writing file
        writefile(x, V, file_name);
    

    }







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