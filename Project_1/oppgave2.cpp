#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <iomanip>
#include <fstream>



// declearing u(x) function and a function used to write a file
std::vector<double> u(std::vector<double> x);

int writefile(std::vector<double> x, std::vector<double> res);



int main() {
    //creating the vector x of size 101
    std::vector<double> x(101);
    //creating stepsize h
    double h = 1/100.; 
    
    //filling the x vector with evenly spaced values betweeen 0 and 1
    for (int i = 0; i <=100; i++){
        x[i] = h*i;
    }

    // finding u(x)
    std::vector<double> res = u(x);
    
    //writing file
    writefile(x, res);
    

    return 0;
}



// function finding u(x)
std::vector<double> u(std::vector<double> x){

    // creating vector res and filling it with values for all x
    std::vector<double> res(x.size());
    for (int i = 0; i < x.size(); i++){
        res[i] = 1 - (1 - std::exp(-10))*x[i] -std::exp(-10*x[i]);
    
    }

    return res; 
}


// function used to write a file with two colums
int writefile(std::vector<double> x, std::vector<double> y) {
    //name of file
    std::string filename = "u_x_output.txt";

    std::ofstream ofile;
    ofile.open(filename);

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