#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <iomanip>
#include <fstream>


std::vector<double> u(std::vector<int> x);

int writefile(std::vector<int> x, std::vector<double> res);

int main() {
    std::vector<int> vec {1,2,3,4,5,6,7,8,9,10};

    std::vector<int> vec2(101);

    for (int i = 0; i <=100; i++){
        vec2[i] = i;
    }

    std::vector<double> res = u(vec2);
    
    writefile(vec2, res);

    return 0;
}




std::vector<double> u(std::vector<int> x){

    int width = 12;
    int prec = 6;
    std::vector<double> res(x.size());
    for (int i = 0; i < x.size(); i++){
        res[i] = 1 - (1 - std::exp(-10))*x[i] -std::exp(-10*x[i]);
        std::cout  << x[i] << "    " <<std::setw(width) << std::setprecision(prec) << std::scientific  << res[i] <<"\n" ;
        //<< std::setw(width) << std::setprecision(prec) << std::scientific<< res[i] <<"\n";
    }

    
    return res; 
}



int writefile(std::vector<int> x, std::vector<double> res) {
    std::string filename = "u_x_output.txt";

    std::ofstream ofile;
    ofile.open(filename);

    int width = 12;
    int prec = 6;

    for (int i = 0; i < x.size(); i++){
        ofile << x[i] << "    " <<std::setw(width) << std::setprecision(prec) << std::scientific  << res[i] <<"\n" ;
        //ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x[i] << std::setw(width) << std::setprecision(prec) << std::scientific  << res[i] <<"\n" ;

    }
    ofile.close();

    return 0;
}