#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <iomanip>
#include <fstream>




std::vector<double> thomas_algorithm(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> g);

int writefile(std::vector<double> x, std::vector<double> res, std::string name);

int main(){

    std::vector<int> n{10,100,1000,100000};
    for (int i = 0; i<4; i++){

        std::vector<double> a (n[i]-2,-1);
        std::vector<double> b (n[i]-1,2);
        std::vector<double> c (n[i]-2,-1);
        
        std::vector<double> x(n[i]+1 );
        std::vector<double> g(n[i]-1);

        double h = 1./n[i];
   
        for (int j = 0; j<=n[i]; j++){
             x[j] = h*j;
        }

        int v_0 = 0;
        int v_n = 0;

        g[0] = 100*std::exp(-10*x[1])*h*h + v_0;
        for (int j = 1; j<n[i]-1; j++){
            g[j] = 100*std::exp(-10*x[j+1])*h*h;
        }

        g[n[i]-2] = 100*std::exp(-10*x[n[i]-1])*h*h + v_n;

        std::vector <double> V = thomas_algorithm(a, b, c, g); 

        std::string file_name = "P7_" + std::to_string(n[i]) + ".txt";

        V.push_back(0.);

        V.insert(V.begin(), 0);

        writefile(x, V, file_name);
    

    }







    return 0;
}



std::vector<double> thomas_algorithm(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> g){
    

    std::vector<double> G(g.size());
    std::vector<double> B(b.size());
    std::vector<double> V(g.size());

    int n = g.size();
    B.at(0) = b.at(0);
    for (int i=1; i<n; i++){

        double a_B = a.at(i-1)/B.at(i-1);//a[i-1] because the first value in row-2 is at index 0

        G.at(i) = g.at(i) - a_B*g.at(i-1);
        B.at(i) = b.at(i) - a_B*c.at(i-1);

    }
    
    
    V.at(n-1) = G.at(n-1)/B.at(n-1);
    for (int i = 2; i<=n ; i++){

        V.at(n-i) = (G.at(n-i) - G.at(n-i+1)*c.at(n-i))/B.at(n-i);
    }

    return V;

}

int writefile(std::vector<double> x, std::vector<double> res, std::string name) {
    std::string filename = name;

    std::ofstream ofile;
    ofile.open(filename);

    int width = 12;
    int prec = 6;

    for (int i = 0; i < (int) x.size(); i++){
        ofile << x[i] << "    " <<std::setw(width) << std::setprecision(prec) << std::scientific  << res[i] <<"\n" ;
    }
    ofile.close();

    return 0;

}