

We use Windows and thus need to have "-larmadillo" at the end of build line. For example "g++ Problem5.cpp -o p5.exe -larmadillo"



Problem 2:
    g++ Problem2.cpp -o p2.exe
    ./p2.exe

running p2.exe will print matrix A, armadillo solutions for eigenvalues and eigenvectors and analytical solutions for eigenvalues and eigenvectors

  



Problem 3:
    g++ Problem3.cpp -o p3.exe
    ./p3.exe
    
running p3.exe will print the matrix A and the largest off diagonal value.





Problem 4:
    g++ Problem4.cpp -o p4.exe
    ./p4.exe

running p4.exe will print the egienvalues and eigenvectors found using both the jacobirotation algorithm and the analytical solution.
We also print the number of iterations ain the jacobi rotation algorithm and whether or not the jacobi rotation algorithm converged.
 
 
 
 

Problem 5:
    g++ Problem5.cpp -o p5.exe
    ./p5.exe
    py Problem5.py

Running problem5.py will give a plot.





Problem 6:
    g++ Problem6.cpp -o p6.exe
    ./p6.exe
    py Problem6.py
    
running Priblem6.py will give 6 plots
