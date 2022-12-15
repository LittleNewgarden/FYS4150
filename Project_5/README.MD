In this project we have solved the unitless time-dependent Schrödinger equation for a non-relativistic particle with a mass with a time independant potential. 
The paricle is initially defined with an 2D gaussian wave function u^0 at time step n=0, we solved for u^(n+1) using Au^(n+1) = Bu^n


Our proganised ito two folders, include, src. In each of these there is a single file named all_functions.hpp and all_functions.cpp respecfully.
In all_functions.cpp we have defined all the function called appon in the main file : p5.cpp. The declearations for the functions are in all_functions.hpp.

We also have the folder "figures" in this folder we have all the plots we used in our paper.


We where for some reason unable to make use of the python armadillo library. Because of this we have used raw_ascii when saving the data. 


To reproduce our results you need to first compile the program using:

    g++ p5.cpp src/all_functions.cpp -o p5.exe -larmadillo -Wall -O2


The different comand line argument that the program take in are: step size in x and y "h", time step "dt", end time "T", 
center of 2D gaussian in "xc" and "yc", the sigma of the gaussian "sigma_x" and "sigma_y", the momentum "px" and "py", 
the potential value in the wall "v0", the name of the output file "name", the type of potential wanted "slit"
and finally if we are going to save all u^n "save_u".
 
  
Then to solve the tasks you need to first run the program for the different needed senarios, then you need to run the python script to get the results.

To produce the data needed for the different tasks run the following comands:

  for task 7 :

    ./p5.exe 0.005 0.000025 0.008 0.25 0.5 0.05 0.05 200 0 0 task7_1 no_slit false 
    ./p5.exe 0.005 0.000025 0.008 0.25 0.5 0.05 0.1 200 0 10000000000 task7_2 double false

  for task 8 :
    
    ./p5.exe 0.005 0.000025 0.002 0.25 0.5 0.05 0.2 200 0 10000000000 task8 double true

  for task 9 :
    
    ./p5.exe 0.005 0.000025 0.002 0.25 0.5 0.05 0.2 200 0 10000000000 task9_single single false
    ./p5.exe 0.005 0.000025 0.002 0.25 0.5 0.05 0.2 200 0 10000000000 task9_triple triple false



The double slit P file from task 8 is also used in task 9


Then run the python script:
  
    py p5.py
  

