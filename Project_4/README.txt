In order to simulate a ferromagnet, we have implemented the Ising model using the Markov Chain Monte Carlo method together 
with a modified version of the Metropolis algorithm

As the project is divided into problems, with problem  4 - 9 being based on numerical implementation, we have
made different code python files for solving the particular problems, but we only have one main file called p4_all_in_one.cpp
which gives us all the files used in the python scripts if given the coroect command line arguments.

In the folder figs/ we have given the figures and in the folder utils/ we have one hpp file and one cpp file. The cpp file contains all the functions we 
use in the main file and the hpp file contains their declerations.




In order to get the parallelized program we run the following comand in the command line:  
    g++ p4_all_in_one.cpp utils/mcmc.cpp -o p4_aio_O2.exe -larmadillo -Wall -O2 -fopenmp



In order to get the non-parallelized program we run the following comand in the command line:  
    g++ p4_all_in_one.cpp utils/mcmc.cpp -o p4_aio_O2.exe -larmadillo -Wall -O2
    
    
    
 when running the program ./p4_aio_O2.exe, there will come text up in the command line asking for imputs.
 
 To solve problem 4 run the following comands in the command line:
    ./p4_aio_O2.exe
      Latice side lenght :2
      percentage of spin -1 :.5
      Start temperature :1
      End temoerature :1
      Number of temperatures :1
      number of MC sycles :1000000
      Write all E and M values to file? :1
      Write average values of E and M to file? :0
      Name of file :2x2_n=10^6
    
    py p4.py
    
    
 To solve problem 5 run the following comands in the command line:
    ./p4_aio_O2.exe
      Latice side lenght :20
      percentage of spin -1 :.5
      Start temperature :1
      End temoerature :2.4
      Number of temperatures :2
      number of MC sycles :1000000
      Write all E and M values to file? :1
      Write average values of E and M to file? :0
      Name of file :20x20_per=50_n=10^6
      
    ./p4_aio_O2.exe
      Latice side lenght :20
      percentage of spin -1 :1
      Start temperature :1
      End temoerature :2.4
      Number of temperatures :2
      number of MC sycles :1000000
      Write all E and M values to file? :1
      Write average values of E and M to file? :0
      Name of file :20x20_per=100_n=10^6
      
     py p5.py
      
      
Problem 6 needs no new text files:
    py p6.py
    
    
Problem 7 needs no new text files:
    py p7.py
    
    
 To solve problem 8-9 run the following comands in the command line:
    ./p4_aio_O2.exe
      Latice side lenght :40
      percentage of spin -1 :.5
      Start temperature :2.1
      End temoerature :2.4
      Number of temperatures :100
      number of MC sycles :1000000
      Write all E and M values to file? :0
      Write average values of E and M to file? :1
      Name of file :40x40_n=10^6_Tn=100
      
      
    ./p4_aio_O2.exe
      Latice side lenght :60
      percentage of spin -1 :.5
      Start temperature :2.1
      End temoerature :2.4
      Number of temperatures :100
      number of MC sycles :1000000
      Write all E and M values to file? :0
      Write average values of E and M to file? :1
      Name of file :60x60_n=10^6_Tn=100
      
    ./p4_aio_O2.exe
      Latice side lenght :80
      percentage of spin -1 :.5
      Start temperature :2.1
      End temoerature :2.4
      Number of temperatures :100
      number of MC sycles :1000000
      Write all E and M values to file? :0
      Write average values of E and M to file? :1
      Name of file :80x80_n=10^6_Tn=100
      
    ./p4_aio_O2.exe
      Latice side lenght :100
      percentage of spin -1 :.5
      Start temperature :2.1
      End temoerature :2.4
      Number of temperatures :100
      number of MC sycles :1000000
      Write all E and M values to file? :0
      Write average values of E and M to file? :1
      Name of file :100x100_n=10^6_Tn=100
      
    py p8-9.py
