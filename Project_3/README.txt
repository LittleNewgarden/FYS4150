


We split the project into two main files. One for problem 8 and one for problem 9. main_t8 runs very fast, but main_t9 takes ~2 hours on our machine. 
The python file does the visualisation and calculation conected with both tasks. 


Problem 8:
  g++-11 main_task8.cpp Particle.cpp PenningTrap.cpp -o main_t8.exe -larmadillo -O2
 
  ./main_t8.exe
  
 
  
 Problem 9: 
  g++-11 main_task9.cpp Particle.cpp PenningTrap.cpp -o main_t9.exe -larmadillo -O2
  
  ./main_t9.exe
   

Visualisation and calculation of specific values:
  py penningtrap.py
   
