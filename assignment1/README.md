# README

### Instructions for running the script
1. Clone the repository
2. Run `make` in the command line while in the `assignment1` folder
   - This will run the makefile which runs everything we need
   - I have already included the output files in the project_outputs folder in case there are any issues with running.
     - Feel free to remove the files from the outputs folder to get newly generated files
    
### File Breakdown
root/  
│── include/ - These are the header files 
│   ├── atom.hpp  
│   ├── cluster.hpp  
│   ├── derivative_approximation.hpp  
│   ├── lennard_jones.hpp  
│── src/  
│   ├── atom.cpp                        - Atom class  
│   ├── cluster.cpp                     - System of Atoms - includes steepest descent methods  
│   ├── derivative_approximation.cpp    - Forward and Central Difference  
│   ├── lennard_jones.cpp               - Lennard Jones Calculations  
│   ├── main.cpp                        - Main  
│   │── plot_generation.ipynb           - Generates plot  

### Important things to know
1. Output files are generated. In case this does not work with your system, I left commented out code for each portion of the questions to print to the terminal.
2. I added the `plot_generation.ipynb` file to generate the plots for round-off and truncation error. Please know that I input the data that my script generates for this. I would've loved to link it together, but that would be a lot of work that I couldn't prioritize.

### Plots discussion
Based on the step size changes, the truncation error starts to increase exponentially, however when step size is smaller, then round-off error starts to have the higher error percentage, although very small. Based on this, the optimal step size would be one where both the on the smaller side 0.00001 or lower. However, I believe the round-off error should be increasing as we have smaller step sizes due to the amount of precision we need for our numbers. To my understanding, there must be something that is affecting my script's ability to generate more precise numbers. Usually the smaller the step size, the greater the round-off error, so I would've expected to see the round-off error getting higher and higher. This way the intersection of the round-off and truncation error would most likely be our the optimal step size.

### Known Issues - This needs review and if I have the time this weekend, then I will come back to fix this.
1. As I've been transforming my code, when I incorporated armadillo vectors and matrices, this seemed to reduce my precision. I tried fixing this by setting the precision when outputting values. I think the issue might be with precision when doing calculations. The precision might already be truncated/reduced by the time the calculations are made. IMPORTANTLY, for question 1.2 I changed my outputs to output the force vector in the atom class because we could actually see better differences with forward and central difference, HOWEVER if it wasn't for the low precision for my matrices, I would've used those. I did output coordinate matrices for the rest of the assignment though.
2. Originally I was going to develop a steepest descent .cpp and .hpp file, however the intertwining with cluster.cpp felt necessary for the steepest descent portion of this homework. Unfortunately, I did not have enough time to be able to fix this and create a totally well-organized repo. I started off very strong, and then many problems came about with the last question.
3. The standard_SD and SD_with_line_search portions of this are not working. I implemented gradient descent, bracketing, and golden section search, and it does not work. I don't know why, I've stared at it for a long time and it feels like too many moving parts. I included code directly from the numerical recipes book and I've tried other methods, which is why there are duplicate methods in there. I want to go back and revamp the code.
4. Also, I apologize for the gradescope submission. I thought it was just the repo that would be used to review my code.
