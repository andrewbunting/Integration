#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/*
 
 This produces the coordinates to sample a circle of radius 1 with each point corresponding to an equal area of the circle (and therefore an equal projected area if you're looking at a sphere)
 
 N is the number of rings of points
 
 beta_0 is the number of points on the innermost ring (we require beta_0 > 1, in order to make it make sense at the centre, and if beta_0 \approx Pi, then the areas will be approximately square (so betq_0 = 3 is a good choice))
 
 alpha counts which ring you are on (it runs from 0 (at the centre) to N-1 (at the outer edge))
 
 j counts how far around in phi you are on a ring (it runs from 0 to (beta_0*(2*alpha + 1) - 1))
 
 */


using namespace std;

int main()
{
    int alpha, j, N, beta_0, beta_alpha;
    
    double x, y, theta, phi, rho, alpha_doub, N_doub, beta_0_doub, j_doub,pi;
    
    pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148;
    
    // number of points in the innermost ring
    beta_0 = 3;
    
    beta_0_doub = (double) (beta_0);
    
    // number of rings
    N = 100;
    
    N_doub = (double) (N);
    
    // file with all of the sample locations
    ofstream outfile;
    
    outfile.open("Output/disk_points.dat", ios::out);
    
    outfile.precision(10);
    
    cout << "N = \t" << N << "\n";
    cout << "beta_0 = \t" << beta_0 << "\n\n";
    
    outfile << "x" << "\t\t\t\t\t" << "y" << "\t\t\t\t\t" << "theta" << "\t\t\t\t\t" << "phi" << "\n";
    
    for (alpha=0; alpha < N; alpha = alpha+1) {
        
        alpha_doub = (double) (alpha);
        
        beta_alpha = beta_0*((2*alpha) + 1);
        
        theta = asin((alpha_doub + 0.5)/N_doub);
        
        for (j=0; j<beta_alpha; j=j+1) {
            
            j_doub = (double) (j);
            
            rho = ((alpha_doub + 0.5)/N_doub);
            
            phi = (pi/beta_0_doub) * ( ((2*j_doub) + 1.0)/((2*alpha_doub) + 1.0) );
            
            x = rho * cos(phi);
            
            y = rho * sin(phi);
            
            outfile << x << "\t\t\t\t\t" << y << "\t\t\t\t\t" << theta << "\t\t\t\t\t" << phi << "\n";
            
        }
        
    }
    
    
}



