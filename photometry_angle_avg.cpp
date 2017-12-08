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
    int alpha, j, N, beta_alpha, k, kmax;
    
    double x, y, theta, phi, rho, alpha_doub, N_doub, j_doub, pi, k_doub, kmax_doub, beta_alpha_doub;
    
    // Here I define the variables to be used in the photometry side of things
    double L0, Lprime, f0, f_r, f_i, delta_f, df_dr, xi_r, xi_i, c_x, u_x, R, r, omegat, theta_0, phi_0, C_1, phase, nprime_dot_n_obs, h_x, dA, dA0, dL, r_dot_n_obs, dL_f, dL_dot;
    
    // to do with the change in area of a surface element
    C_1 = 1.0;
    
    // to do with normalising the limb-darkening, but will cancel out in this case because we'll look at the ratio of perturbed luminosity to equilibrium luminosity
    c_x = 1.0;
    
    
    // Thisvalue determines the kind of limb-darkening that you're using (0.6 for Eddington limb-darkening)
    u_x = 0.6;

    
    
    
    // Here the necessary physical values are being input, although this will need to be changed in order to take these from a file automatically
    
    f0 = 64274443019.8;
    
    f_r = 283124.324612705;
    
    f_i = 1327054.54713986;
    
    df_dr = -1.80824868386493;

    R = 71090420513.5;
    
    xi_r = -1152.30948016823;
    
    xi_i = -1218.44658087104;
    
  
    
    
    
    
    
    pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148;
    
    
    
    // Here you set the angle from which the observer is viewing the star, according to the coordinates in which the perturber is at the equator
    
    
    theta_0 = 0.5*pi;
    
    phi_0 = 0.0;
    
    
    // this sets the number of time-steps taken to cover omegat going from 0 to (just under) 2*pi
    kmax = 100;
    
    kmax_doub = (double) (kmax);
    
    
    // number of rings
    N = 4000;
    
    N_doub = (double) (N);
    
    // file with all of the sample locations
    ofstream outfile;
    
    outfile.open("Output/luminosities_angle_avg.dat", ios::out);
    
    outfile.precision(15);
    
    cout << "N = \t" << N << "\n";
    
    cout << "kmax = \t" << kmax << "\n";
    
    outfile <<  "omega*t" << "\t\t\t\t\t" << "L0, analytically" << "\t\t\t\t\t" << "L0, numerically" << "\t\t\t\t\t" << "Lprime" << "\t\t\t\t\t" << "( Lprime - L0 ) / L0" << "\t\t\t\t\t" << "dL_f/L0" << "\t\t\t\t\t" << "dL_dot/L0"  << "\n";
    
    for (k=0; k < kmax; k = k+1) {
        
       
        k_doub = (double) (k);
        
        if ((100*k/kmax)%10 == 0) {
            cout << 100*k/kmax << "% complete\n";
        }
        
        omegat = 2.0*pi*k_doub/kmax_doub;
        
        // This resets the L0 and Lprime values so they can be re-summed at this new time
        
        L0 = 0.0;
        
        Lprime = 0.0;
        
        dL_f = 0.0;
        
        dL_dot = 0.0;
        
        
        
        for (alpha=0; alpha < N; alpha = alpha+1) {
            
            alpha_doub = (double) (alpha);
            
            theta = pi*(alpha_doub + 0.5)/N_doub;
            
            beta_alpha = 2*N*sin(theta);
            
            beta_alpha_doub = (double) (beta_alpha);
            
            for (j=0; j<beta_alpha; j=j+1) {
                
                j_doub = (double) (j);
                
                phi = (2.0*pi/beta_alpha_doub) * ( j_doub + 0.5 );
                
                
            // Here is where the evaluation of the bits in the sum actually happens
                
                
                // preliminaries
                
                phase = 2.0*( omegat - phi );
                
                r = R + 3.0*sin(theta + theta_0)*sin(theta + theta_0)*( xi_r*cos(phase) - xi_i*sin(phase) );
                
                r_dot_n_obs = (sin(theta)*cos(phi)*sin(theta_0)*cos(phi_0)) + (sin(theta)*sin(phi)*sin(theta_0)*sin(phi_0)) + (cos(theta)*cos(theta_0));
                
                
                // As we integrate over the whole surface, but only want to include bits that we can see, we only include the surface elements with a positive dot product with the observer's vector
                
                if (r_dot_n_obs > 0.0) {
                    
                    // flux
                    
                    delta_f = 3.0*sin(theta)*sin(theta)*( (f_r + df_dr*xi_r)*cos(phase) - (f_i + df_dr*xi_i)*sin(phase) );
                    
                    
                    // normal
                    
                    nprime_dot_n_obs = ( 6.0*(sqrt(1.0 - (r_dot_n_obs*r_dot_n_obs)))*sin(theta)*cos(theta)*( xi_r*cos(phase) - xi_i*sin(phase) ) )/R;
                    
                    
                    // limb-darkening (in filter x)
                    
                    h_x = c_x * ( 1.0 - u_x*( 1.0 - r_dot_n_obs ) );
                    
                    
                    // dA
                    
                    dA = (R*R) * sin(theta) * ( (2.0*pi*pi) / (N_doub * beta_alpha_doub) );
                    
                    
                    // Here's the actual integral bit
                    
                    // for the perturbed star as a whole
                    
                    dL = ( f0*r_dot_n_obs + delta_f*r_dot_n_obs + f0*nprime_dot_n_obs) * h_x * dA;
                    
                    Lprime = Lprime + dL;
                    
                    // for the equilibrium star
                    
                    dL = f0 * r_dot_n_obs * h_x * dA;
                    
                    L0 = L0 + dL;
                    
                    // for the individual perturbing bits
                    
                    // change in flux
                    dL = delta_f*r_dot_n_obs * h_x * dA;
                    
                    dL_f = dL_f + dL;
                    
                    // change in normal
                    dL = f0*nprime_dot_n_obs * h_x * dA;
                    
                    dL_dot = dL_dot + dL;
                    
                    
                    
                }
                
                
                
                
                
            }
            
        }
        
        
        // This outputs the luminosities
        
        outfile << omegat << "\t\t\t\t\t" << pi*R*R*f0*c_x*(1.0 - (u_x/3.0)) << "\t\t\t\t\t" << L0 << "\t\t\t\t\t" << Lprime << "\t\t\t\t\t" << (Lprime - L0) / L0 << "\t\t\t\t\t" << dL_f/L0 << "\t\t\t\t\t" << dL_dot/L0<< "\n";
        
        
        
    }
    

    
    
}



