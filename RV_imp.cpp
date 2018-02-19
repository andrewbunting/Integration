#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>


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
    int alpha, j, N, beta_alpha, k, kmax, bins, i;
    
    double x, y, theta, phi, rho, alpha_doub, N_doub, j_doub, pi, k_doub, kmax_doub, beta_alpha_doub;
    
    // Here I define the variables to be used in the photometry side of things
    double L0, f0, f_r, f_i, xi_r, xi_i, c_x, u_x, R, r, omegat, theta_0, phi_0, C_1, phase, h_x, dA, dA0, dL, r_dot_n_obs;
    
    // Here I define the variables needed on top of the photometry variables in order to do the RV side of things
    double omega, m, p_prime_r, p_prime_i, rho_0, f_grav, theta_dot_n_obs, phi_dot_n_obs, v_RV_max, v_RV_r, v_RV_theta, v_RV_phi, v_RV_tot, bins_doub, delta_v, i_doub, v_bin;

    
    // This deals with how the ata is output - a different file for each point in time
    
    string prefix = "Output/RV_";
    string suffix = ".dat";
    string filename;
    string k_str;


    
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
    
    //df_dr = -1.80824868386493;

    R = 71090420513.5;
    
    xi_r = -1152.30948016823;
    
    xi_i = -1218.44658087104;
    
    
    
    m = 2.0;
    
    omega = 1.719741e-5;
    
    p_prime_r = -4.44911111;
    
    p_prime_i = 54808.82553*(-8.727223541e-5);
    
    rho_0 = 1.48995335133e-07;
    
    f_grav = -7.0482251937123e-14;
  
    
    
    v_RV_max = 6.0*omega*sqrt( (xi_r*xi_r) + (xi_i*xi_i) + (8.0/((R*m*m*omega*omega)*(R*m*m*omega*omega)))*(  ( (p_prime_r/rho_0) + f_grav*R*R )*( (p_prime_r/rho_0) + f_grav*R*R )    +    ( p_prime_i/rho_0 )*( p_prime_i/rho_0 )  ) );
    
    cout << "v_RV_max =\t" << v_RV_max << "\n";
    
    
    
    
    
    pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148;
    
    
    
    // Here you set the angle from which the observer is viewing the star, according to the coordinates in which the perturber is at the equator
    
    
    theta_0 = 0.5*pi;
    
    phi_0 = 0.0;
    
    
    // This sets the number of bins for RV which will be used
    
    bins = 500;
    
    bins_doub = (double) (bins);
    
    
    
    // Here we define the array which will be used as the bins
    double B[bins];
    
    for (i=0; i<bins; i=i+1) {
        
        B[i] = 0.0;
        
    }
    
    
    
    // this sets the number of time-steps taken to cover omegat going from 0 to (just under) 2*pi
    kmax = 30;
    
    kmax_doub = (double) (kmax);
    
    
    // number of rings
    N = 1000;
    
    N_doub = (double) (N);
    
/*    // file with all of the sample locations
    ofstream outfile;
    
    outfile.open("Output/RV.dat", ios::out);
    
    outfile.precision(15);
    
    outfile <<  "omega*t" << "\t\t\t\t\t" << "L0, analytically" << "\t\t\t\t\t" << "L0, numerically" << "\t\t\t\t\t" << "Lprime" << "\t\t\t\t\t" << "( Lprime - L0 ) / L0" << "\t\t\t\t\t" << "dL_f/L0" << "\t\t\t\t\t" << "dL_dot/L0"  << "\n";
 */
    
    
    cout << "N = \t" << N << "\n";
    
    cout << "kmax = \t" << kmax << "\n";
    
    
    
    
    
    for (k=0; k < kmax; k = k+1) {
        
       
        k_doub = (double) (k);
        
        if ((100*k/kmax)%10 == 0) {
            cout << 100*k/kmax << "% complete\n";
        }
        
        omegat = 2.0*pi*k_doub/kmax_doub;
        
        // This resets the L0 and Lprime values so they can be re-summed at this new time
        
        L0 = 0.0;
        
        // This resets the bin values
        for (i=0; i<bins; i=i+1) {
            
            B[i] = 0.0;
            
        }
        

        
        
        
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
                
                theta_dot_n_obs = cos(theta)*sin(theta_0)*( cos(phi)*cos(phi_0) + sin(phi)*sin(phi_0) )  -  sin(theta)*cos(theta_0);
                
                phi_dot_n_obs = sin(theta_0)*( cos(phi)*sin(phi_0) - sin(phi)*cos(phi_0) );
                
                // As we integrate over the whole surface, but only want to include bits that we can see, we only include the surface elements with a positive dot product with the observer's vector
                
                if (r_dot_n_obs > 0.0) {
                    
                    
                    // limb-darkening (in filter x)
                    
                    h_x = c_x * ( 1.0 - u_x*( 1.0 - r_dot_n_obs ) );
                    
                    
                    // dA
                    
                    dA = (R*R) * sin(theta) * ( (2.0*pi*pi) / (N_doub * beta_alpha_doub) );
                    
                    
                    // Here's the actual integral bit
                    
                    // Here are the individual radial velocity components being computed
                    
                    v_RV_r = -6.0*omega*sin(theta)*sin(theta) * r_dot_n_obs * ( (xi_r*sin(phase)) + (xi_i*cos(phase)) );
                    
                    v_RV_theta = ( (-12.0*sin(theta)*cos(theta))/(R*m*m*omega) ) * theta_dot_n_obs * (  ( (p_prime_r/rho_0) + f_grav*R*R )*sin(phase)  +  ( (p_prime_i/rho_0)*cos(phase) )  );
                    
                    v_RV_phi = ( (12.0*sin(theta))/(R*m*m*omega) ) * phi_dot_n_obs * (  ( (p_prime_r/rho_0) + f_grav*R*R )*cos(phase)  -  ( (p_prime_i/rho_0)*sin(phase) )  );
                    
                    
                    v_RV_tot = v_RV_r + v_RV_theta + v_RV_phi;
                    
                    
                    // This is where you got up to!!
                    dL = (h_x/(2.0*pi)) * f0 * r_dot_n_obs * r_dot_n_obs * dA;
                    
                    
                    for (i=0; i<bins; i=i+1) {
                        
                        delta_v = 2.0*v_RV_max/bins_doub;
                        
                        i_doub = (double) (i);
                        
                        if ( v_RV_tot > -v_RV_max + delta_v*i_doub && v_RV_tot < -v_RV_max + delta_v*(i_doub + 1.0)) {
                            
                            B[i] = B[i] + dL;
                            
                        }
                        
                    }
                    
                    
                    
                    
                    
                    
                }
                
                
                
                
                
            }
            
        }
        
        
        // This sets the name of the file which this time's data will be written to
        
        k_str = to_string(k);
        filename = prefix + k_str + suffix;
        
       
        // file with all of the sample locations
        ofstream outfile;
        
        outfile.open(filename, ios::out);
        
        outfile.precision(15);
        
        outfile <<  "k" << "\t\t\t\t\t" << "omega * t" << "\t\t\t\t\t" << "v_RV" << "\t\t\t\t\t" << "dL" << "\n";
        
        
        for (i=0; i<bins; i=i+1) {
            
            i_doub = (double) (i);
            
            v_bin = -v_RV_max + delta_v*(i_doub + 0.5);
            
            outfile << k << "\t\t\t\t" << omegat << "\t\t\t\t" << v_bin << "\t\t\t\t\t" << B[i] << "\n";
            
        }

        
        
        outfile.close();
        
        // This outputs the luminosities
        
        //outfile << omegat << "\t\t\t\t\t" << pi*R*R*f0*c_x*(1.0 - (u_x/3.0)) << "\t\t\t\t\t" << L0 << "\t\t\t\t\t" << Lprime << "\t\t\t\t\t" << (Lprime - L0) / L0 << "\t\t\t\t\t" << dL_f/L0 << "\t\t\t\t\t" << dL_dot/L0<< "\n";
        
        
        
    }
    

    
    
}



