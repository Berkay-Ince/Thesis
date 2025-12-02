// Run for radially symmetric GPE
#include <iostream>
#include <vector>
#include "Methods.h"
#include "parameters.h"
#include "Common_functions.h"

int main()
{
    //Parameters & Initials
    cGPE::params par;
    double dr = cGPE::Inf_spat_step(par.R, cGPE::N_default); 
    double dt = cGPE::Inf_time_step(dr);
    std::vector<double> r = cGPE::grid(par.R, dr, cGPE::N_default);
    std::vector<double> initial_psi = cGPE::Initial_wavefunction(r,par.R);
    std::vector<double> V = cGPE::Potential_function(r);
    //Time Evolution using Forward Euler Method
    /*std::vector<double> psi_FE = initial_psi; // Copy initial wavefunction
    for (int iter = 0; iter < par.max_iter; ++iter) {
        cGPE::FE_Method(psi_FE, V, r, dr, dt, par.g);
        cGPE::normalize_symmetric(psi_FE, dr);
    }
    //Output final result
    double E_final_FE = cGPE::energy_expectation(psi_FE, V, r, par.g, dr);
    cGPE::save_r_and_psi(r, psi_FE, "psi_FE.txt");
    std::cout << "Final Energy (FE Method): " << E_final_FE << std::endl;*/
    // Time Evolution using TSSP Method
    std::vector<double> psi_TSSP = initial_psi; // Copy initial
    for (int iter = 0; iter < par.max_iter; ++iter) {
        cGPE::TSSP_Step(psi_TSSP, V, r, dr, dt, par.g);
        cGPE::normalize_symmetric(psi_TSSP, dr);
    }
    // Output final result
    double E_final_TSSP = cGPE::energy_expectation(psi_TSSP, V, r, par.g, dr);
    cGPE::save_r_and_psi(r, psi_TSSP, "psi_TSSP.txt");
    std::cout << "Final Energy (TSSP Method): " << E_final_TSSP << std::endl;
    return 0;
}