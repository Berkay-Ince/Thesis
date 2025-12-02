#pragma once
#include <vector>
#include <cmath>
#include "parameters.h"
#include "Common_functions.h"
#include <fftw3.h>

namespace cGPE {

    // Forward Euler Time Stepping Method for Radially Symmetric GPE

inline void FE_Method(std::vector<double> &psi,
                                             const std::vector<double> &V,
                                             const std::vector<double> &r, 
                                             double dr, double dt, double g)
    {
        std::size_t N = psi.size();
        std::vector<double> laplacian_psi = cGPE::laplacian(psi, dr);
        for (int i = 0; i <N; i++){
            psi[i] -= dt * (-0.5 * laplacian_psi[i] + V[i] * psi[i] + g * psi[i] *psi[i] * psi[i] / (r[i]*r[i]+ 1e-30));
        }
        psi [0] = 0.0; // Boundary condition at r=0
        psi [N-1] = 0.0; // Boundary condition at r=R
    }

    // Time-Splitting Spectral Method for Radially Symmetric GPE 

    // Kinetic Part

inline void TSSP_Kinetic(std::vector<double> &psi,
                                    double dr, double dt)
    {
    
    const std::size_t N = psi.size();
    
    //Allocate input (real) and output (complex)
    double* in = (double*) fftw_malloc(sizeof(double) * N); // fftw_malloc is used for alignment it can be written as std::vector but malloc is specilized for fftw
    fftw_complex* psi_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2 + 1)); 
    
    // std::memcpy(in, psi.data(), N * sizeof(double)); there is thing like that but I do not now about memcpy thats why for a while I will use for loop.
    for (int i = 0; i < N; ++i)
        in[i] = psi[i]; 
    
    // FFTW_ESTIMATE is for fast planning, FFTW_MEASURE is for more accurate planning but slower (doing benchmarking),
    // FFTW_PATIENT is very accurate but very slow planning (testing everything, could take minutes to hours)

    fftw_plan plan = fftw_plan_dft_r2c_1d((int)N, in, psi_k, FFTW_ESTIMATE); 
    
    fftw_execute(plan);

    double dk = 2 * cGPE::pi / (dr * N);

    // Apply kinetic evolution in Fourier space

    for (int m = 0; m < N/2 + 1;++m) {
        double k = m * dk;
        double damping = std::exp(-0.5 * k * k * dt);
        psi_k[m][0] *= damping; // Real part
        psi_k[m][1] *= damping; // Imaginary part
    }

    // Inverse FFT to get back to real space

    fftw_plan inv_plan = fftw_plan_dft_c2r_1d((int)N, psi_k, in, FFTW_ESTIMATE);
    fftw_execute(inv_plan);
    for (int i = 0; i < N; ++i) {
        psi[i] = in[i] / (double)N; // Normalize the inverse FFT output
    }
    // Clean up
    fftw_destroy_plan(plan);
    fftw_destroy_plan(inv_plan);
    fftw_free(in);
    fftw_free(psi_k);
}

inline void TSSP_Step(std::vector<double> &psi,
                                    const std::vector<double> &V,
                                    const std::vector<double> &r,
                                    double dr, double dt, double g)
{
    const std::size_t N = psi.size();

    // First step
    for(int i = 0; i<N; i++){
        double damping_part = std::exp(-dt * (V[i] + g * psi[i] * psi[i] / (r[i]*r[i]+ 1e-30))*0.5);
        psi[i] *= damping_part;
    }

    // Kinetic evolution
    TSSP_Kinetic(psi,dr,dt);

    // Second step
    for(int i = 0; i<N; i++){
        double damping_part = std::exp(-dt * (V[i] + g * psi[i] * psi[i] / (r[i]*r[i]+ 1e-30))*0.5);
        psi[i] *= damping_part;
    }
}

}