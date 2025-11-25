// Common Functions for Radially Symmetric GPE
#pragma once
#include <vector>
#include <cmath>
#include "parameters.h"

namespace cGPE {

    // Infenitesimal spatial step size

inline double Inf_spat_step(double R, int N)
    {
        return R/(N-1);
    }

    // Infenitesimal temporal step size

inline double Inf_time_step(double dr)
    {
        return std::pow(dr,2)*0.5;
    }

    // Grid

inline std::vector<double> grid(double R, double dr, int N)
    {
        std::vector<double> r(N);
        for (int i=0; i<N; i++){
            r[i] = i*dr;
        }
        return r;
    }

    // Initial Wave Function

inline std::vector<double> Initial_wavefunction(const std::vector<double> &r, double R)
    {
        std::vector<double> psi(r.size());
        for (int i=0; i<r.size(); i++){
            psi[i] = r[i]*std::exp(-std::pow(r[i],2));
        }
        return psi;
    }

    // Potential Function

inline std::vector<double> Potential_function(const std::vector<double> &r)
    {
        std::vector<double> V(r.size());
        for (int i=0; i<r.size(); i++){
            V[i] = 0.5 * std::pow(r[i],2);
        }
        return V;
    }

    // Normalization for 3-dim radially symmetric mathematical wave function (/psi = r*phi)

inline void normalize_symmetric(std::vector<double> &psi, double dr)
    {
        double norm = 0.0;
        for (int i=0; i<psi.size(); i++){
            norm += psi[i]*psi[i]*dr;
        }
        norm = std::sqrt(4 * cGPE::pi*norm);
        for (int i=0; i<psi.size(); i++){
            psi[i] = psi[i]/norm;
        }
    }

    // Expectation value of Energy

inline double energy_expectation(const std::vector<double> &psi,
                             const std::vector<double> &V, 
                             const std::vector<double>& r,
                             double g, double dr)
    {
        double E_kin = 0.0;
        double E_pot = 0.0;
        double E_int = 0.0;

        for(int i =1; i<psi.size()-1; i++){
            E_kin += 0.5 * 4 * cGPE::pi * std::pow((psi[i+1]-psi[i-1])/(dr*2),2) * dr;
            E_pot += 4 * cGPE::pi * V[i] * std::pow(psi[i],2) * dr;
            if (r[i] != 0){
                E_int += 2 * cGPE::pi * g * std::pow(psi[i],4) / (r[i]*r[i]) * dr; // we do not include interaction at r=0
            }
        }
        return E_kin + E_pot + E_int;
    }

}