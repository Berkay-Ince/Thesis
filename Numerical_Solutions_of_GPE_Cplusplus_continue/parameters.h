// Parameters for Radially Symmetric GPE
#pragma once
#include <vector>

namespace cGPE {

struct params {
    
    int N = 4096;          // Number of spatial grid points
    double R = 6.0;       // Boundary radius
    double g = 0.0;      // Interaction strength
    double tol = 1e-10;     // Tolerance for numerical methods
    int max_iter = 5e5; // Maximum number of iterations
    
};

inline constexpr double pi = 3.14159265358979323846;

struct sweep_parameters
{
    std::vector<double> R_values{2.0,6.0,10.0}; // List of boundary radius
    std::vector<double> g_values{0.0,5.0,100.0}; // List of interaction strengths
};

}

