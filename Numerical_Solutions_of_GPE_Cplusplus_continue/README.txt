parameters.h includes   N_default,   R,  g,  tol,   max_iter,   pi
Common_function.h includes  Inf_spat_step, Inf_time_step,  grid, Initial_wavefunction, Potential_function,    normalize_symmetric, energy_expectation, laplacian, save_r_and_psi(this part can be better for the implement different parameters.)
Methods.h includes FE_Method, TSSP_Kinetic (here for fft, fftw library is used. Later, I am planning to get outside the fftw_plan function. For now, I was just looking at consistency with python code.), TSSP_Step
plot.ipynb it will be delevoped for more than one parameters and methods. 
main.cpp it is there for running methods.
