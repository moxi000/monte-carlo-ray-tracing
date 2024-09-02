# Monte Carlo Ray Tracing for Large Dust Aggregate Particles

This project investigates the radiative properties of large dust aggregate particles using the Monte Carlo ray tracing method. It simulates the interaction of light with aggregated dust particles, considering scattering, absorption, and reflection.

## Project Structure

- `ParticleMC.cpp`: Main program file implementing the Monte Carlo ray tracing method.
- `Function.h`: Header file containing auxiliary functions used in the simulation.
- `Vector.h`: Header file providing vector operations required for simulations.
- `shape_sphere.txt`: Input file containing the parameters of the spheres used in the simulation.
- `position_range.txt`: Input file specifying the position range or boundaries for the simulation.
- `MC.sln`, `MC.vcxproj`, etc.: Visual Studio solution and project files for building the application.
- Output data files (`F_angle.txt`, `f_function.txt`, `f_position.txt`, `L_I.txt`, `q_angle.txt`, `result.txt`): These files store the results of the simulation, such as scattering angles, phase functions, positions, and intensity distributions.

## Requirements

- **Compiler**: The project is set up using Visual Studio. Ensure you have Visual Studio installed with C++ support.
- **Libraries**: The code uses OpenMP for parallel processing. Ensure OpenMP support is enabled in your compiler settings.


## Authors

- **Liu Xiaochuan**
- **Tang Yanxia**
- **Zhu Keyong**
- **Huang Yong**

School of Aeronautic Science and Engineering, Beihang University (Beijing University of Aeronautics and Astronautics), Beijing, China.

Contact Email: huangy@buaa.edu.cn
