#include "ISOP2P1.h"

void ISOP2P1::config(std::string _config_file)
{
    std::string trash;
    std::ifstream input(_config_file.c_str());
    input >> trash >> trash >> mesh_file;
    input >> trash >> trash >> l_tol;
    input >> trash >> trash >> l_Euler_tol;
    input >> trash >> trash >> n_tol;
    input >> trash >> trash >> viscosity;
    input >> trash >> trash >> t0;
    input >> trash >> trash >> t1;
    input >> trash >> trash >> t;
    input >> trash >> trash >> dt;
    input >> trash >> trash >> CFL;
    input >> trash >> trash >> n_method;
    input >> trash >> trash >> time_step_control;
    input >> trash >> trash >> Stokes_init;
    input >> trash >> trash >> NS_init;
    input >> trash >> trash >> scheme;
    input >> trash >> trash >> body_force;
    input >> trash >> trash >> angle;
    input.close();
};

