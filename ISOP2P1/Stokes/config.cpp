#include "ISOP2P1.h"

void ISOP2P1::config(std::string _config_file)
{
    std::string trash;
    std::ifstream input(_config_file.c_str());
    input >> trash >> trash >> mesh_file;
    input >> trash >> trash >> viscosity;
    input >> trash >> trash >> body_force;
    input >> trash >> trash >> angle;
    input.close();
};

