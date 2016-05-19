
#ifndef GENERIC_TOOLS_H
#define GENERIC_TOOLS_H

void checkParameterSanity(IncomprFlowParam<T> param)
{
	// Check lattice viscosity [0.01, 0.45]
    if(param.getLatticeNu() < 0.01 || param.getLatticeNu() > 0.45)
        pcout << "(WARNING) lattice viscosity is not in the stable range for LBM [0.01, 0.45]!" << std::endl;

    // Check for lattice velocity to ensure low Courant number (LBM is explicit afterall...)
    if(param.getLatticeU() > 0.1)
        pcout << "(WARNING) lattice velocity is too high [>0.1]!" << std::endl;
}

#endif // GENERIC_TOOLS_H