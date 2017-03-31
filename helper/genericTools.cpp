#ifndef GENERIC_TOOLS_CPP
#define GENERIC_TOOLS_CPP

#include "hemocell.h"

void checkParameterSanity(T nu_lbm, T u_max_lbm)
{
	// Check lattice viscosity [0.01, 0.45]
    if(nu_lbm < 0.01 || nu_lbm > 0.45)
        pcout << "(WARNING!!!) lattice viscosity [" << nu_lbm << "] is not in the stable range for LBM [0.01, 0.45]!" << std::endl;

    // Check for lattice velocity to ensure low Courant number (LBM is explicit afterall...)
    if(u_max_lbm > 0.1)
        pcout << "(WARNING!!!) lattice velocity [" << u_max_lbm << "] is too high [>0.1]!" << std::endl;
}

void printHeader()
{
	pcout << " _   _   ____   __  __   _____    ___   ____   __     __ " << endl;
	pcout << "( )_( ) ( ___) (  \\/  ) (  _  )  / __) ( ___) (  )   (  )" << endl;
	pcout << " ) _ (   )__)   )    (   )(_)(  ( (__   )__)   )(__   )(__ " << endl;
	pcout << "(_) (_) (____) (_/\\/\\_) (_____)  \\___) (____) (____) (____) " << endl;
	pcout << "                         v." << VERSION_MAJOR << "." << VERSION_MINOR << endl;
	pcout << endl;
}

#endif // GENERIC_TOOLS_H
