/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab 
in the University of Amsterdam. Any questions or remarks regarding this library 
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "genericFunctions.h"

void weakScaling(int Nx, int Ny, int Nz, int numberOfProcesses, vector<int> & newNxNyNz) {
    int fmod = int(log2(numberOfProcesses))%3;
    int fdiv = int((log2(numberOfProcesses)/3));
    int ffactor = pow(2.0,fdiv);
    newNxNyNz.clear();
    newNxNyNz.push_back(Nx * (ffactor * (1 + fmod>1) + int(pow(2.0,fdiv))));
    newNxNyNz.push_back(Ny * (ffactor * (1 + fmod>2) + int(pow(2.0,fdiv))));
    newNxNyNz.push_back(Nz * (ffactor * (1 + fmod>0)));
}

int renameFileToDotOld(std::string fName) {
    int renameStatus = 0;
    if (file_exists(fName)) {
        int renameStatus = rename(fName.c_str(), (fName + ".old").c_str());
        if (renameStatus != 0) {
            pcout << fName << " error." << std::endl;
        }
    }
    return renameStatus;
}



int do_mkdir(const char *path, mode_t mode)
{
    Stat            st;
    int             status = 0;

    if (stat(path, &st) != 0)
    {
        /* Directory does not exist. EEXIST for race condition */
        if (mkdir(path, mode) != 0 && errno != EEXIST)
            status = -1;
    }
    else if (!S_ISDIR(st.st_mode))
    {
        errno = ENOTDIR;
        status = -1;
    }

    return(status);
}

/**
** mkpath - ensure all directories in path exist
** Algorithm takes the pessimistic view and works top-down to ensure
** each directory in path exists, rather than optimistically creating
** the last element and working backwards.
*/
// mkpath(argv[i], 0777);
int mkpath(const char *path, mode_t mode)
{
    char           *pp;
    char           *sp;
    int             status;
    char           *copypath = strdup(path);

    status = 0;
    pp = copypath;
    while (status == 0 && (sp = strchr(pp, '/')) != 0)
    {
        if (sp != pp)
        {
            /* Neither root nor double slash in path */
            *sp = '\0';
            status = do_mkdir(copypath, mode);
            *sp = '/';
        }
        pp = sp + 1;
    }
    if (status == 0)
        status = do_mkdir(path, mode);
    free(copypath);
    return (status);
}

std::string zeroPadNumber(int num)
{
    std::ostringstream ss;
    ss << std::setw( 8 ) << std::setfill( '0' ) << num;
    return ss.str();
}

void checkParameterSanity(T nu_lbm, T u_max_lbm)
{
	// Check lattice viscosity [0.01, 0.45]
    if(nu_lbm < 0.01 || nu_lbm > 0.45)
        hlog << "(WARNING!!!) lattice viscosity [" << nu_lbm << "] is not in the stable range for LBM [0.01, 0.45]!" << std::endl;

    // Check for lattice velocity to ensure low Courant number (LBM is explicit afterall...)
    if(u_max_lbm > 0.1)
        hlog << "(WARNING!!!) lattice velocity [" << u_max_lbm << "] is too high [>0.1]!" << std::endl;
}

void printHeader()
{
	hlog << " _   _   ____   __  __   _____    ___   ____   __     __ " << endl;
	hlog << "( )_( ) ( ___) (  \\/  ) (  _  )  / __) ( ___) (  )   (  )" << endl;
	hlog << " ) _ (   )__)   )    (   )(_)(  ( (__   )__)   )(__   )(__ " << endl;
	hlog << "(_) (_) (____) (_/\\/\\_) (_____)  \\___) (____) (____) (____) " << endl;
	hlog << "                         v." << VERSION_MAJOR << "." << VERSION_MINOR << endl;
	hlog << endl;
}