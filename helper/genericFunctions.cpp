#ifndef FCN_GENERIC_FUNCTIONS_HH
#define FCN_GENERIC_FUNCTIONS_HH
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



inline bool file_exists (const std::string& name) {
    /* Checks if a file exists */
// http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
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

#endif // FCN_GENERIC_FUNCTIONS_HH
