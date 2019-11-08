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

/*
Lees-Edwards boundary condition class which makes use of two data processors
for calculating and writing of the populations according to the Lees-Edwards algorithm.

Populations are read and writen between two static 2D arrays (topPopulations & bottomPopulations).
The memory address for the current Lees-Edwards displacement is shared with the hemoCell class.

The 2D arrays and Lees-Edwards displacement are static in order to guarantee data consistency
after subdivision of the code over the individual lattice blocks.

The velocity in the system comes from setting the macroscopic velocity in the boundary nodes

NOTE: currently the Lees-Edwards algorithm works only in one direction (z-direction)

* @author Daan van Ingen
*/

#ifndef LEESEDWARDSBC
#define LEESEDWARDSBC

#include "hemocell.h"
#include "palabos3D.h"
#include "palabos3D.hh"

namespace hemo
{

/*  Calculate the streaming of the populations according to the Lees-Edwards algorithm
    Inherits from data processing functional which handles the subdivision of the code
    onto individual lattice blocks
*/
template <typename T, template <typename U> class Descriptor>
class LeesEdwardsBCGetPopulations : public plb::BoxProcessingFunctional3D_L<T, Descriptor>
{
public:
    plint nx; // dimension in x direction
    plint nz; // dimension in z direction
    double * LEcurrentDisplacement; // Current Lees-Edwards displacement
    std::vector<std::vector<std::vector<T>>> * topPopulations; // Data structure containing the top populations
    std::vector<std::vector<std::vector<T>>> * bottomPopulations; // Data structure containing the bottom populations
    T topVelocity; // Macroscopic top velocity
    T bottomVelocity; // Macroscopic bottom velocity

    LeesEdwardsBCGetPopulations(plint nx, plint nz, T topVelocity, T bottomVelocity, double * LEcurrentDisplacement,
            std::vector<std::vector<std::vector<T>>> * topPopulations,
            std::vector<std::vector<std::vector<T>>> * bottomPopulations)
    {
        this->nx = nx;
        this->nz = nz;
        this->LEcurrentDisplacement = LEcurrentDisplacement;
        this->topVelocity = topVelocity;
        this->bottomVelocity = bottomVelocity;
        this->topPopulations = topPopulations;
        this->bottomPopulations = bottomPopulations;
    }

    /* Convert Palabos array to vector

     * @param populations Palabos array containing lattice node populations
     * @return vector containing lattice node populations
     */
    std::vector<T> populationsArrayToVector(plb::Array<T, Descriptor<T>::q> populations) {
        std::vector<T> pop(Descriptor<T>::q);
        for(plint i = 0; i < Descriptor<T>::q; i++) {
            pop[i] = populations[i];
        }
        return pop;
    }

    /* Non-negative modulo (-a % b gives a negative number in C++)
     */
    plint modNonNegative(plint a, plint b)
    {
        return (a % b + b) % b;
    }

    /* Process function containg the code to be executed onto the individual lattice block.
       Set macroscopic velocity
       Calculate new populations according to the Lees-Edwards algorithm

     * @param domain Domain on which to execute code
     * @param lattice individual block lattice
     */
    virtual void process(plb::Box3D domain, plb::BlockLattice3D<T, Descriptor> &lattice)
    {
        Dot3D absoluteOffset = lattice.getLocation(); // Location of this block lattice {x, y, z}

        // The portion of the population streamed from node S1 or S2 is determined by the overlap with the reference node
        double gfrac = std::fmod((*this->LEcurrentDisplacement), 1.0); 
        plint globalX, globalY;
	plint globalZ0 = absoluteOffset.z + domain.z0;
	plint globalZ1 = absoluteOffset.z + domain.z1;
        Cell<T, Descriptor> curCell;
        Cell<T, Descriptor> s1Cell;
        Cell<T, Descriptor> s2Cell;

        // Variables for setting the macroscopic velocity
        T rhoBar;
        plb::Array<T, Descriptor<T>::d> j;

        if (globalZ1 == this->nz - 1) // Top
        {
            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                globalX = absoluteOffset.x + iX;
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    globalY = absoluteOffset.y + iY;
                    curCell = lattice.get(iX, iY, domain.z1);

                    // Set macroscopic velocity
                    curCell.getDynamics().computeRhoBarJ(curCell, rhoBar, j);
                    curCell.getDynamics().collideExternal(curCell, rhoBar, plb::Array<T,3>(this->topVelocity, 0.0, 0.0), T(), lattice.getInternalStatistics());
                    (*this->topPopulations)[globalX][globalY] = populationsArrayToVector(curCell.getRawPopulations());

                    // Determine reference nodes coordinates according to the current displacement
                    plint s1 = domain.x0 + (plint)modNonNegative(std::ceil((*this->LEcurrentDisplacement) + globalX), this->nx);
                    plint s2 = domain.x0 + (plint)modNonNegative(std::floor((*this->LEcurrentDisplacement) + globalX), this->nx);

                    // Get the nodes
                    s1Cell = lattice.get(s1, iY, domain.z1);
                    s2Cell = lattice.get(s2, iY, domain.z1);

                    // Write resulting LE stream to data structure
                    (*this->topPopulations)[globalX][globalY][3] = gfrac * s1Cell[3] + (1 - gfrac) * s2Cell[3];
                    (*this->topPopulations)[globalX][globalY][6] = gfrac * s1Cell[16] + (1 - gfrac) * s2Cell[16];
                    (*this->topPopulations)[globalX][globalY][8] = gfrac * s1Cell[8] + (1 - gfrac) * s2Cell[8];
                    (*this->topPopulations)[globalX][globalY][16] = gfrac * s1Cell[6] + (1 - gfrac) * s2Cell[6];
                    (*this->topPopulations)[globalX][globalY][18] = gfrac * s1Cell[18] + (1 - gfrac) * s2Cell[18];
                }
            }
        }
        if (globalZ0 == 0) // Bottom
        {
            for (plint iX = domain.x0; iX <= domain.x1; ++iX)
            {
                globalX = absoluteOffset.x + iX;
                for (plint iY = domain.y0; iY <= domain.y1; ++iY)
                {
                    globalY = absoluteOffset.y + iY;
                    curCell = lattice.get(iX, iY, domain.z0);

                    // Set macroscopic velocity
                    curCell.getDynamics().computeRhoBarJ(curCell, rhoBar, j);
                    curCell.getDynamics().collideExternal(curCell, rhoBar, plb::Array<T,3>(this->bottomVelocity, 0.0, 0.0), T(), lattice.getInternalStatistics());
                    (*this->bottomPopulations)[globalX][globalY] = populationsArrayToVector(curCell.getRawPopulations());

                    // Determine reference nodes coordinates according to the current displacement
                    plint s1 = domain.x0 + (plint)modNonNegative(std::floor(-(*this->LEcurrentDisplacement) + globalX), this->nx);
                    plint s2 = domain.x0 + (plint)modNonNegative(std::ceil(-(*this->LEcurrentDisplacement) + globalX), this->nx);

                    // Get the nodes
                    s1Cell = lattice.get(s1, iY, domain.z0);
                    s2Cell = lattice.get(s2, iY, domain.z0);

                    // Write resulting LE stream to data structure
                    (*this->bottomPopulations)[globalX][globalY][7] = gfrac * s1Cell[15] + (1 - gfrac) * s2Cell[15];
                    (*this->bottomPopulations)[globalX][globalY][9] = gfrac * s1Cell[9] + (1 - gfrac) * s2Cell[9];
                    (*this->bottomPopulations)[globalX][globalY][12] = gfrac * s1Cell[12] + (1 - gfrac) * s2Cell[12];
                    (*this->bottomPopulations)[globalX][globalY][15] = gfrac * s1Cell[7] + (1 - gfrac) * s2Cell[7];
                    (*this->bottomPopulations)[globalX][globalY][17] = gfrac * s1Cell[17] + (1 - gfrac) * s2Cell[17];
                }
            }
        }
        // cout << global::mpi().getRank() << " " << bottomPopulations[5][5][10] << endl;
    }

    // Clone class to individual lattice blocks
    virtual LeesEdwardsBCGetPopulations<T, Descriptor> *clone() const
    {
        return new LeesEdwardsBCGetPopulations<T, Descriptor>(*this);
    }

    // Return the type of data that is modified
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
    }

    // Only need the bulk of the domain
    virtual plb::BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulk;
    }
};

/*  Write the streaming of the populations according to the Lees-Edwards algorithm
    Inherits from data processing functional which handles the subdivision of the code
    onto individual lattice blocks
*/
template <typename T, template <typename U> class Descriptor>
class LeesEdwardsBCSetPopulations : public plb::BoxProcessingFunctional3D_L<T, Descriptor>
{
public:
    plint nz; // dimension in z direction
    std::vector<std::vector<std::vector<T>>> * topPopulations; // Data structure containing the top populations
    std::vector<std::vector<std::vector<T>>> * bottomPopulations; // Data structure containing the bottom populations

    LeesEdwardsBCSetPopulations(plint nz,
            std::vector<std::vector<std::vector<T>>> * topPopulations,
            std::vector<std::vector<std::vector<T>>> * bottomPopulations) {
        this->nz = nz;
        this->topPopulations = topPopulations;
        this->bottomPopulations = bottomPopulations;
    }

    /* Convert vector to Palabos array

     * @param populations vector containing lattice node populations
     * @return Palabos array containing lattice node populations
     */
    plb::Array<T, Descriptor<T>::q> populationsVectorToArray(std::vector<T> populations)
    {
        plb::Array<T, Descriptor<T>::q> pop;
        for (plint i = 0; i < DESCRIPTOR<T>::q; i++)
        {
            pop[i] = populations[i];
        }
        return pop;
    }

    /* Process function containg the code to be executed onto the individual lattice block.
       Write new populations to lattice.

     * @param domain Domain on which to execute code
     * @param lattice individual block lattice
     */
    virtual void process(plb::Box3D domain, plb::BlockLattice3D<T, Descriptor> &lattice)
    {
        Dot3D absoluteOffset = lattice.getLocation(); // Location of this block lattice {x, y, z}

        plint globalZ0 = absoluteOffset.z + domain.z0;
        plint globalZ1 = absoluteOffset.z + domain.z1;
        plint globalX, globalY;
        plb::Array<T, Descriptor<T>::q> populations;

        if (globalZ1 == this->nz - 1) { // Top
            for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
                globalX = absoluteOffset.x + iX;
                for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                    globalY = absoluteOffset.y + iY;
                    populations = populationsVectorToArray((*this->topPopulations)[globalX][globalY]);
                    lattice.get(iX, iY, domain.z1).setPopulations(populations);
                }
            }
        }
        if (globalZ0 == 0) { // Bottom
            for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
                globalX = absoluteOffset.x + iX;
                for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                    globalY = absoluteOffset.y + iY;
                    populations = populationsVectorToArray((*this->bottomPopulations)[globalX][globalY]);
                    lattice.get(iX, iY, domain.z0).setPopulations(populations);
                }
            }
        }
    }

    // Clone class to individual lattice blocks
    virtual LeesEdwardsBCSetPopulations<T, Descriptor> *clone() const
    {
        return new LeesEdwardsBCSetPopulations<T, Descriptor>(*this);
    }

    // Return the type of data that is modified
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }

    // Only need the bulk of the domain
    virtual plb::BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulk;
    }
};


/* Lees-Edwards boundary conditions base class responsible for the setup and initialization of the data processors
 */
template<typename T, template<class U> class Descriptor>
class LeesEdwardsBC
{
public:
    plb::MultiBlockLattice3D<T,Descriptor> &lattice;
    plint nx;      // Dimension in x direction
    plint ny;      // Dimension in y direction
    plint nz;      // Dimension in z direction
    // int direction; // TODO direction of the LE boundary x = 0, y = 1, z = 2
    double LEdisplacement; // Displacement per timestep
    static double LEcurrentDisplacement; // Current total displacement
    T dt; // Timestep size
    T topVelocity; // Macroscopic velocity top boundary layer
    T bottomVelocity; // Macroscopic velocity bottom boundary layer
    plint dataProcessorLevel; // level of operation of the Lees-Edwards data processors

    static std::vector<std::vector<std::vector<T>>> topPopulations; // Data structure containing the top populations
    static std::vector<std::vector<std::vector<T>>> bottomPopulations; // Data structure containing the bottom populations

    LeesEdwardsBC(plb::MultiBlockLattice3D<T, Descriptor> &lattice, T shearRate, T dt, double ** hemoCellLEcurrentDisplacement, plint dataProcessorLevel = 1) : lattice(lattice) {
        this->lattice = lattice;
        this->nx = lattice.getNx();
        this->ny = lattice.getNy();
        this->nz = lattice.getNz();
        this->dt = dt;
        this->LEdisplacement = shearRate * dt;
        T vHalf = (nz - 1) * shearRate * 0.5;
        this->topVelocity = -vHalf;
        this->bottomVelocity = vHalf;
        this->dataProcessorLevel = dataProcessorLevel;
        *hemoCellLEcurrentDisplacement = &LEcurrentDisplacement; // Let the hemocell LE displacement point to this memory address
    }

    void initialize()
    {
        // Uses the default periodicity because otherwise a form of bounceback boundary will be initialized per default
        this->lattice.periodicity().toggleAll(true);

        this->topPopulations.resize(this->nx, std::vector<std::vector<T>>(this->ny, std::vector<T>(Descriptor<T>::q)));
        this->bottomPopulations.resize(this->nx, std::vector<std::vector<T>>(this->ny, std::vector<T>(Descriptor<T>::q)));

        // set data processors        
        integrateProcessingFunctional(
            new LeesEdwardsBCGetPopulations<T, DESCRIPTOR>(this->nx, this->nz, this->topVelocity, this->bottomVelocity, &LEcurrentDisplacement, &topPopulations, &bottomPopulations),
            this->lattice.getBoundingBox(), this->lattice, this->dataProcessorLevel);

        integrateProcessingFunctional(
            new LeesEdwardsBCSetPopulations<T, DESCRIPTOR>(this->nz, &topPopulations, &bottomPopulations),
            this->lattice.getBoundingBox(), this->lattice, this->dataProcessorLevel + 1);
        }

    /* Update the current Lees-Edwards displacement

     * @param iter current hemocell iteration counter
     */
    void updateLECurDisplacement(unsigned int iter)
    {
        LEcurrentDisplacement = std::fmod(this->LEdisplacement * iter, (double)this->nx);
    }
};

#ifndef LEBC_STATICS
#define LEBC_STATICS
// Initialize static variables
// Lees-Edwards displacement
template <typename T, template <typename U> class Descriptor>
double LeesEdwardsBC<T, Descriptor>::LEcurrentDisplacement;
// Top populations data storage for separate writing and reading across different data processors
template <typename T, template <typename U> class Descriptor>
std::vector<std::vector<std::vector<T>>> LeesEdwardsBC<T, Descriptor>::topPopulations;
// Bottom populations data storage for separate writing and reading across different data processors
template <typename T, template <typename U> class Descriptor>
std::vector<std::vector<std::vector<T>>> LeesEdwardsBC<T, Descriptor>::bottomPopulations;
#else
#endif

}; // namespace hemo

#endif
