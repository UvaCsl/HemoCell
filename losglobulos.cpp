/* This file is part of the Palabos library.
 * Copyright (C) 2009, 2010 Jonas Latt
 * E-mail contact: jonas@lbmethod.org
 * The most recent release of Palabos can be downloaded at 
 * <http://www.lbmethod.org/palabos/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <limits>

#include "palabos3D.h"
#include "palabos3D.hh"

#include "immersedWallParticle3D.h"
#include "immersedWallParticle3D.hh"
#include "immersedWallParticleFunctional3D.h"
#include "immersedWallParticleFunctional3D.hh"
#include "immersedWallParticleVtk3D.h"
#include "immersedWallParticleVtk3D.hh"
#include "shellModel3D.h"
#include "shellModel3D.hh"

#include "ficsionInit.hh"
#include "immersedCells3D.h"
#include "immersedCells3D.hh"


using namespace plb;
using namespace std;

typedef double T;
typedef Array<T,3> Velocity;
#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor

plint borderWidth     = 1;  // Because Guo acts in a one-cell layer.
// Requirement: margin>=borderWidth.
plint extraLayer      = 0;  // Make the bounding box larger; for visualization purposes
                            //   only. For the simulation, it is OK to have extraLayer=0.
const plint extendedEnvelopeWidth = 2;  // Because Guo needs 2-cell neighbor access.
const plint particleEnvelopeWidth = 6;




/* ************* Class GetTensorFieldFromExternalVectorFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor, int nDim>
class GetTensorFieldFromExternalVectorFunctional3D : public BoxProcessingFunctional3D_LT<T,Descriptor, T, nDim> {
public:
    GetTensorFieldFromExternalVectorFunctional3D (
        int vectorStartsAt_ ) : vectorStartsAt(vectorStartsAt_)
    {
        PLB_ASSERT( vectorStartsAt+nDim <=
        Descriptor<T>::ExternalField::numScalars );
    }
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice, TensorField3D<T,nDim>& tensor) {
        Dot3D offset = computeRelativeDisplacement(lattice, tensor);
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint oX = iX + offset.x;
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint oY = iY + offset.y;
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    plint oZ = iZ + offset.z;
                    Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                    Array<T,nDim> externalVector; 
                    
                    for (plint iD=0; iD<nDim; ++iD) {
                        externalVector[iD] = *cell.getExternal(vectorStartsAt+iD);
                    }
                    tensor.get(oX,oY,oZ) = externalVector;
                }
            }
        }
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }
    virtual GetTensorFieldFromExternalVectorFunctional3D<T,Descriptor,nDim>* clone() const {
        return new GetTensorFieldFromExternalVectorFunctional3D<T,Descriptor,nDim>(*this);
    }
    
private:
    int vectorStartsAt;
};

void iniLattice( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                 IncomprFlowParam<T> const& parameters,
                 OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    
    Box3D top    = Box3D(0,    nx-1, 0, ny-1, 0,    0);
    Box3D bottom = Box3D(0,    nx-1, 0, ny-1, nz-1, nz-1);

    Box3D left   = Box3D(0, nx-1, 0,    0,    1, nz-2);
    Box3D right  = Box3D(0, nx-1, ny-1, ny-1, 1, nz-2);

    Box3D inlet  = Box3D(0,    0,    1,    ny-2, 1, nz-2);
    Box3D outlet = Box3D(nx-1, nx-1, 1,    ny-2, 1, nz-2);

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, outlet );

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right );

    setBoundaryVelocity(lattice, inlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));
    setBoundaryVelocity(lattice, outlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));

    setBoundaryVelocity(lattice, top, Array<T,3>(0.0,0.0,0.0));
    setBoundaryVelocity(lattice, bottom, Array<T,3>(0.0,0.0,0.0));
    setBoundaryVelocity(lattice, left, Array<T,3>(0.0,0.0,0.0));
    setBoundaryVelocity(lattice, right, Array<T,3>(0.0,0.0,0.0));

    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), SquarePoiseuilleDensityAndVelocity<T>(parameters, NMAX));
    
    setExternalVector( lattice, lattice.getBoundingBox(), 
                       DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T,DESCRIPTOR<T>::d>(0.0,0.0,0.0));

    lattice.initialize();
}


//template<typename T, template<typename U> class Descriptor>
//void createImmersedWallParticles (
//        MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
//        TriangleBoundary3D<T>& boundary, plint tag, plint numPartsPerBloodCell )
//{
//    boundary.pushSelect(0,1);
//    std::vector<MultiBlock3D*> particleArg;
//    particleArg.push_back(&particleField);
//    applyProcessingFunctional (
//        new CreateTaggedImmersedWallParticle3D<T,Descriptor>(boundary,tag,numPartsPerBloodCell),
//        particleField.getBoundingBox(), particleArg );
//    boundary.popSelect();
//}
//
//template<typename T, template<typename U> class Descriptor>
//void deleteBloodCell(MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
//                   const Box3D &outlet, std::vector<plint> &numParts,
//                   std::vector<plint> &tags, TriangleBoundary3D<T> &bloodCells,
//                   std::vector<Array<T,3> > &centers, std::vector<plint > &radii )
//{
//    bool erased = false;
//    for (pluint iA = 0; iA < tags.size(); ++iA) {
//        // count all particles of a certain tag in a buffer zone
//        plint numPartsPerTag = countParticles(particleField,outlet,tags[iA]);
//        // if all the particle of a certain tag are in a buffer zone
//        if (numPartsPerTag == numParts[iA] && numPartsPerTag > 0) {
//            // then delete these particles in all the buffer zones
//            plint before = countParticles(particleField,outlet,tags[iA]);
//            std::vector<MultiBlock3D*> particleArg;
//            particleArg.push_back(&particleField);
//
//            applyProcessingFunctional (
//                new AbsorbTaggedParticlesFunctional3D<T,DESCRIPTOR>(tags[iA]),
//                    outlet, particleArg );
//
//            plint after = countParticles(particleField,outlet,tags[iA]);
//
//            pcout << "erased particles = " << before << ", " << after << std::endl;
//
////             delete bloodCells[iA];                     // delete the iA-th pointer mesh
////             numParts.erase(numParts.begin()+iA);        // erase the iA-th number of particles
////             tags.erase(tags.begin()+iA);                // erase the iA-th tag
////             bloodCells.erase(bloodCells.begin()+iA);  // erase the iA-th mesh
////             centers.erase(centers.begin()+iA);          // delete the iA-th center
////             radii.erase(radii.begin()+iA);              // delete the iA-th radius
////             --iA;
//
//            erased = true;
//        }
//    }
//
//    if (erased) {
//        pcout << "Particles absorbed : Number of particles per tag." << std::endl;
//        for (pluint iA = 0; iA < centers.size(); ++iA) {
//            pcout << tags[iA] << " , " <<  countParticles(particleField, particleField.getBoundingBox(), tags[iA]) << std::endl;
//        }
//    }
//}
//
//template<typename T, template<typename U> class Descriptor>
//bool generateBloodCells(MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
//                      const Box3D &inlet, std::vector<plint> &tags, TriangleBoundary3D<T> &bloodCells,
//                      plint numPartsPerBloodCell, plint numOfBloodCellsPerInlet, plint &slice )
//{
//    bool created = false;
//    plint numPartsPerTag = 0;
//    for (pluint iA = 0; iA < tags.size(); ++iA) {
//        // count all particles of a certain tag in a buffer zone
//        numPartsPerTag += countParticles(particleField,inlet,tags[iA]);
//    }
//
//    std::vector<plint> newTags;
//    for (plint iA = slice*numOfBloodCellsPerInlet; iA < (slice+1)*numOfBloodCellsPerInlet; ++iA) {
//        newTags.push_back(tags[iA]);
//    }
//
//    if (numPartsPerTag == 0) {
//        createBloodCells(bloodCells, newTags, numPartsPerBloodCell, particleField);
//        created = true;
//        ++slice;
//    }
//
//    if (created) {
//        pcout << "Particles created : Number of particles per tag." << std::endl;
//        for (pluint iA = 0; iA < newTags.size(); ++iA) {
//            pcout << newTags[iA] << " , " <<  countParticles(particleField, particleField.getBoundingBox(), newTags[iA]) << std::endl;
//        }
//    }
//    return created;
//}
//
//void createBloodCells(TriangleBoundary3D<T> &bloodCells,
//                    const std::vector<plint> &tags,
//                    plint numPartsPerBloodCell,
//                    MultiParticleField3D<DenseParticleField3D<T,DESCRIPTOR> > &immersedParticles)
//{
//    for (pluint iA = 0; iA < tags.size(); ++iA) {
//        createImmersedWallParticles(immersedParticles, bloodCells, tags[iA], numPartsPerBloodCell);
//    }
//}



template<class BlockLatticeT>
void writeVTK(BlockLatticeT& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    
    MultiTensorField3D<T,3> force(lattice);
    
    applyProcessingFunctional(new GetTensorFieldFromExternalVectorFunctional3D<T,DESCRIPTOR,3>(
        DESCRIPTOR<T>::ExternalField::forceBeginsAt), lattice.getBoundingBox(), lattice, force);
                              
    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<3,float>(force, "force", (T)1);
//     vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);
//     vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
}

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    global::IOpolicy().setStlFilesHaveLowerBound(true);
    global::IOpolicy().setLowerBoundForStlFiles(-1.);

    T shellDensity = 0.;
    T k_rest = 0.;
    T k_stretch = 0.;
    T k_shear = 0.;
    T k_bend = 0.;

    string paramXmlFileName;
    global::argv(1).read(paramXmlFileName);

    XMLreader document(paramXmlFileName);
    document["shellDensity"].read(shellDensity);
    document["k_rest"].read(k_rest);
    document["k_stretch"].read(k_stretch);
    document["k_shear"].read(k_shear);
    document["k_bend"].read(k_bend);

    IncomprFlowParam<T> parameters(
            0.01, // u 
            100., // Re
            60,   // N
            10.,        // lx
            1.,        // ly
            1.         // lz
    );

    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNz();
    
    MultiBlockLattice3D<T, DESCRIPTOR> lattice(
        defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz,
                                                            extendedEnvelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
        new GuoExternalForceBGKdynamics<T,DESCRIPTOR>(parameters.getOmega()));
    

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    pcout << std::endl << "Initializing lattice" << std::endl;
    iniLattice(lattice, parameters, *boundaryCondition);

    MultiBlockManagement3D const& latticeManagement(lattice.getMultiBlockManagement());
    MultiBlockManagement3D particleManagement (
            latticeManagement.getSparseBlockStructure(),
            latticeManagement.getThreadAttribution().clone(),
            particleEnvelopeWidth,
            latticeManagement.getRefinementLevel() );

    MultiParticleField3D<DenseParticleField3D<T,DESCRIPTOR> > immersedParticles (
            particleManagement,
            defaultMultiBlockPolicy3D().getCombinedStatistics() );

    plint x0 = 22; plint x1 = nx-x0;

    Box3D inlet(0,  x0,  0, ny-1, 0, nz-1);
    Box3D outlet(x1, nx-1,0, ny-1, 0, nz-1);

    std::vector<plint> pos;
    pos.push_back(30); pos.push_back(10);pos.push_back(50);
    
    std::vector<Array<T,3> > centers;
    std::vector<plint > radii;
    plint nMax = 1;
    for (plint iN = 0; iN < nMax; ++iN) { // create 40 * inlet amount of particles
        for (pluint iA = 0; iA < pos.size(); ++iA) {
            for (pluint iB = 0; iB < pos.size(); ++iB) {
                centers.push_back(Array<T,3>(10+iN,pos[iA],pos[iB]));
                radii.push_back(3.91);
            }
        }
    }

    plint numOfBloodCellsPerInlet = pos.size() * pos.size(); // number used for the generation of bloodCells at inlet

    std::vector<plint> tags;
    plint numPartsPerBloodCell = 0; plint slice = 0; // number of particles per tag and number of slice of created particles
    TriangleBoundary3D<T> bloodCells = createCompleteMesh(centers, radii, tags, numPartsPerBloodCell);
    bool created = generateBloodCells(immersedParticles, inlet, tags, bloodCells, numPartsPerBloodCell, numOfBloodCellsPerInlet, slice );

    std::vector<plint> numParts(tags.size());
    for (pluint iA = 0; iA < tags.size(); ++iA) {
        numParts[iA] = countParticles(immersedParticles, immersedParticles.getBoundingBox(), tags[iA]);
    }
    
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&immersedParticles);

    std::vector<MultiBlock3D*> particleLatticeArg;
    particleLatticeArg.push_back(&immersedParticles);
    particleLatticeArg.push_back(&lattice);

    SpringModel3D<T> springModel(shellDensity, k_rest, k_stretch, k_shear, k_bend );

    plint maxIter = 100000;
    plint imageIter = 10;
    pcout << std::endl << "Starting simulation" << std::endl;
    global::timer("sim").start();
    for (plint i=0; i<maxIter; ++i) {

        applyProcessingFunctional ( // advance particles in time according to a velocity, acceleration, ...
            new AdvanceParticlesFunctional3D<T,DESCRIPTOR>,
            immersedParticles.getBoundingBox(), particleArg );

        applyProcessingFunctional ( // copy fluid velocity on particles
            new FluidVelocityToImmersedWall3D<T,DESCRIPTOR>(),
            immersedParticles.getBoundingBox(), particleLatticeArg);

        bloodCells.pushSelect(0,1);
        applyProcessingFunctional ( // update mesh position
            new CopyParticleToVertex3D<T,DESCRIPTOR>(bloodCells.getMesh()),
            immersedParticles.getBoundingBox(), particleArg);
        bloodCells.popSelect();

        applyProcessingFunctional ( // compute force applied on the particles by springs
            new ComputeImmersedElasticForce3D<T,DESCRIPTOR> (
                bloodCells, springModel.clone() ), // used because pushSelect is not used
            immersedParticles.getBoundingBox(), particleArg );
        
        setExternalVector( lattice, lattice.getBoundingBox(), 
                           DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T,DESCRIPTOR<T>::d>(0.0,0.0,0.0));
        
        applyProcessingFunctional ( // compute force applied on the particles by springs
            new ForceToFluid3D<T,DESCRIPTOR> (), // used because pushSelect is not used
                immersedParticles.getBoundingBox(), particleLatticeArg );
        
        lattice.collideAndStream();
        
//         applyProcessingFunctional (
//             new ComputeFluidForceOnParticle3D<T,DESCRIPTOR> (
//                 boundary, springModel.getDensity(), flowType ),
//             dynamicMeshParticles.getBoundingBox(), fluidForceArg );

        deleteBloodCell(immersedParticles, outlet, numParts, tags, bloodCells, centers, radii );
        if (slice < 1) {
            bool created = generateBloodCells(immersedParticles, inlet, tags, bloodCells, numPartsPerBloodCell, numOfBloodCellsPerInlet, slice );
        }

        if (i%imageIter==0) {
            pcout << "totParticles = " << countParticles(immersedParticles, immersedParticles.getBoundingBox()) << std::endl;
            pcout << "Write Particle VTK" << std::endl;
            std::vector<std::string> force_scalarNames;
            force_scalarNames.push_back("pressure");
            force_scalarNames.push_back("wss");
            std::vector<std::string> velocity_scalarNames;
            std::vector<std::string> force_vectorNames;
            force_vectorNames.push_back("force");
            std::vector<std::string> velocity_vectorNames;
            velocity_vectorNames.push_back("velocity");
//             for (pluint iA = 0; iA < centers.size(); ++iA) {
                bool dynamicMesh = true;
                plint tag = -1; // Take all triangles.

                // serialize the particle information to write them.
                // a correspondance between the mesh and the particles is made.
                writeImmersedSurfaceVTK (
                        bloodCells,
                        *getParticlePosAndVelocity(immersedParticles),
                        velocity_scalarNames, velocity_vectorNames,
                        global::directories().getOutputDir()+createFileName("RBC",i,6)+".vtk", dynamicMesh, tag );
    //             bool printHeader = true;
    //             VtkImageOutput3D<T> vtkOut(createFileName("volume",i,6), boundary.getDx(), boundary.getPhysicalLocation());
    //             vtkOut.writeData<float>(*boundaryCondition->computePressure(), "p", 1.);
    //             vtkOut.writeData<float>(*boundaryCondition->computeVelocityNorm(), "uNorm", 1.);
//             }

            writeVTK(lattice, parameters, i);
        }
    }

}
