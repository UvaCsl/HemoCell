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

//#include "immersedWallParticle3D.h"
//#include "immersedWallParticle3D.hh"
//#include "immersedWallParticleFunctional3D.h"
//#include "immersedWallParticleFunctional3D.hh"
//#include "immersedWallParticleVtk3D.h"
//#include "immersedWallParticleVtk3D.hh"
//#include "shellModel3D.h"
//#include "shellModel3D.hh"

#include "ficsionInit.hh"
#include "immersedCells3D.h"
#include "immersedCells3D.hh"
#include "immersedCellsFunctional3D.h"
#include "immersedCellsFunctional3D.hh"


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

//    ImageWriter<T> imageWriter("leeloo");
//    add(force, forceScalar, force.getBoundingBox());
//    imageWriter.writeScaledPpm(scalarField, createFileName("PPM", iter, 6));

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

    T u = 0.01;
    T Re = 100;
    plint N = 20;
    T lx = 5.;
    T ly = 0.5;
    T lz = 0.5;
    string paramXmlFileName;
    global::argv(1).read(paramXmlFileName);

	XMLreader document(paramXmlFileName);
	pcout << "reading.." <<std::endl;
	document["cell"]["shellDensity"].read(shellDensity);
	document["cell"]["k_rest"].read(k_rest);
	document["cell"]["k_stretch"].read(k_stretch);
	document["cell"]["k_shear"].read(k_shear);
	document["cell"]["k_bend"].read(k_bend);

	document["parameters"]["u"].read(u);
	document["parameters"]["Re"].read(Re);
	document["parameters"]["N"].read(N);
	document["parameters"]["lx"].read(lx);
	document["parameters"]["ly"].read(ly);
	document["parameters"]["lz"].read(lz);
	plint forceToFluid = 0, shape = 0;
	T radius = 0.0;
	document["ibm"]["forceToFluid"].read(forceToFluid);
	document["ibm"]["shape"].read(shape);
	document["ibm"]["radius"].read(radius);


    IncomprFlowParam<T> parameters(
            u, // u
            Re, // Re
            N,   // N
            lx,        // lx
            ly,        // ly
            lz         // lz
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
    plint x0 = 10; plint x1 = nx-x0;

    Box3D inlet(0,  x0,  0, ny-1, 0, nz-1);
    Box3D outlet(x1, nx-1,0, ny-1, 0, nz-1);

    std::vector<plint> pos;
//    pos.push_back(ny/2); //pos.push_back(30);pos.push_back(50);
    pos.push_back(ny/3);pos.push_back((2*ny)/3);
    
    std::vector<Array<T,3> > centers;
    std::vector<plint > radii;
    plint nMax = 1;
    for (plint iN = 0; iN < nMax; ++iN) { // create 40 * inlet amount of particles
        for (pluint iA = 0; iA < pos.size(); ++iA) {
            for (pluint iB = 0; iB < pos.size(); ++iB) {
                centers.push_back(Array<T,3>(5+iN,pos[iA],pos[iB]));
                radii.push_back(radius);
            }
        }
    }

    plint numOfBloodCellsPerInlet = pos.size() * pos.size(); // number used for the generation of bloodCells at inlet

    std::vector<plint> tags;
    plint numPartsPerBloodCell = 0; plint slice = 0; // number of particles per tag and number of slice of created particles
    TriangleBoundary3D<T> bloodCells = createCompleteMesh(centers, radii, tags, numPartsPerBloodCell, shape);
	generateBloodCells(immersedParticles, inlet, tags, bloodCells, numPartsPerBloodCell, numOfBloodCellsPerInlet, slice);

    std::vector<plint> numParts(tags.size());
    for (pluint iA = 0; iA < tags.size(); ++iA) {
        numParts[iA] = countParticles(immersedParticles, immersedParticles.getBoundingBox(), tags[iA]);
    }
    
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&immersedParticles);

    std::vector<MultiBlock3D*> particleLatticeArg;
    particleLatticeArg.push_back(&immersedParticles);
    particleLatticeArg.push_back(&lattice);

    CellModel3D<T> cellModel(shellDensity, k_rest, k_stretch, k_shear, k_bend, k_shear, k_shear);

    plint maxIter = 100000;
    plint imageIter = 10;
    pcout << std::endl << "Starting simulation" << std::endl;
    global::timer("sim").start();
    applyProcessingFunctional ( // copy fluid velocity on particles
        new FluidVelocityToImmersedWall3D<T,DESCRIPTOR>(),
        immersedParticles.getBoundingBox(), particleLatticeArg);
    plint itotParticles = countParticles(immersedParticles, immersedParticles.getBoundingBox());

    for (plint i=0; i<maxIter; ++i) {

    	applyProcessingFunctional ( // compute force applied on the particles by springs
            new ComputeImmersedElasticForce3D<T,DESCRIPTOR> (
                bloodCells, cellModel.clone() ), // used because pushSelect is not used
            immersedParticles.getBoundingBox(), particleArg );
        if (forceToFluid != 0) {
			setExternalVector( lattice, lattice.getBoundingBox(),
							   DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T,DESCRIPTOR<T>::d>(0.0,0.0,0.0));
			applyProcessingFunctional ( // compute force applied on the fluid by the particles
				new ForceToFluid3D<T,DESCRIPTOR> (),
					immersedParticles.getBoundingBox(), particleLatticeArg );
        }
        lattice.collideAndStream();

        applyProcessingFunctional ( // copy fluid velocity on particles
            new FluidVelocityToImmersedWall3D<T,DESCRIPTOR>(),
            immersedParticles.getBoundingBox(), particleLatticeArg);

        applyProcessingFunctional ( // advance particles in time according to a velocity, acceleration, ...
            new AdvanceParticlesFunctional3D<T,DESCRIPTOR>,
            immersedParticles.getBoundingBox(), particleArg );

        bloodCells.pushSelect(0,1);
        applyProcessingFunctional ( // update mesh position
            new CopyParticleToVertex3D<T,DESCRIPTOR>(bloodCells.getMesh()),
            immersedParticles.getBoundingBox(), particleArg);
        bloodCells.popSelect();


        
//        deleteBloodCell(immersedParticles, outlet, numParts, tags, bloodCells, centers, radii );
//        if (slice < 1) {
//            bool created = generateBloodCells(immersedParticles, inlet, tags, bloodCells, numPartsPerBloodCell, numOfBloodCellsPerInlet, slice );
//            pcout << "Used \n";
//        }

        if (i%imageIter==0) {
        	plint totParticles = countParticles(immersedParticles, immersedParticles.getBoundingBox());
            pcout << i << " totParticles = " << totParticles << std::endl;
            PLB_ASSERT(itotParticles == totParticles);
            pcout << "Write Particle VTK. " ;
            std::vector<std::string> force_scalarNames;
            force_scalarNames.push_back("pressure");
            force_scalarNames.push_back("wss");
            std::vector<std::string> velocity_scalarNames;
            std::vector<std::string> force_vectorNames;
            force_vectorNames.push_back("force");
            std::vector<std::string> velocity_vectorNames;
            velocity_vectorNames.push_back("velocity");
			bool dynamicMesh = true;
			plint tag = -1; // Take all triangles.
//			bloodCells.writeAsciiSTL("tmp.stl");
			// serialize the particle information to write them.
			// a correspondance between the mesh and the particles is made.
			writeImmersedSurfaceVTK (
					bloodCells,
					*getParticlePosAndVelocity(immersedParticles),
					velocity_scalarNames, velocity_vectorNames,
					global::directories().getOutputDir()+createFileName("RBC",i,6)+".vtk", dynamicMesh, tag );

            writeVTK(lattice, parameters, i);
            pcout << "Timer: " << global::timer("sim").stop() << endl;
            global::timer("sim").restart();
        }
    }

}
