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
#include "readPositionsBloodCells.h"
#include "tools/packCells/geometry.h"
#include <sstream>
#include "logfile.h"
#include "hemocell.h"
#include "preInlet.h"

namespace hemo {
  
inline void meshRotation (TriangularSurfaceMesh<T> * mesh, hemo::Array<T,3> rotationAngles) {
    plb::Array<T,2> xRange, yRange, zRange;
    mesh->computeBoundingBox (xRange, yRange, zRange);
    plb::Array<T,3> meshCenter = plb::Array<T,3>(xRange[1] + xRange[0], yRange[1] + yRange[0], zRange[1] + zRange[0]) * (T)0.5;
    mesh->translate((T)-1.0 * meshCenter);
    mesh->rotateXYZ(rotationAngles[0], rotationAngles[1], rotationAngles[2]);
    mesh->translate(meshCenter);
}

inline void positionCellInParticleField(HEMOCELL_PARTICLE_FIELD& particleField, BlockLattice3D<T,DESCRIPTOR>& fluid,
                                            TriangularSurfaceMesh<T> * mesh, hemo::Array<T,3> startingPoint, plint cellId, pluint celltype) {
    plint nVertices=mesh->getNumVertices();
    Box3D fluidbb = fluid.getBoundingBox();
    Dot3D rrelfluidloc; // used later;
    int denyLayerSize; // used later;
    HemoCellParticle to_add_particle; // used_later

    for (plint iVertex=0; iVertex < nVertices; ++iVertex) {
        hemo::Array<T,3> vertex = startingPoint + mesh->getVertex(iVertex);
        
        //If we cannot place it in the particle field continue
        if (!particleField.isContainedABS(vertex,particleField.getBoundingBox())) { continue; }
        
        //Custom contains function for generic blocklattice3D structure
        Dot3D absfluidloc(int(vertex[0]+0.5),int(vertex[1]+0.5),int(vertex[2]+0.5));
        Dot3D relfluidloc = absfluidloc - fluid.getLocation();
        
        
        if (!(relfluidloc.x < fluidbb.x0 || relfluidloc.x > fluidbb.x1 ||
            relfluidloc.y < fluidbb.y0 || relfluidloc.y > fluidbb.y1 ||
            relfluidloc.z < fluidbb.z0 || relfluidloc.z > fluidbb.z1)) {
              // Test when particle is inside if in a boundary -> dont add this particle
            if (fluid.get(relfluidloc.x,relfluidloc.y,relfluidloc.z).getDynamics().isBoundary()) {
        	goto no_add;
            }
        } 
            
      //bool neighboringBoundariesAnywhere = false;  
      // Deny particles that are in the outer most layer, aka. the "shear layer"    
      denyLayerSize = ((*particleField.cellFields)[celltype]->minimumDistanceFromSolid*1e-6)/param::dx;
      for (int px = -denyLayerSize; px <= denyLayerSize; ++px) {  
        for (int py = -denyLayerSize; py <= denyLayerSize; ++py) { 
          for (int pz = -denyLayerSize; pz <= denyLayerSize; ++pz) {
            rrelfluidloc = {relfluidloc.x+px,relfluidloc.y+py,relfluidloc.z+pz};
            if (!(rrelfluidloc.x < fluidbb.x0 || rrelfluidloc.x > fluidbb.x1 ||
                  rrelfluidloc.y < fluidbb.y0 || rrelfluidloc.y > fluidbb.y1 ||
                  rrelfluidloc.z < fluidbb.z0 || rrelfluidloc.z > fluidbb.z1)) {
              if (fluid.get(rrelfluidloc.x,rrelfluidloc.y,rrelfluidloc.z).getDynamics().isBoundary()) {
        	goto no_add;
              }
            }
          }  
        }  
      }
      
      to_add_particle = HemoCellParticle(vertex,cellId,iVertex,celltype);
      particleField.addParticle(&to_add_particle);
no_add:;
    }
}

int getTotalNumberOfCells(HemoCellFields & cellFields){
  int nCells, totalCells = 0;
  for (pluint j = 0; j < cellFields.size(); j++) {
    fstream fIn;
    fIn.open(cellFields[j]->name + ".pos", fstream::in);
    fIn >> nCells;
    totalCells += nCells;
    fIn.close();
  }
  return totalCells;
  
}

bool first_time =true;
void getReadPositionsBloodCellsVector(Box3D realDomain,
                                            std::vector<TriangularSurfaceMesh<T>* > & meshes,
                                            std::vector<plint> & Np,
                                            std::vector<std::vector<hemo::Array<T,3> > > & positions,
                                            std::vector<std::vector<plint> > & cellIds,
                                            std::vector<std::vector<hemo::Array<T,3> > > & randomAngles,
                                            T dx, Config & cfg, HemoCellFields & cellFields,
                                            HemoCellParticleField & particleField)
{

    


    vector<vector3> packPositions[cellFields.size()];
    vector<vector3> packAngles[cellFields.size()];
    vector<plint> cellIdss[cellFields.size()];

    Np.resize(cellFields.size());

    int cellid = 0;
    // TODO: proper try-catch
    for(pluint j = 0; j < Np.size(); j++) {

        // Reading data from file
        fstream fIn;
        fIn.open(cellFields[j]->name + ".pos", fstream::in);

        if(!fIn.is_open())
        {
            cout << "*** WARNING! particle positions input file " << cellFields[j]->name << ".pos does not exist!" << endl;
        }

        fIn >> Np[j];
        if (first_time) {
          hlog << "(readPositionsBloodCels) Particle count in file (" << cellFields[j]->name << "): " << Np[j] << "." << endl;
        }
        packPositions[j].resize(Np[j]); packAngles[j].resize(Np[j]);cellIdss[j].resize(Np[j]);
        int less = 0;
        for (plint i = 0; i < Np[j]; i++) {
            fIn >> 
              packPositions[j][i-less][0] >> packPositions[j][i-less][1] >> packPositions[j][i-less][2] >> 
              packAngles[j][i-less][0] >> packAngles[j][i-less][1] >> packAngles[j][i-less][2];
            packAngles[j][i-less] *= PI/180.0; // Deg to Rad
            packAngles[j][i-less] *= -1.0;  // Right- to left-handed coordinate system
            cellIdss[j][i-less] = cellid;

            if (cellFields.hemocell.preInlet->initialized) {
              //Translate system to preInlet Location (set 0,0,0 point)
              packPositions[j][i-less][0] += cellFields.hemocell.preInlet->location.x0/dx;
              packPositions[j][i-less][1] += cellFields.hemocell.preInlet->location.y0/dx;
              packPositions[j][i-less][2] += cellFields.hemocell.preInlet->location.z0/dx;
            }
            
            
            //Check if it actually fits (mostly) in this atomic block
            if (packPositions[j][i-less][0]*dx < realDomain.x0 - cfg["domain"]["particleEnvelope"].read<int>() ||
                packPositions[j][i-less][0]*dx > realDomain.x1 + cfg["domain"]["particleEnvelope"].read<int>() ||
                packPositions[j][i-less][1]*dx < realDomain.y0 - cfg["domain"]["particleEnvelope"].read<int>() ||
                packPositions[j][i-less][1]*dx > realDomain.y1 + cfg["domain"]["particleEnvelope"].read<int>() ||
                packPositions[j][i-less][2]*dx < realDomain.z0 - cfg["domain"]["particleEnvelope"].read<int>() ||
                packPositions[j][i-less][2]*dx > realDomain.z1 + cfg["domain"]["particleEnvelope"].read<int>()) {
              less ++;
            }
            cellid++;

        }
        //Destroy unwanted particles and adjust list size;
        Np[j] -= less;
        packPositions[j].resize(Np[j]);
        packAngles[j].resize(Np[j]);
        cellIdss[j].resize(Np[j]);

        fIn.close();
    }
    
    cellFields.number_of_cells = getTotalNumberOfCells(cellFields);
    
    // Copy results of packing to appropriate hemo::Arrays for ficsion
    positions.clear();	positions.resize(Np.size());
    randomAngles.clear(); randomAngles.resize(Np.size());
    cellIds.clear();	cellIds.resize(Np.size());


    for (pluint i = 0; i < Np.size(); ++i)
    {
        positions[i].resize(Np[i]);
        cellIds[i].resize(Np[i]);
        randomAngles[i].resize(Np[i]);

        for (plint j = 0; j < Np[i]; ++j)
        {
            // Store mesh positions and rotations
            //randomAngles[i][j] = hemo::Array<T, 3>(1.0,0.0,0.0);
            randomAngles[i][j] = hemo::Array<T, 3>({packAngles[i][j][0], packAngles[i][j][1], packAngles[i][j][2]});
            positions[i][j] =hemo::Array<T,3>({packPositions[i][j][0], packPositions[i][j][1], packPositions[i][j][2]});
            cellIds[i][j] = cellIdss[i][j];

        }
    }
    first_time = false;
}



/* ******** ReadPositionMultipleCellField3D *********************************** */
void ReadPositionsBloodCellField3D::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{  
    int numberOfCellFields = blocks.size() -1;
    //T ratio;
    BlockLattice3D<T,DESCRIPTOR>& fluid =
            *dynamic_cast<BlockLattice3D<T,DESCRIPTOR>*>(blocks[0]);
    std::vector<HEMOCELL_PARTICLE_FIELD* > particleFields(numberOfCellFields);
    std::vector<T> volumes(numberOfCellFields);
    std::vector<TriangularSurfaceMesh<T>* > meshes(numberOfCellFields);
    std::vector<ElementsOfTriangularSurfaceMesh<T> > emptyEoTSM(numberOfCellFields);
    std::vector<plint> Np(numberOfCellFields);

    T totalVolumeFraction=0;

    //T Vdomain = domain.getNx() * domain.getNy() * domain.getNz();

    Dot3D fLocation(fluid.getLocation());
    Box3D realDomain(
            domain.x0 + fLocation.x, domain.x1 + fLocation.x,
            domain.y0 + fLocation.y, domain.y1 + fLocation.y,
            domain.z0 + fLocation.z, domain.z1 + fLocation.z );
    plb::Array<T,2> xRange, yRange, zRange;

    for (pluint iCF = 0; iCF < cellFields.size(); ++iCF) {

        TriangularSurfaceMesh<T> * mesh = copyTriangularSurfaceMesh(cellFields[iCF]->getMesh(), emptyEoTSM[iCF]);
        mesh->computeBoundingBox (xRange, yRange, zRange);
        mesh->translate(plb::Array<T,3>(-(xRange[0]+xRange[1])/2.0, -(yRange[0]+yRange[1])/2.0, -(zRange[0]+zRange[1])/2.0));
        meshes[iCF] = mesh;
        volumes[iCF] = MeshMetrics<T>(*mesh).getVolume();

        totalVolumeFraction += volumes[iCF];
        particleFields[iCF] = ( dynamic_cast<HEMOCELL_PARTICLE_FIELD*>(blocks[iCF+1]) );
        particleFields[iCF]->removeParticles(particleFields[iCF]->getBoundingBox());
    }
   

    std::vector<std::vector<hemo::Array<T,3> > > positions;
    std::vector<std::vector<plint> > cellIds;
    std::vector<std::vector<hemo::Array<T,3> > > randomAngles;

    // Note: this method uses the center of the particles for location
    T posRatio = 1e-6/dx;
    getReadPositionsBloodCellsVector(realDomain, meshes, Np, positions, cellIds, randomAngles, posRatio, cfg, cellFields,*particleFields[0]);

    // Change positions to match dx (it is in um originally)
    T wallWidth = 0; // BB wall in [lu]. Offset to count in width of the wall in particle position (useful for pipeflow, not necessarily useful elswhere)
    
    for (pluint iCF = 0; iCF < positions.size(); ++iCF)
    {
        for (pluint c = 0; c < positions[iCF].size(); ++c)
        {
            ElementsOfTriangularSurfaceMesh<T> emptyEoTSMCopy;
    	    TriangularSurfaceMesh<T> * meshCopy = copyTriangularSurfaceMesh(*meshes[iCF], emptyEoTSMCopy);
    	    
            meshRotation (meshCopy, randomAngles[iCF][c]);
            
            positionCellInParticleField(*(particleFields[iCF]), fluid,
                                         meshCopy, positions[iCF][c]*posRatio+wallWidth, cellIds[iCF][c], iCF); 
			delete meshCopy;
        }

        
        //SYNC THEM ENVELOPES
        //cellFields.syncEnvelopes();
        // DELETE CELLS THAT ARE NOT WHOLE
        particleFields[iCF]->deleteIncompleteCells(iCF,false);
        std::vector<HemoCellParticle*> particles;
        particleFields[iCF]->findParticles(particleFields[iCF]->getBoundingBox(), particles, iCF);
        
//delete meshes[iCF];
    }
    //cout << "Atomic Block ID: " << particleFields[0]->atomicBlockId;
    //cout    << " Total complete cells (with periodicity): " << particleFields[0]->get_lpc().size() << std::endl;
}


ReadPositionsBloodCellField3D* ReadPositionsBloodCellField3D::clone() const {
    return new ReadPositionsBloodCellField3D(*this);
}

void ReadPositionsBloodCellField3D::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // Fluid field.
    for (unsigned int iField = 1; iField < modified.size(); ++iField) {
        modified[iField] = modif::nothing; // Particle fields.
    }
}

void ReadPositionsBloodCellField3D::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = false; // Fluid field.
    for (unsigned int iField = 1; iField < isWritten.size(); ++iField) {
        isWritten[iField] = false; // Particle fields.
    }

}


BlockDomain::DomainT ReadPositionsBloodCellField3D::appliesTo() const {
    return BlockDomain::bulk;
}



void readPositionsBloodCellField3D(HemoCellFields & cellFields, T dx, Config & cfg) {
    std::vector<MultiBlock3D *> fluidAndParticleFieldsArg;

    fluidAndParticleFieldsArg.push_back(cellFields.lattice);

    for (pluint icf = 0; icf < cellFields.size(); ++icf) {
        fluidAndParticleFieldsArg.push_back(cellFields[icf]->getParticleField3D());
    }
    hlog << "(readPositionsBloodCels) Reading particle positions..." << std::endl;
    
    if (cellFields.hemocell.preInlet->initialized && !cellFields.hemocell.partOfpreInlet) { } else {
    applyProcessingFunctional(
            new ReadPositionsBloodCellField3D(cellFields, dx, cfg),
            cellFields.lattice->getBoundingBox(), fluidAndParticleFieldsArg);
    hlogfile << "Mpi Process: " << global::mpi().getRank()  << " Completed loading particles" << std::endl;
    }
    cellFields.number_of_cells = getTotalNumberOfCells(cellFields);
}

}