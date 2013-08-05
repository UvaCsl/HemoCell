/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011 FlowKit Sarl
 * Avenue de Chailly 23
 * 1012 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef IMMERSED_PARTICLE_VTK_3D_HH
#define IMMERSED_PARTICLE_VTK_3D_HH

#include "palabos3D.h"
#include "palabos3D.hh"

#include "immersedCellParticleVtk3D.h"
#include "immersedCellParticle3D.h"
#include "immersedCellParticle3D.hh"

namespace plb {

template<typename T, template<typename U> class Descriptor>
void writeImmersedSurfaceVTK( TriangleBoundary3D<T> const& boundary,
                      MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particles,
                      std::vector<std::string> const& scalars,
                      std::vector<std::string> const& vectors,
                      std::string const& fName, bool dynamicMesh, plint tag )
{
    writeImmersedSurfaceVTK( boundary, particles, scalars, vectors, fName, dynamicMesh, tag,
                     std::vector<T>(), std::vector<T>() );
}

template<typename T, template<typename U> class Descriptor>
void writeImmersedSurfaceVTK( TriangleBoundary3D<T> const& boundary,
                      MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particles,
                      std::vector<std::string> const& scalars,
                      std::vector<std::string> const& vectors,
                      std::string const& fName, bool dynamicMesh, plint tag,
                      std::vector<T> const& scalarFactor, std::vector<T> const& vectorFactor )
{
    SparseBlockStructure3D blockStructure(particles.getBoundingBox());
    blockStructure.addBlock(particles.getBoundingBox(), 0);
    plint envelopeWidth=1;
    MultiBlockManagement3D serialMultiBlockManagement (
            blockStructure, new OneToOneThreadAttribution, envelopeWidth );

    MultiParticleField3D<DenseParticleField3D<T,Descriptor> > multiSerialParticles (
            serialMultiBlockManagement,
            defaultMultiBlockPolicy3D().getCombinedStatistics() );

    copy ( particles, particles.getBoundingBox(), multiSerialParticles, particles.getBoundingBox() );
    if (global::mpi().isMainProcessor()) {
        ParticleField3D<T,Descriptor>& atomicSerialParticles =
            dynamic_cast<ParticleField3D<T,Descriptor>&>(multiSerialParticles.getComponent(0));

        std::vector<Particle3D<T,Descriptor>*> found;
        SmartBulk3D oneBlockBulk(serialMultiBlockManagement, 0);
        atomicSerialParticles.findParticles(oneBlockBulk.toLocal(particles.getBoundingBox()), found);
        if (found.size() > 0) {
            vtkForImmersedVertices(found, boundary, scalars, vectors, fName, dynamicMesh, tag, scalarFactor, vectorFactor);
        }
        else {
            pcout << "No particles found inside the domain" << std::endl;
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void vtkForImmersedVertices(std::vector<Particle3D<T,Descriptor>*> const& particles,
                    TriangleBoundary3D<T> const& boundary,
                    std::vector<std::string> const& scalars,
                    std::vector<std::string> const& vectors,
                    std::string fName, bool dynamicMesh, plint tag,
                    std::vector<T> const& scalarFactor, std::vector<T> const& vectorFactor )
{
    PLB_ASSERT( scalarFactor.empty() || scalarFactor.size()==scalars.size() );
    PLB_ASSERT( vectorFactor.empty() || vectorFactor.size()==vectors.size() );
    TriangularSurfaceMesh<T> const& mesh = boundary.getMesh();
    // If this assertion fails, a likely explanation is that the margin of your sparse block
    // structure is too small, and one of the particles was outside the allocated domain.
//     PLB_PRECONDITION((plint)particles.size() == mesh.getNumVertices());

    ImmersedCellParticle3D<T,Descriptor>* p0 =
        dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (particles[0]);
    std::vector<std::string> vectorNames, scalarNames;
    for (plint iVector = 0; iVector < p0->getVectorsNumber(); ++iVector) {
        vectorNames.push_back(p0->getVectorName(iVector));
    }
    for (plint iScalar = 0; iScalar < p0->getScalarsNumber(); ++iScalar) {
    	scalarNames.push_back(p0->getScalarName(iScalar));
    }

    std::ofstream ofile(fName.c_str());
    ofile << "# vtk DataFile Version 3.0\n";
    ofile << "Surface mesh created with Palabos\n";
    ofile << "ASCII\n";
    ofile << "DATASET UNSTRUCTURED_GRID\n";

    ofile << "POINTS " << particles.size()
          << (sizeof(T)==sizeof(double) ? " double" : " float")
          << "\n";

    std::vector<std::vector<T> > scalarData(scalarNames.size());
    for (pluint iScalar=0; iScalar<scalarNames.size(); ++iScalar) {
        scalarData[iScalar].resize(particles.size());
    }
    std::vector<std::vector<Array<T,3> > > vectorData(vectorNames.size());
    for (pluint iVector=0; iVector<vectorNames.size(); ++iVector) {
        vectorData[iVector].resize(particles.size());
    }
    std::vector<Array<T,3> > posVect(particles.size());
    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
        plint iVertex = particles[iParticle]->getTag();
        posVect[iVertex] = particles[iParticle]->getPosition();
        posVect[iVertex] *= boundary.getDx();
        //posVect[iVertex] += boundary.getPhysicalLocation();
        for (pluint iScalar=0; iScalar<scalarNames.size(); ++iScalar) {
            T scalar;
            ImmersedCellParticle3D<T,Descriptor>* iparticle =
                dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (particles[iParticle]);
            iparticle->getScalar(iScalar, scalar);
            if (!scalarFactor.empty()) {
                scalar *= scalarFactor[iScalar];
            }
            scalarData[iScalar][iVertex] = scalar;
        }
        for (pluint iVector=0; iVector<vectorNames.size(); ++iVector) {
            Array<T,3> vector;
            ImmersedCellParticle3D<T,Descriptor>* iparticle =
                dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (particles[iParticle]);
            iparticle->getVector(iVector, vector);
            if (!vectorFactor.empty()) {
                vector *= vectorFactor[iVector];
            }
            vectorData[iVector][iVertex] = vector;
        }
    }

    for (pluint iVertex=0; iVertex<particles.size(); ++iVertex) {
        ofile << posVect[iVertex][0] << " "
              << posVect[iVertex][1] << " "
              << posVect[iVertex][2] << "\n";
    }
    ofile << "\n";

    plint numWallTriangles=0;
    for (plint iTriangle=0; iTriangle<mesh.getNumTriangles(); ++iTriangle) {
        if ( tag<0 || boundary.getTag(iTriangle)==tag )
        {
            ++numWallTriangles;
        }
    }

    ofile << "CELLS " << numWallTriangles
          << " " << 4*numWallTriangles << "\n";

    for (plint iTriangle=0; iTriangle<mesh.getNumTriangles(); ++iTriangle) {
        plint i0 = mesh.getVertexId(iTriangle, 0);
        plint i1 = mesh.getVertexId(iTriangle, 1);
        plint i2 = mesh.getVertexId(iTriangle, 2);
        if ( tag<0 || boundary.getTag(iTriangle)==tag )
        {
            ofile << "3 " << i0 << " " << i1 << " " << i2 << "\n";
        }
    }
    ofile << "\n";

    ofile << "CELL_TYPES " << numWallTriangles << "\n";
    for (plint i=0; i<numWallTriangles; ++i) {
        ofile << "5\n";
    }
    ofile << "\n";

    ofile << "POINT_DATA " << particles.size() << "\n";
    for (pluint iVector=0; iVector<vectorNames.size(); ++iVector) {
        ofile << "VECTORS " << vectorNames[iVector]
              << (sizeof(T)==sizeof(double) ? " double" : " float")
              << "\n";
        for (plint iVertex=0; iVertex<(plint)particles.size(); ++iVertex) {
            ofile << vectorData[iVector][iVertex][0] << " "
                  << vectorData[iVector][iVertex][1] << " "
                  << vectorData[iVector][iVertex][2] << "\n";
        }
        ofile << "\n";
    }

    for (pluint iScalar=0; iScalar<scalarNames.size(); ++iScalar) {
        ofile << "SCALARS " << scalarNames[iScalar]
              << (sizeof(T)==sizeof(double) ? " double" : " float")
              << " 1\n"
              << "LOOKUP_TABLE default\n";
        for (plint iVertex=0; iVertex<(plint)particles.size(); ++iVertex) {
            ofile << scalarData[iScalar][iVertex] << "\n";
        }
        ofile << "\n";
    }
}


template<typename T>
void writeImmersedPointsVTK(TriangleBoundary3D<T> const& boundary, std::vector<plint> const& vertices, T const& dx,
                    std::string fName)
{
    TriangularSurfaceMesh<T> const& mesh = boundary.getMesh();
    std::ofstream ofile(fName.c_str());
    ofile << "# vtk DataFile Version 3.0\n";
    ofile << "Surface point created with Palabos and ficsion\n";
    ofile << "ASCII\n";
    ofile << "DATASET POLYDATA\n";

    ofile << "POINTS " << vertices.size()
          << (sizeof(T)==sizeof(double) ? " double" : " float")
          << "\n";
    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T,3> vertex = mesh.getVertex(vertices[i]);
        vertex *= dx;
        ofile << vertex[0] << " " <<vertex[1] << " " <<vertex[2] << "\n";
    }
    ofile << "POINT_DATA " << vertices.size() << "\n";
    ofile << "SCALARS ID float 1 \n";
    ofile << "LOOKUP_TABLE default \n";
    for (pluint i = 0; i < vertices.size(); ++i) {
        ofile << i*1.0 << "\n";
    }
    ofile.close();
}


}  // namespace plb

#endif  // IMMERSED_PARTICLE_VTK_3D_HH
