#ifndef SHAPE_MEMORY_MODEL_3D_H
#define SHAPE_MEMORY_MODEL_3D_H

#include "cellModel3D.h"
#include "config.h"


#ifndef KBT__
#define KBT__
const double kB_p = 1.3806503e-23; // In SI, m2 kg s-2 K-1 (or J/K)
const double kBT_p = 4.100531391e-21; // In SI, m2 kg s-2 (or J) for T=300
double kB=0, kBT=0;
double dNewton=0;
#endif  // KBT__

#ifndef PI__
#define PI__
const double pi = 4.*atan(1.);
#endif  // PI__

namespace plb {


template<typename T>
void getCellShapeQuantitiesFromMesh(TriangularSurfaceMesh<T> const& dynMesh,
                            vector<T> & eqAreaPerTriangle, map<plint,T> & eqLengthPerEdge, map<plint,T> & eqAnglePerEdge,
                            plint cellNumTriangles, plint cellNumPartsPerCell);



class ShapeMemoryModel3D :  public ShellModel3D<double>
{
public:
    static ShapeMemoryModel3D* PlateletShapeMemoryModel3D(Config* cfg, 
                       double dx_, double dt_, double dm_,
                       TriangularSurfaceMesh<double> const& meshElement);
    /* All input should be in dimensionless units */
    static ShapeMemoryModel3D* RBCShapeMemoryModel3D(Config* cfg, 
                       double dx_, double dt_, double dm_,
                       TriangularSurfaceMesh<double> const& meshElement);
private:
    ShapeMemoryModel3D(Config* config,
                       TriangularSurfaceMesh<double> const& meshElement);
    void Initialize(TriangularSurfaceMesh<double> const& meshElement);

public:
    //ShapeMemoryModel3D(ShapeMemoryModel3D const& rhs); //Let the compiler auto
    //generate the copy function

    ~ShapeMemoryModel3D() { } ;

    Array<double,3> computeElasticForce (TriangleBoundary3D<double> const& boundary, plint iVertex )
    { 
      return Array<double,3>(0.0, 0.0, 0.0); 
    }

    virtual void computeCellForce (Cell3D<double,DESCRIPTOR> * cell);
    
    // TODO: these are assumptions made on normal RBC behaviour, might not be OK for every cell type!
    // TODO: these should be moved to the user end
    virtual plint getMaximumEdgeExtensionLengthLU() { return maximum(ceil(maxLength + 0.5),4); };
    virtual plint getMaxCellDiameterLU() { return maximum(ceil(eqLengthRatio*cellRadiusLU),4); };

    virtual double getDx() { return dx; };
    virtual double getDt() { return dt; };
    virtual double getDm() { return dm; };

    virtual SyncRequirements & getSyncRequirements() {return syncRequirements;} ;
    virtual SyncRequirements const& getSyncRequirements() const {return syncRequirements;} ;

private:
    inline void computeCellForceHighOrder (Cell3D<double,DESCRIPTOR> *cell);
    inline void computeCellForceSuresh (Cell3D<double,DESCRIPTOR> * cell);

    plint getTriangleId(plint iTriangle);
    plint getEdgeId(plint iVertex, plint jVertex);
private:
    MeshMetrics<double> meshmetric;
    SyncRequirements syncRequirements;
    double maxLength, cellRadiusLU;
    double k_rest, k_shear, k_bend, k_stretch, k_inPlane, k_elastic, k_surface, k_volume, k_WLC;
    double C_elastic;
    double eta_m, gamma_T, gamma_C;
    double eqLength, eqArea, eqAngle;
    vector<double> eqAreaPerTriangle;
    map<plint,double> eqLengthPerEdge, eqAnglePerEdge;
    double eqVolume, eqSurface, eqTileSpan;
    double persistenceLengthCoarse, eqLengthRatio;
    double dx, dt, dm;
    double persistenceLengthFine;
    pluint cellNumTriangles, cellNumVertices;
    pluint materialModel;
public:
    /* Computes the equilibrium quantities to correspond to the an inflated cell with
     * 		eqVolume=ratio*eqVolume.
     * Can also be used for deflation. */
    virtual void inflate(double ratio) {
        eqVolume = meshmetric.getVolume() * ratio*ratio*ratio;
        eqSurface = meshmetric.getSurface() * ratio*ratio;

        getCellShapeQuantitiesFromMesh(meshmetric.getMesh(), eqAreaPerTriangle, eqLengthPerEdge, eqAnglePerEdge, cellNumTriangles, cellNumVertices);
        persistenceLengthCoarse = persistenceLengthFine/dx * sqrt( (cellNumVertices-2.0) / (23867-2.0)) * ratio;

    	for (pluint i=0; i < eqAreaPerTriangle.size() ; i++) {
    		eqAreaPerTriangle[i] *= ratio*ratio;
    	}
    	typename std::map<plint,double>::iterator iter;
        for (iter = eqLengthPerEdge.begin(); iter != eqLengthPerEdge.end(); ++iter) {
        	eqLengthPerEdge[iter->first] *= ratio;
        }
    	eqLength *= ratio;
    }
public:
    /* Coefficients */
    /* Coefficients */
    virtual double& getRestingStiffness() { return k_rest; }
    virtual double& getBendingStiffness() { return k_bend; }
    virtual double& getStretchingStiffness() { return k_stretch; }
    virtual void setRestingStiffness(double value) { k_rest = value; }
    virtual void setBendingStiffness(double value) { k_bend = value; }
    virtual void setStretchingStiffness(double value) { k_stretch = value; }

    /* TODO: Fix change of coefficients */

    // Units are N/m
    virtual double getMembraneShearModulus() {
        double Lmax = eqLength*eqLengthRatio;
        double x0 = 1.0/eqLengthRatio;
        // k_inPlane = k_WLC_ * kBT /(4.0*persistenceLengthCoarse);
        double kP =  k_inPlane;
        double k_rep = -1.0*((eqLength*eqLength)* kP * x0 * (-6 + (9 - 4*x0)*x0)/( (x0-1)*(x0-1) ));
        double oneMinusX0 = 1-x0;
        double oneMinusX0Square = oneMinusX0*oneMinusX0;
        double oneMinusX0Cube = oneMinusX0*oneMinusX0*oneMinusX0;
        // Fedosov et al.  doi:10.1016/j.bpj.2010.02.002
        // A Multiscale Red Blood Cell Model with Accurate Mechanics, Rheology, and Dynamics  
        double ansFedosov = kP * sqrt(3) / (Lmax*x0) * (x0/(2.0*oneMinusX0Cube) - 1.0/(4.0*oneMinusX0Square) + 1.0/4.0) + 
                3 * sqrt(3) * k_rep/(4*eqLength*eqLength*eqLength);
        // Reasor et al. doi:10.1002/fld.2534
        // Coupling the lattice-Boltzmann and spectrin-link methods for the direct numerical simulation of cellular blood flow
        double ansReasor = sqrt(3) * kP / (Lmax*x0)*( 3.0/( 4*(1-x0)*(1-x0) ) - 3.0/4.0 +4*x0 + x0/(2.0 * (1-x0)*(1-x0)*(1-x0)));
        return ansFedosov;
    }
    // Units are N/m
    virtual double getMembraneElasticAreaCompressionModulus() {
        double Lmax = eqLength*eqLengthRatio;
        double x0 = 1.0/eqLengthRatio;
        // T kP =  kBT /(4.0*persistenceLengthCoarse);
        double kP =  k_inPlane;
        return sqrt(3) * kP / (Lmax*(1-x0)*(1-x0))*(1.5*(6-9*x0+4*x0*x0) + (1+2*(1-x0)*(1-x0)*(1-x0))/(1-x0)  );
    }
    // Units are N/m
    virtual double getYoungsModulus() {
        double mu0 = getMembraneShearModulus();
        double K = getMembraneElasticAreaCompressionModulus();
        return (4*K*mu0)/(K+mu0);
    }
    // Dimensionless number
    virtual double getPoissonRatio() {
        double mu0 = getMembraneShearModulus();
        double K = getMembraneElasticAreaCompressionModulus();
        return (K-mu0)/(K+mu0);
    }
    // Units are N s/m
    virtual double& getMembraneShearViscosity() { return eta_m; }
    virtual void setMembraneShearViscosity(double value) {
        eta_m = value;
        gamma_T = (eta_m * 12.0/(13.0 * sqrt(3.0)));
        gamma_C = (gamma_T/3.0);
    }

    virtual double& getDissipativeParameterT() { return gamma_T; }
    virtual void setDissipativeParameterT(double value) {
        gamma_T = value;
    }
    virtual double& getDissipativeParameterC() { return gamma_C; }
    virtual void setDissipativeParameterC(double value) {
        gamma_C = value;
    }
    /* Equilibrium parameters */
    virtual double& getEquilibriumLinkLength() { return eqLength; }
    virtual double& getEquilibriumTriangleArea() { return eqArea; }
    virtual double& getEquilibriumAngle() { return eqAngle; }
    virtual double& getEquilibriumVolume() { return eqVolume; }
    virtual double& getEquilibriumSurface() { return eqSurface; }
    virtual double& getEquilibriumTileSpan() { return eqTileSpan; }
    virtual void setEquilibriumLinkLength(double value) { pcout << "setEquilibriumLinkLength not available for ShapeMemoryModel3D."; }
    virtual void setEquilibriumTriangleArea(double value) { pcout << "setEquilibriumTriangleArea not available for ShapeMemoryModel3D."; }
    virtual void setEquilibriumAngle(double value) { pcout << "setEquilibriumAngle not available for ShapeMemoryModel3D."; }
    virtual void setEquilibriumVolume(double value) { eqVolume = value; }
    virtual void setEquilibriumSurface(double value) { eqSurface = value; }
    virtual void setEquilibriumTileSpan(double value) { eqTileSpan = value; }

    /* State parameters */
    virtual pluint& getNumberOfVertices() { return cellNumVertices; }
};

}

#include "shapeMemoryModel3D.hh"
#endif  // SHAPE_MEMORY_MODEL_3D_H
