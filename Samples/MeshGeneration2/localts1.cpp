#include "../../GFD/Types/Random.hpp"
#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Mesh/PartMesh.hpp"
#include "../../GFD/Discrete/Dec.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include <filesystem>
#include <iostream>
using namespace std;
using namespace gfd;

const uint RANDOMSEED=123;
const double MINNODEDST=0.01;
const double EDGELENGTH=0.1;
const uint MESHPOINTS=EDGELENGTH*EDGELENGTH*EDGELENGTH*1000;
const double BOUNDARYLENGTH=0.1;


void createMesh(PartMesh &pm){
    BuilderMesh msh(3);
    msh.createGrid(Vector4(0,0,0,0), EDGELENGTH*Vector4(1,1,1,0), BOUNDARYLENGTH*Vector4(1,1,1,1));
    Random rnd(RANDOMSEED);
    for(uint i=0; i<MESHPOINTS; i++){
        Vector4 pt(rnd.getUniform(),rnd.getUniform(),rnd.getUniform(),0);
        if((msh.getNodePosition(msh.searchNode(pt))-pt).lensq()>MINNODEDST)
            msh.insertNode(pt,0,0,false);
    }
    msh.fillBoundaryFlags(1);
    pm.swap(msh);
    Text t;
    pm.writeStatistics(t);
    t.save("build\\meshstats.txt");
}

void calculateOperators(Dec &dec, Derivative &d1 ,Diagonal<double> &h2,Diagonal<double> &h1i, PartMesh &pm){
    dec.integrateDerivative(fg_prim1, d1);
    //just guessing HodgeUnit3 -> calculations in 3 dimensions
    //lover values leave some points in the diagonal unaffected
    //higher values seem to access unmapped memory registers though buffer -.-
    dec.integrateHodge(HodgeUnit3, 0, fg_prim2, h2);
    dec.integrateHodge(HodgeUnit3, 0, fg_prim1, h1i).invert();//integrating over fg_dual2 doesnt work? matrix is wrong size
    // sanity checks
    // h2.printShape();
    // h1i.printShape();
    Buffer<double> h2b = h2.getBuffer();
    cout << "h2 size: " << h2b.size() << " expected: " << pm.getFaceSize() << endl;
    for(uint i=0; i<h2b.size(); i++)
        if(abs(h2b[i]-pm.getFaceHodge(i))>1e-6)
            cout << "unexpected value: " << h2b[i] << " "<< pm.getFaceHodge(i) << "\n";
    Buffer<double> h1ib = h1i.getBuffer();
    cout << "h1i size: "<< h1ib.size() << " expected: "<<pm.getEdgeSize() << endl;
    for(uint i=0; i<h1ib.size(); i++)
        if(abs(h1ib[i]- 1.0/pm.getEdgeHodge(i))>1e-6)
            cout << "unexpected value: " << h1ib[i] << " "<< 1.0/pm.getEdgeHodge(i) << "\n";
}

double largestEig(Sparse<double> m){
    m.printData();
    cout << m.m_val.size();
    return 0.0;
}

int main() {
    //create random delaunay triangulation of a cube
    PartMesh pm(0,1,3);
    cout << "meshing\n";
    createMesh(pm);
    //calculate opertators
    Dec dec(pm, 0,pm.getDimension());
    Derivative d1;
    Diagonal<double> h2,h1i;
    cout << "making operators\n";
    calculateOperators(dec, d1,h2,h1i, pm);
    //calculate timesteps
    Sparse<double> systemMatrix;
    systemMatrix = h1i*transpose(d1)*h2*d1;
    double dt = 1.99/sqrt(largestEig(systemMatrix));
    //iterate over time and record energy norm
    
}