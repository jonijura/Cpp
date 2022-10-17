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
const double EDGELENGTH=2;
const uint MESHPOINTS=EDGELENGTH*EDGELENGTH*EDGELENGTH*1000;
const double BOUNDARYLENGTH=0.1;

const uint TIMESTEPS = 100;
const double CFLCONST = 1.99;


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

void calculateOperators(Dec &dec, Derivative &d1 ,Diagonal<double> &h2, Diagonal<double> &h1i, PartMesh &pm){
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

// Calculate largest eigenvalue through power iteration
// https://en.wikipedia.org/wiki/Power_iteration
double largestEig(Sparse<double> m, uint iterc){
    Random rnd(321);
    Buffer<double> v(m.m_height);
    for(uint i=0; i<m.m_height; i++)
        v[i]=rnd.getUniform();
    Column<double> cd(v,0);
    for(uint i=0; i<iterc; i++){
        cd.scale(1.0/sqrt(cd.getLensq()));
        cd = m*cd;
    }
    cd.scale(1.0/sqrt(cd.getLensq()));
    return cd.getDot(m*cd);
}

Random rndam;
double one(const Buffer<double> &q){
    return rndam.getUint()/100.0;
}

void saveMatrix(Sparse<double> &m){
    Text txt;
    // uint lc = 0;
    // for(uint i=0; i<m.m_width; i++){
    //     for(uint j=0; j<m.m_height; j++){
    //         if(lc<m.m_col.size() && j==m.m_col[lc])
    //             txt << m.m_val[lc++] << " ";
    //         else
    //             txt << "0 ";
    //     }
    //     txt << "\n";
    // }
    for(uint i=0; i<m.m_val.size(); i++)
        txt << m.m_val[i] << " ";
    txt << endl;
    for(uint i=0; i<m.m_col.size(); i++)
        txt << m.m_col[i] << " ";;
    txt << endl;
    for(uint i=0; i<m.m_beg.size()-1; i++){
        for(int j=0; j<m.m_beg[i+1]-m.m_beg[i]; j++)
            txt << i << " ";
    }
    for(uint i=0; i<m.m_val.size()-m.m_beg[m.m_beg.size()-1]; i++)
        txt << m.m_beg.size()-1 << " ";
    txt.save("SystemMatrix.txt");
}

int main() {
    //create random delaunay triangulation of a cube
    PartMesh pm(0,1,3);
    cout << "meshing\n";
    createMesh(pm);
    //calculate opertators
    Dec dec(pm, 0,pm.getDimension());
    Derivative d1;
    Diagonal<double> h2,h1i,h2i,h1;
    dec.integrateHodge(HodgeUnit3, 0, fg_prim2, h2i).invert();
    dec.integrateHodge(HodgeUnit3, 0, fg_prim1, h1);
    cout << "making operators\n";
    calculateOperators(dec, d1,h2,h1i, pm);
    //calculate timesteps
    Sparse<double> systemMatrix;
    systemMatrix = h1i*transpose(d1)*h2*d1;
    saveMatrix(systemMatrix);
    double le = largestEig(systemMatrix, 100);
    cout << "largest eigenvalue " << le;
    double dt = CFLCONST/sqrt(le);
    cout << " timestep size " << dt << endl;
    //iterate over time and record energy norm
    Column<double> e(0.0);
	dec.integrateForm(one, 10, fg_prim1, e);
	Column<double> h(0.0);
	dec.integrateZeroForm(fg_prim2, h);
    Sparse<double> A,B;
    A.setScale(dt,h1i*transpose(d1));
    B.setScale(dt,-h2*d1);
    Text sol;
    for(uint i=0; i<TIMESTEPS; i++){
        h+=B*e;
        e+=A*h;
        double p=0.5*(e.getDot(h1*e) + h.getDot(h2i*h));
        sol << p << "\n";
    }
    sol.save("sol.txt");
}