#include "../../GFD/Types/Random.hpp"
#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Mesh/PartMesh.hpp"
#include "../../GFD/Discrete/Dec.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include <filesystem>
#include <iostream>
using namespace std;
using namespace gfd;

const uint RANDOMSEED=222;

const double MINNODEDST=0.0001;
const double EDGELENGTH=1;
const uint MESHPOINTS=EDGELENGTH*EDGELENGTH*EDGELENGTH*1000;
const double BOUNDARYLENGTH=0.1;

const uint TIMESTEPS = 20000;
const double CFLCONST = 1.9;

void savePicture(Mesh &mesh){
    Picture pc(500,500);
    MeshDrawer md;
    md.initPicture(&pc);
    md.initPosition(Vector4(2,2,2,0), Vector4(0.5,0.5,0.5,0), Vector4(.2,-.2,0,0), Vector4(0,.2,-.2,0));
    UintSet ui(1);
    md.drawPrimalEdges(mesh,Vector3(1,0,0),ui);
    pc.save("build\\pic.bmp", true);
}

void createMesh(PartMesh &pm){
    BuilderMesh msh(3);
    double n_b = round(EDGELENGTH/BOUNDARYLENGTH);
    double dxyz = EDGELENGTH/n_b;
    double bl = EDGELENGTH;
    Random rnd(RANDOMSEED);
    for(double i=0; i<=n_b; i++){
        for(double j=0; j<=n_b; j++){
            msh.insertNode(Vector4(bl,i*dxyz, j*dxyz,0),0,0,false);
            msh.insertNode(Vector4(i*dxyz,0, j*dxyz,0),0,0,false);
            msh.insertNode(Vector4(i*dxyz,bl, j*dxyz,0),0,0,false);
            msh.insertNode(Vector4(i*dxyz, j*dxyz,0,0),0,0,false);
            msh.insertNode(Vector4(i*dxyz, j*dxyz,bl,0),0,0,false);
            msh.insertNode(Vector4(0,i*dxyz, j*dxyz,0),0,0,false);
        }
    }
    for(uint i=0; i<MESHPOINTS; i++){
        Vector4 pt(rnd.getUniform(),rnd.getUniform(),rnd.getUniform(),0);
        pt*=EDGELENGTH;
        if((msh.getNodePosition(msh.searchNode(pt))-pt).lensq()>MINNODEDST)
            msh.insertNode(pt,0,0,false);
    }
    msh.fillBoundaryFlags(1);
    pm.swap(msh);
    Text t;
    pm.writeStatistics(t);
    t.save("build\\meshstats.txt");
    savePicture(pm);
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
double rnd(const Buffer<double> &q){
    return 0;
    return rndam.getUniform();
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
        for(uint j=0; j<m.m_beg[i+1]-m.m_beg[i]; j++)
            txt << i << " ";
    }
    for(uint i=0; i<m.m_val.size()-m.m_beg[m.m_beg.size()-1]; i++)
        txt << m.m_beg.size()-1 << " ";
    txt.save("build\\SystemMatrix.txt");
}

int main() {
    //create random delaunay triangulation of a cube
    PartMesh pm(0,1,3);
    cout << "meshing\n";
    createMesh(pm);
    //calculate opertators
    Dec dec(pm, 0, pm.getDimension());
    Derivative d1;
    Diagonal<double> h2,h1i,h2i,h1;
    cout << "making operators\n";
    dec.integrateDerivative(fg_prim1, d1);
    dec.integrateHodge(HodgeUnit3, 0, fg_prim1, h1);
    dec.integrateHodge(HodgeUnit3, 0, fg_prim2, h2);
    dec.integrateHodge(HodgeUnit3, 0, fg_prim1, h1i).invert();
    dec.integrateHodge(HodgeUnit3, 0, fg_prim2, h2i).invert();
    for(uint i=0; i<pm.getEdgeSize(); i++)
        h1i.m_val[i]=pm.getEdgeFlag(i)>0 ? 0 : h1i.m_val[i];   
    // calculateOperators(dec, d1,h2,h1i, pm);
    //calculate timesteps
    Sparse<double> systemMatrix;
    systemMatrix = h1i*transpose(d1)*h2*d1;
    saveMatrix(systemMatrix);
    double le = largestEig(systemMatrix, 100);
    cout << "largest eigenvalue " << le;
    double dt = CFLCONST/sqrt(le);
    cout << " timestep size " << dt << endl;
    //iterate over time and record energy norm
    Column<double> e(0.0), h(0.0), hsync(0.0);
	dec.integrateForm(rnd, 10, fg_prim1, e);
	dec.integrateZeroForm(fg_prim2, h);
    Sparse<double> A,B,C;
    A.setScale(dt,h1i*transpose(d1));
    B.setScale(dt,-h2*d1);
    C.setScale(0.5*dt,-h2*d1);
    uint source = pm.findNode(Vector4(EDGELENGTH/2.0,EDGELENGTH/2.0,EDGELENGTH/2.0,0),0.1);
    cout << source <<  " " << pm.getNodePosition3(source).x;
    const double fm = 5.0;
    const double lev = 20/(PI*fm);  
    Text sol;
    for(uint i=0; i<TIMESTEPS; i++){
        h+=B*e;
        // double p = 0.5*((e+A*h).getDot(h1*e) +  h.getDot(h2i*h));
        e+=A*h;
        e.m_val[source]+=exp(-(i*dt-3*lev)*(i*dt-3*lev)/(lev*lev))*sin(PIx2*fm*(i*dt-3*lev));
        double p=0.5*(e.getDot(h1*e) + h.getDot(h2i*h));
        // hsync = h + C*e; //has very little effect? not even the right way to calculate energy, but gives the same result
        // double p=0.5*(e.getDot(h1*e) + hsync.getDot(h2i*hsync));
        sol << p << "\n";
    }
    sol.save("build\\sol.txt");
}