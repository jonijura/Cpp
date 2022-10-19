#include "../../GFD/Types/Random.hpp"
#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Mesh/PartMesh.hpp"
#include "../../GFD/Discrete/Dec.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include <filesystem>
#include <iostream>
using namespace std;
using namespace gfd;

const uint RANDOMSEED=27;

const double MINNODEDST=0.0001;
const double EDGELENGTH=0.3;
const uint MESHPOINTS=EDGELENGTH*EDGELENGTH*EDGELENGTH*1000;
const double BOUNDARYLENGTH=0.1;

const uint TIMESTEPS = 50;
const double CFLCONST = 1.9;
const double NONUNIFORMC = 3.9;

const uint BOUNDARYFLAG = 1;

Derivative d1;
Diagonal<double> h2,h1i,h2i,h1;
PartMesh pm(0,1,3);

void savePicture(){
    Picture pc(500,500);
    MeshDrawer md;
    md.initPicture(&pc);
    md.initPosition(Vector4(2,2,2,0), Vector4(0.5,0.5,0.5,0), Vector4(.2,-.2,0,0), Vector4(0,.2,-.2,0));
    UintSet ui(1);
    md.drawPrimalEdges(pm,Vector3(1,0,0),ui);
    pc.save("build\\pic.bmp", true);
}

void createMesh(){
    BuilderMesh msh(3);
    double n_b = round(EDGELENGTH/BOUNDARYLENGTH);
    double dxyz = EDGELENGTH/n_b;
    double bl = EDGELENGTH;
    Random rnd(RANDOMSEED);
    for(double i=0; i<=n_b; i++){//boundary nodes
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
    msh.fillBoundaryFlags(BOUNDARYFLAG);

    // cout << "USING CUBE MESH" <<endl;
    // msh.clear();
    // msh.createGrid(Vector4(0,0,0,0), EDGELENGTH*Vector4(1,1,1,0), BOUNDARYLENGTH*Vector4(1,1,1,1));

    pm.swap(msh);
    Text t;
    pm.writeStatistics(t);
    t.save("build\\meshstats.txt");
    savePicture();
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

void createOperations(Dec &dec){
    dec.integrateDerivative(fg_prim1, d1);
    dec.integrateHodge(HodgeUnit3, 0, fg_prim1, h1);
    dec.integrateHodge(HodgeUnit3, 0, fg_prim2, h2);
    dec.integrateHodge(HodgeUnit3, 0, fg_prim1, h1i).invert();
    dec.integrateHodge(HodgeUnit3, 0, fg_prim2, h2i).invert();
}

double calculateGlobalTimestep(){
    Sparse<double> systemMatrix;
    systemMatrix = h1i*transpose(d1)*h2*d1;
    saveMatrix(systemMatrix);
    double le = largestEig(systemMatrix, 100);
    cout << "largest eigenvalue " << le;
    double dt = CFLCONST/sqrt(le);
    cout << " timestep size " << dt << endl;
    return dt;
}

void iterateGlobal(Dec &dec){
    double dt = calculateGlobalTimestep();
    Buffer<double> edgeFlags(pm.getEdgeSize());
    for(uint i=0; i<pm.getEdgeSize(); i++)
        edgeFlags[i]=pm.getEdgeFlag(i)==BOUNDARYFLAG ? 0 : 1;
    Diagonal<double> diriclet(edgeFlags,0.0);
    Sparse<double> A,B;
    A.setScale(dt,diriclet*h1i*transpose(d1));
    B.setScale(dt,-h2*d1);
    Column<double> e(0.0), h(0.0);
	dec.integrateForm(rnd, 10, fg_prim1, e);
	dec.integrateZeroForm(fg_prim2, h);
    uint source = pm.findNode(Vector4(EDGELENGTH/2.0,EDGELENGTH/2.0,EDGELENGTH/2.0,0),0.1);
    cout << "source position: " << source <<  " " << pm.getNodePosition3(source).x << endl;
    const double fm = 5.0;
    const double lev = 20/(PI*fm);
    Text sol;
    Text pulse;
    for(double i=0; i<TIMESTEPS; i++){
        h+=B*e;
        // double p = 0.5*((e+A*h).getDot(h1*e) +  h.getDot(h2i*h));
        e+=A*h;
        double j = exp(-(i*dt-3*lev)*(i*dt-3*lev)/(lev*lev))*sin(PIx2*fm*(i*dt-3*lev));
        e.m_val[source]+=j;
        double p=0.5*(e.getDot(h1*e) + h.getDot(h2i*h));
        sol << p << "\n";
        pulse << j << "\n";
    }
    pulse.save("build\\globalTimestepEnergy.txt");
    sol.save("build\\sourcePulse.txt");
}

void calculateLocalTimesteps(Buffer<double> &et, Buffer<double> &ht){
    uint en = pm.getEdgeSize();
    Column<double> maxet(en,0.0);
    for(uint i=0; i<en; i++)
        maxet.m_val[i] = sqrt(fabs(h1i.m_val[i]));
    Sparse<double> systeme(0.0);
    systeme = h1i*transpose(d1)*h2*d1;
    for(uint i=0; i<systeme.m_val.size(); i++)
        systeme.m_val[i] = fabs(systeme.m_val[i]);
    systeme = systeme*maxet;
    et.resize(en);
    for(uint i=0; i<en; i++)
        et[i] = sqrt(NONUNIFORMC*maxet.m_val[i]/systeme.m_val[i]);

    // h1.m_val.print();
    // maxet.m_val.print();
    // systeme.m_val.print();
    // et.print();

    uint ef = pm.getFaceSize();
    Column<double> maxht(ef,0.0);
    for(uint i=0; i<ef; i++)
        maxht.m_val[i] = sqrt(fabs(h2.m_val[i]));
    Sparse<double> systemh(0.0);
    systemh = h2*d1*h1i*transpose(d1);
    for(uint i=0; i<systemh.m_val.size(); i++)
        systemh.m_val[i] = fabs(systemh.m_val[i]);
    systemh = systemh*maxht;
    ht.resize(ef);
    for(uint i=0; i<ef; i++)
        ht[i] = sqrt(NONUNIFORMC*maxht.m_val[i]/systemh.m_val[i]);
}

void iterateLocal1(Dec &dec){
    Buffer<double> et,ht;
    calculateLocalTimesteps(et,ht);
    cout << "e max timestep: " << et.max() << " min: " << et.min();
    cout << "\nh max timestep: " << ht.max() << " min: " << ht.min() << endl;
}

int main() {
    //create random delaunay triangulation of a cube
    cout << "meshing\n\n";
    createMesh();
    //calculate opertators
    cout << "making operators\n\n";
    Dec dec(pm, 0, pm.getDimension());
    createOperations(dec);
    //iterate over time and record energy norm
    cout << "timestepping\n\n";
    iterateGlobal(dec);
    iterateLocal1(dec);
}