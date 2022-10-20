#include "../../GFD/Types/Random.hpp"
#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Mesh/PartMesh.hpp"
#include "../../GFD/Discrete/Dec.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include <algorithm>
#include <iomanip>
#include <filesystem>
#include <iostream>
using namespace std;
using namespace gfd;

/**
 * Repeat experiments on asyncronous timesteps by Jukka and AVI
 * alustaminen satunnaisfunktiolla reunaehdon mukaisesti, reunoilla nollaa?
 */
const uint RANDOMSEED=27;

const double MINNODEDST=0.0001;
const double EDGELENGTH=1.0;
const uint MESHPOINTS=EDGELENGTH*EDGELENGTH*EDGELENGTH*1000;
const double BOUNDARYLENGTH=0.1;

const double SIMULATIONPERIOD = 30.0;
const double CFLCONST = 1.9;
const double NONUNIFORMC = 1.9;

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
    cout << "\tlargest eigenvalue " << le;
    double dt = CFLCONST/sqrt(le);
    cout << "\n\ttimestep size " << dt << " iterations " << int(SIMULATIONPERIOD/dt);
    cout << "\n\toperations per unit time: " << (pm.getEdgeSize()+pm.getFaceSize())/dt << endl;
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
    e=diriclet*e;
	dec.integrateZeroForm(fg_prim2, h);
    uint source = pm.findNode(Vector4(EDGELENGTH/2.0,EDGELENGTH/2.0,EDGELENGTH/2.0,0),0.1);
    cout << "\tsource position: " << source <<  " " << pm.getNodePosition3(source).x << endl;
    const double fm = 5.0;
    const double lev = 20/(PI*fm);
    Text sol;
    Text pulse;
    double timesteps = SIMULATIONPERIOD/dt;
    for(double i=0; i<timesteps; i++){
        h+=B*e;
        double p = 0.5*((e+A*h).getDot(h1*e) +  h.getDot(h2i*h));
        e+=A*h;
        double j = exp(-(i*dt-3*lev)*(i*dt-3*lev)/(lev*lev))*sin(PIx2*fm*(i*dt-3*lev));
        // e.m_val[source]+=j;
        // double p=0.5*(e.getDot(h1*e) + h.getDot(h2i*h));
        sol << p << "\n";
        pulse << j << "\n";
    }
    sol.save("build\\globalTimestepEnergy.txt");
    pulse.save("build\\sourcePulse.txt");
}

void calculateLocalTimestepLimits(Buffer<double> &et, Buffer<double> &ht){
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
    
    cout << "\te max timestep: " << et.max() << " min: " << et.min();
    cout << "\n\th max timestep: " << ht.max() << " min: " << ht.min();
}

void calculateDivisionCoeficcients(Buffer<double> &et, Buffer<double> &ht, double dt, Buffer<uint> &ecoeffs, Buffer<uint> &hcoeffs){
    int count[10] = {};

    uint ne = et.size();
    ecoeffs.resize(ne);
    for(uint i=0; i<ne; i++){
        Buffer<uint> nbf = pm.getEdgeFaces(i);
        double mindt = et[i];
        for(uint j=0; j<nbf.size(); j++)
            mindt = min(mindt, ht[nbf[j]]);
        int dc = 1;
        int pow = 0;
        while(mindt*dc<dt){
            dc*=3;
            pow++;
        }
        count[pow]++;
        ecoeffs[i]=dc;
    }

    uint nf = ht.size();
    hcoeffs.resize(nf);
    for(uint i=0; i<nf; i++){
        Buffer<uint> nbe = pm.getFaceEdges(i);
        double mindt = ht[i];
        for(uint j=0; j<nbe.size(); j++)
            mindt = min(mindt, et[nbe[j]]);
        int dc = 1;
        int pow = 0;
        while(mindt*dc<dt){
            dc*=3;
            pow++;
        }
        count[pow]++;
        hcoeffs[i]=dc;
    }

    int factor = 1;
    double avgTs;
    int ops = 0;
    cout<< "\tfactor: ";
    for(int i=0; i<10; i++){
        if(count[i]!=0)
            cout << factor << ": " << count[i] << " ";
        avgTs += count[i]*dt/factor;
        ops += count[i]*factor;
        factor*=3;
    }
    cout << "\n\taverage timestep size: " << avgTs/(ne+nf) << "\n\toperations per unit time: " << ops/dt << endl;
}

void arrangeUpdates(vector<tuple<bool,uint,double>> &order, Buffer<uint> &dce, Buffer<uint> &dch){
    int size = dce.sum()+dch.sum();
    order.resize(size);
    int loc = 0;
    for(uint i=0; i<dce.size(); i++)
        for(uint j=0; j<dce[i]; j++)
            order[loc++]=make_tuple(true,i,(j+0.5)/dce[i]);
    for(uint i=0; i<dch.size(); i++)
        for(uint j=0; j<dch[i]; j++)
            order[loc++]=make_tuple(false,i,(j+0.0)/dch[i]);
    sort(order.begin(), order.end(), [](const tuple<bool,uint,double> & a, const tuple<bool,uint,double> & b) -> bool
    { 
        return get<2>(a) < get<2>(b); 
    });
}

void printLocalOps(Sparse<double> &A, vector<Column<double>> &oA){
    if(A.m_height>=30){
        cout << "refusing to print something this big!\n";
        return; 
    }
    uint lc = 0;
    for(uint i=0; i<A.m_height; i++){
        for(uint j=0; j<A.m_width; j++){
            if(lc<A.m_col.size() && j==A.m_col[lc])
                cout << " " << setprecision(3) << setfill(' ') << setw(5) << A.m_val[lc++];
            else
                cout << "     0";
        }
        cout << "\n";
    }
    cout << "\n\n";
    for(auto a : oA){
        lc=0;
        for(uint j=0; j<a.m_height; j++){
            if(lc<a.m_row.size() && j==a.m_row[lc])
                cout << " " << setprecision(3) << setfill(' ') << setw(5) << a.m_val[lc++];
            else
                cout << "     0";
        }
        cout << "\n";
    }
}

void makeLocalOperators(vector<Column<double>> &oA, vector<Column<double>> &oB){
    Buffer<double> edgeFlags(pm.getEdgeSize());
    for(uint i=0; i<pm.getEdgeSize(); i++)
        edgeFlags[i]=pm.getEdgeFlag(i)==BOUNDARYFLAG ? 0 : 1;
    Diagonal<double> diriclet(edgeFlags,0.0);//TODO: make sparse to save time
    Sparse<double> A,B;
    A=diriclet*h1i*transpose(d1);
    B=-h2*d1;
    for(uint i=0; i<A.m_height; i++){
        Buffer<pair<uint, double>> vals;
        uint end = i+1==A.m_height ? A.m_val.size() : A.m_beg[i+1];
        for(uint j=A.m_beg[i]; j<end; j++)
            vals.push_back(make_pair(A.m_col[j],A.m_val[j]));
        Column<double> row(A.m_width, vals, 0.0);
        oA.push_back(row);
    }
    for(uint i=0; i<B.m_height; i++){ 
        Buffer<pair<uint, double>> vals;
        uint end = i+1==B.m_height ? B.m_val.size() : B.m_beg[i+1];
        for(uint j=B.m_beg[i]; j<end; j++)
            vals.push_back(make_pair(B.m_col[j],B.m_val[j]));
        Column<double> row(B.m_width, vals, 0.0);
        oB.push_back(row);
    }
}

void iterateLocal1(Dec &dec){
    Buffer<double> et,ht;
    calculateLocalTimestepLimits(et,ht);
    double dt = 0.15*(et.max()+ht.max());//81.0*min(et.min(),ht.min());//
    double steps = SIMULATIONPERIOD/dt; 
    cout << "\n\tsyncronous timestep size: " << dt  << " iteration count " << int(steps) << endl;
    Buffer<uint> dce, dch;
    calculateDivisionCoeficcients(et,ht,dt, dce, dch);
    vector<tuple<bool,uint,double>> updateOrder;
    arrangeUpdates(updateOrder, dce, dch);
    vector<Column<double>> A,B;
    makeLocalOperators(A,B); 
    Column<double> e(0.0), h(0.0);
	dec.integrateForm(rnd, 10, fg_prim1, e);
	dec.integrateZeroForm(fg_prim2, h);
    Buffer<double> edgeFlags(pm.getEdgeSize());
    for(uint i=0; i<pm.getEdgeSize(); i++)
        edgeFlags[i]=pm.getEdgeFlag(i)==BOUNDARYFLAG ? 0 : 1;
    Diagonal<double> diriclet(edgeFlags,0.0);
    e=e*diriclet;
    uint source = pm.findNode(Vector4(EDGELENGTH/2.0,EDGELENGTH/2.0,EDGELENGTH/2.0,0),0.1);
    cout << "\tsource position: " << source <<  " " << pm.getNodePosition3(source).x << endl;
    uint n;
    Text sol;
    for(double i=0; i<steps; i++){
        for(auto a : updateOrder){
            n = get<1>(a);
            if(get<0>(a))
                e.m_val[n]+=dt/dce[n]*A[n].getDot(h);
            else
                h.m_val[n]+=dt/dch[n]*B[n].getDot(e);
        }
        double p=0.5*(e.getDot(h1*e) + h.getDot(h2i*h));
        sol << p << "\n";
    }
    sol.save("build\\localTimestepEnergy.txt");
}

int main() {
    //create random delaunay triangulation of a cube
    cout << "meshing ... \n";
    createMesh();
    //calculate opertators
    cout << "making operators ... \n";
    Dec dec(pm, 0, pm.getDimension());
    createOperations(dec);
    //iterate over time and record energy norm
    cout << "timestepping with global timestep ... \n";
    iterateGlobal(dec);
    cout << "timestepping with local timestep ... \n";
    iterateLocal1(dec);
}