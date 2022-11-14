#include "../../GFD/Types/Random.hpp"
#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Mesh/PartMesh.hpp"
#include "../../GFD/Discrete/Dec.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include <algorithm>
#include <iomanip>
#include <filesystem>
#include <iostream>
#include <chrono>
using namespace std;
using namespace gfd;

/**
 * Repeat experiments on asyncronous timesteps by Jukka and AVI
 */
const uint RANDOMSEED=127;

const double MINNODEDST=0.0001;
const double EDGELENGTH=0.3;
const uint MESHPOINTS=EDGELENGTH*EDGELENGTH*EDGELENGTH*1000;
const double BOUNDARYLENGTH=0.1;

const double SIMULATIONPERIOD = 100.0;
const double CFLCONST = 1.9;
const double NONUNIFORMC = 3.5;

const uint BOUNDARYFLAG = 1;

Derivative d1;
Diagonal<double> h2,h1i,h2i,h1, diriclet;
PartMesh pm(0,1,3);

void savePicture(){
    Picture pc(500,500);
    MeshDrawer md;
    md.initPicture(&pc);
    md.initPosition(Vector4(2,2,2,0), Vector4(0.5,0.5,0.5,0), Vector4(.2,-.2,0,0), Vector4(0,.2,-.2,0));
    UintSet ui(1);
    md.drawPrimalEdges(pm,Vector3(1,0,0),ui);
    pc.save("pic.bmp", true);
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
    t.save("meshstats.txt");
    savePicture();
}

// Calculate largest eigenvalue through power iteration
// https://en.wikipedia.org/wiki/Power_iteration
double largestEig(Sparse<double> m, uint iterc=100){
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

/**
 * save matrix M in sparce row format:
 * first line: nonzero values (vals)
 * second line: column indices (cols)
 * third line: row indices (rows)
 * such that M[rows[i]][cols[i]] = vals[i]
 */
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
    txt.save("SystemMatrix.txt");
}

void createOperations(Dec &dec){
    dec.integrateDerivative(fg_prim1, d1); 
    dec.integrateHodge(HodgeUnit3, 0, fg_prim1, h1);
    dec.integrateHodge(HodgeUnit3, 0, fg_prim2, h2);
    dec.integrateHodge(HodgeUnit3, 0, fg_prim1, h1i).invert();
    dec.integrateHodge(HodgeUnit3, 0, fg_prim2, h2i).invert();
    // Buffer<pair<uint,double>> edgeFlags(pm.getEdgeSize());//index out of bounds for sparse diag*sparce
    // for(uint i=0; i<pm.getEdgeSize(); i++)
    //     if(pm.getEdgeFlag(i)!=BOUNDARYFLAG)
    //         edgeFlags.push_back(make_pair(i,1.0));
    // diriclet.setSparse(pm.getEdgeSize(), edgeFlags);
    Buffer<double> edgeFlags(pm.getEdgeSize());
    for(uint i=0; i<pm.getEdgeSize(); i++)
        edgeFlags[i] = pm.getEdgeFlag(i)==BOUNDARYFLAG ? 0.0 : 1.0;
    diriclet.setFull(edgeFlags);
}

/**
 * maximal global timestep, CFL condition:
 * dt < 2/sqrt(lambda), where lambda is the largest eigenvalue
 * of the system matrix. see refs
 */
double calculateGlobalTimestep(){
    Sparse<double> systemMatrix;
    systemMatrix = h1i*transpose(d1)*h2*d1;
    // saveMatrix(systemMatrix);
    double le = largestEig(systemMatrix);
    cout << "\tlargest eigenvalue " << le;
    double dt = CFLCONST/sqrt(le);
    cout << "\n\ttimestep size " << dt << " iterations " << int(SIMULATIONPERIOD/dt);
    cout << "\n\toperations per unit time: " << (pm.getEdgeSize()+pm.getFaceSize())/dt << endl;
    return dt;
}

/**
 * Iterate with global timestep and record energynorm
 */
void iterateGlobal(Dec &dec){
    double dt = calculateGlobalTimestep();
    Sparse<double> A,B;
    A.setScale(dt,diriclet*h1i*transpose(d1));
    B.setScale(dt,-h2*d1);
    Column<double> e(0.0), h(0.0);
	dec.integrateForm(rnd, 10, fg_prim1, e);
	dec.integrateZeroForm(fg_prim2, h);
    e=diriclet*e;//initial state needs to comply with boundary conditions
    Text sol;
    uint refp = pm.findNode(Vector4(EDGELENGTH/2, EDGELENGTH/2, EDGELENGTH/2,0),0.1);
    double timesteps = SIMULATIONPERIOD/dt;
    cout << "\titerating";
    double goal = timesteps/10.0;
    for(double i=0; i<timesteps; i++){
        h+=B*e;
        double p = 0.5*((e+A*h).getDot(h1*e) +  h.getDot(h2i*h));//energy at syncronized timesteps
        e+=A*h;
        // double p=0.5*(e.getDot(h1*e) + h.getDot(h2i*h));
        sol << p << "\n";
        if(i>goal){
            cout << "*";
            goal+=timesteps/10.0;
        }
    }
    sol.save("globalTimestepEnergy.txt");
}

/**
 * Jukka chp 7.1
 * some dual areas were negative, taking absolute value was enough to fix this.
 */
void calculateLocalTimestepLimits(Buffer<double> &max_ts_e, Buffer<double> &max_ts_h){
    uint num_edges = pm.getEdgeSize();
    Column<double> e_max_val(num_edges,0.0);
    for(uint i=0; i<num_edges; i++)
        e_max_val.m_val[i] = sqrt(fabs(h1i.m_val[i]));
    Sparse<double> system_e(0.0);
    system_e = h1i*transpose(d1)*h2*d1;
    for(uint i=0; i<system_e.m_val.size(); i++)
        system_e.m_val[i] = fabs(system_e.m_val[i]);
    system_e = system_e*e_max_val;
    max_ts_e.resize(num_edges);
    for(uint i=0; i<num_edges; i++)
        max_ts_e[i] = sqrt(NONUNIFORMC*e_max_val.m_val[i]/system_e.m_val[i]);

    uint num_faces = pm.getFaceSize();
    Column<double> h_max_val(num_faces,0.0);
    for(uint i=0; i<num_faces; i++)
        h_max_val.m_val[i] = sqrt(fabs(h2.m_val[i]));
    Sparse<double> system_h(0.0);
    system_h = h2*d1*h1i*transpose(d1);
    for(uint i=0; i<system_h.m_val.size(); i++)
        system_h.m_val[i] = fabs(system_h.m_val[i]);
    system_h = system_h*h_max_val;
    max_ts_h.resize(num_faces);
    for(uint i=0; i<num_faces; i++)
        max_ts_h[i] = sqrt(NONUNIFORMC*h_max_val.m_val[i]/system_h.m_val[i]);
    
    cout << "\te max timestep: " << max_ts_e.max() << " min: " << max_ts_e.min();
    cout << "\n\th max timestep: " << max_ts_h.max() << " min: " << max_ts_h.min();
}

/**
 * jukka chp 7.2 find minimal timestep divisors in {1,3,9,27,...}
 */ 
void calculateDivisionCoeficcients(Buffer<double> &et, Buffer<double> &ht, double dt, Buffer<uint> &ecoeffs, Buffer<uint> &hcoeffs){
    int count[10] = {};

    uint num_edges = et.size();
    ecoeffs.resize(num_edges);
    for(uint i=0; i<num_edges; i++){
        Buffer<uint> edge_faces = pm.getEdgeFaces(i);
        double mindt = et[i];
        for(uint j=0; j<edge_faces.size(); j++)
            mindt = min(mindt, ht[edge_faces[j]]);
        int ts_divisor = 1;
        int pow = 0;
        while(mindt*ts_divisor<dt){
            ts_divisor*=3;
            pow++;
        }
        count[pow]++;
        ecoeffs[i]=ts_divisor;
    }

    uint num_faces = ht.size();
    hcoeffs.resize(num_faces);
    for(uint i=0; i<num_faces; i++){
        Buffer<uint> face_edges = pm.getFaceEdges(i);
        double mindt = ht[i];
        for(uint j=0; j<face_edges.size(); j++)
            mindt = min(mindt, et[face_edges[j]]);
        int ts_divisor = 1;
        int pow = 0;
        while(mindt*ts_divisor<dt){
            ts_divisor*=3;
            pow++;
        }
        count[pow]++;
        hcoeffs[i]=ts_divisor;
    }
    //statistics
    int factor = 1;
    double avg_ts;
    int operations_per_dt = 0;
    cout<< "\tfactor: ";
    for(int i=0; i<10; i++){
        if(count[i]!=0)
            cout << factor << ": " << count[i] << " ";
        avg_ts += count[i]*dt/factor;
        operations_per_dt += count[i]*factor;
        factor*=3;
    }
    cout << "\n\taverage timestep size: " << avg_ts/(num_edges+num_faces) << "\n\toperations per unit time: " << operations_per_dt/dt << endl;
}

/**
 * sort e and h update times, Jukka chp 7.2
 * order<tuple< is e field?, indice, update order>>
 */
void arrangeUpdates(vector<tuple<bool,uint,double>> &order, Buffer<uint> &e_ts_div_coeffs, Buffer<uint> &h_ts_div_coeffs){
    int size = e_ts_div_coeffs.sum()+h_ts_div_coeffs.sum();
    order.resize(size);
    int loc = 0;
    for(uint i=0; i<e_ts_div_coeffs.size(); i++)
        for(uint j=0; j<e_ts_div_coeffs[i]; j++)
            order[loc++]=make_tuple(true,i,(j+0.5)/e_ts_div_coeffs[i]);
    for(uint i=0; i<h_ts_div_coeffs.size(); i++)
        for(uint j=0; j<h_ts_div_coeffs[i]; j++)
            order[loc++]=make_tuple(false,i,(j+0.0)/h_ts_div_coeffs[i]);
    sort(order.begin(), order.end(), [](const tuple<bool,uint,double> & a, const tuple<bool,uint,double> & b) -> bool
    { 
        return get<2>(a) < get<2>(b); 
    });
}

/**
 * debuggin, print sparce matrix and a vector of its columns
 * code that I might need again
 */
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

/**
 * Split the matrices into vector of rows
 */
void makeLocalOperators(vector<Column<double>> &oA, vector<Column<double>> &oB){
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

/**
 * Iterate using local timestepping and record energy norm, Jukka chp7
 */
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
    Column<double> e(0.0), h(0.0), hsync(0.0);
	dec.integrateForm(rnd, 10, fg_prim1, e);
	dec.integrateZeroForm(fg_prim2, h);
    e=e*diriclet;
    Text sol;
    cout << "\titerating";
    double goal = steps/10.0;
    for(double i=0; i<steps; i++){
        for(auto a : updateOrder){
            uint n = get<1>(a);
            if(get<0>(a))
                e.m_val[n]+=dt/dce[n]*A[n].getDot(h);
            else
                h.m_val[n]+=dt/dch[n]*B[n].getDot(e);
        }
        hsync = h;
        for(uint i=0; i<h.m_val.size(); i++){
            hsync.m_val[i] += 0.5*dt/dch[i]*B[i].getDot(e);
        }
        double p=0.5*(e.getDot(h1*e) + hsync.getDot(h2i*hsync));
        sol << p << "\n";
        if(i>goal){
            cout << "*";
            goal+=steps/10.0; 
        }
    }
    sol.save("localTimestepEnergy.txt");
}

std::chrono::_V2::system_clock::time_point starttime;
void toc(){
    cout << "\tcomplete in: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << endl;
    starttime = chrono::system_clock::now();
}

int main() {
    starttime = chrono::system_clock::now();
    //create random delaunay triangulation of a cube
    cout << "meshing ... \n";
    createMesh(); 
    toc(); 
    //calculate opertators
    cout << "making operators ... \n";
    Dec dec(pm, 0, pm.getDimension());
    createOperations(dec);
    toc();
    //iterate over time and record energy norm
    cout << "timestepping with global timestep ... \n";
    iterateGlobal(dec); 
    toc();
    cout << "timestepping with local timestep ... \n";
    iterateLocal1(dec);
    toc();
}