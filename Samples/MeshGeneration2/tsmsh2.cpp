#include "../../GFD/Types/Random.hpp"
#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Discrete/Dec.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include <filesystem>
#include <iostream>
#include <queue>
using namespace std;
using namespace gfd;

PartMesh mesh(0,1,3);
const uint BOUNDARYFLAG = 1;

/**
 * interpolate discrete 1-forms on mesh vertices
*/
Vector3 getField(const Buffer<double> &val, const uint node) {
	const Buffer<uint> par = mesh.getNodeEdges(node);
	Matrix3 A(0,0,0,0,0,0,0,0,0);
	Vector3 b(0,0,0);
	for(uint i=0; i<par.size(); i++) {
		const Vector3 v = mesh.getEdgeVector3(par[i]);
		A += v.outerProduct();
		b += v * val[par[i]];
	}
	return A.inverse() * b;
}

void savePicture(){
    Picture pc(500,500);
    MeshDrawer md;
    md.initPicture(&pc);
    Vector4 campos = Vector4(PI/2.0, -3, PI/2,0); 
    md.initPosition(campos, Vector4(PI/2.0, PI/2.0, PI/2,0), 
                    Vector4(.1,0,0,0), Vector4(0,.01,.1,0));
    md.drawPrimalEdges(mesh,Vector3(1,0,0), UintSet(BOUNDARYFLAG));
    md.drawPrimalEdges(mesh,Vector3(0,1,0), UintSet(3));
    md.drawPrimalEdges(mesh,Vector3(0,0,1), UintSet(4));
    md.drawPrimalEdges(mesh,Vector3(0,1,1), UintSet(5));
    md.drawPrimalEdges(mesh,Vector3(1,0,1), UintSet(6));
    pc.save("ts2.bmp", true);
}

/**
 * update the list and detail of edges that can be solved via dF=0 or d*F=0
*/
void updateComplete(uint edg, Buffer<pair<uint, uint>> &marks, queue<uint> &updates){
    Buffer<uint> faces = mesh.getEdgeFaces(edg);
    for(auto a : faces){
        Buffer<uint> faceEdges = mesh.getFaceEdges(a);
        uint missing = 0;
        for(auto b : faceEdges)
            if(marks[b].first!=3)missing++;
        if(missing==1){
            for(auto b : faceEdges){
                if(marks[b].first==0){
                    updates.push(b);
                    marks[b]={2,a};
                    // mesh.setEdgeFlag(b,4);
                }
            }
        }
    }

    Buffer<uint> nodes = mesh.getEdgeNodes(edg);
    for(auto a : nodes){
        if(mesh.getNodeFlag(a)==1){//using d*F=0 on the boundary assumes homogeneous von neumann boundary conditions.
            continue;
        }
        Buffer<uint> nodeEdges = mesh.getNodeEdges(a);
        uint missing = 0;
        for(auto b : nodeEdges)
            if(marks[b].first!=3)missing++;
        if(missing==1){
            for(auto b : nodeEdges){
                if(marks[b].first==0){
                    updates.push(b);
                    marks[b]={1,a};
                    // mesh.setEdgeFlag(b,4);
                }
            }
        }
    }
}

int main() {
    BuilderMesh bm;
    double tmax = PI+2;
    const double dxy = 0.3;
    const double dt = 0.3;
    // bm.createBccGrid(Vector3(0,0,0), Vector3(PI,PI,PI), 0.4);
    bm.createGrid(Vector4(0,0,0,0), Vector4(PI,PI,tmax,0), Vector4(dxy,dxy,dt,1));
    // bm.createGrid(Vector4(0,0,0,0), Vector4(1,1,1,0), Vector4(0.5,0.5,0.5,1));
    bm.setMetric(SymMatrix4(1,0,1,0,0,-1,0,0,0,0));
    bm.fillBoundaryFlags(BOUNDARYFLAG);
    mesh.swap(bm);
    //status.first: 1,2,3 = node, face, solved
    //second: facenum/nodenum
    Buffer<pair<uint, uint>> marks(mesh.getEdgeSize(), {0,0});
    queue<uint> updates;
    Column<double> sol(mesh.getEdgeSize(), 0.0);
    //edges at the bottom row
    Buffer<uint> initialValues;
    for(uint i = 0; i<mesh.getEdgeSize(); i++){
        auto a = mesh.getEdgePosition3(i);
        if(a.z == 0 || a.x==0 || a.y==0 || a.x==PI || a.y==PI){
            initialValues.push_back(i);
            marks[i].first = 3;
            // mesh.setEdgeFlag(i,3);
        }
        else if(abs(a.z-dt/2)<dt/4){//first timelike edges/ initial dual values
            marks[i].first = 3;
            // mesh.setEdgeFlag(i,4);
            auto nds = mesh.getEdgeNodes(i);
            auto end = mesh.getNodePosition3(nds[0]);
            auto st = mesh.getNodePosition3(nds[1]);
            sol.m_val[i] = (int)mesh.getEdgeIncidence(i,nds[0])*sin(end.x)*sin(end.y)*sin(sqrt(2.0)*end.z) / sqrt(2.0)
                        +(int)mesh.getEdgeIncidence(i,nds[1])*sin(st.x)*sin(st.y)*sin(sqrt(2.0)*st.z) / sqrt(2.0);
        }
    }
    savePicture();
    for(auto a : initialValues)
        updateComplete(a, marks, updates);
    cout << updates.size();
    Dec dec(mesh, 0, mesh.getDimension());
    Derivative d1, d0T;
    d1 = dec.integrateDerivative(fg_prim1,d1);
    d0T = dec.integrateDerivative(fg_prim0,d0T);
    d0T = d0T.setTranspose(d0T);
    //integrate hodge is too difficult to use, this seems to account for negative hodge-elements
    Buffer<double> h1val(mesh.getEdgeSize());
    for(uint i = 0 ; i<mesh.getEdgeSize(); i++){
        h1val[i] = mesh.getEdgeHodge(i);
        // if(h1val[i]==0.125)mesh.setEdgeFlag(i, 3);
        // else if(h1val[i]==0.25)mesh.setEdgeFlag(i, 4);
        // else if(h1val[i]==-0.125)mesh.setEdgeFlag(i, 5);
        // else if(h1val[i]==-0.25)mesh.setEdgeFlag(i, 6);
    }
    Diagonal<double> h1(h1val, 0.0);
    //timestepping
    while(updates.size()!=0){
        uint ind = updates.front();
        uint j = marks[ind].second;
        updates.pop();
        if(marks[ind].first==2){
            marks[ind].first=3;
            sol.m_val[ind] = -(int)mesh.getFaceIncidence(j, ind)*(d1*sol).m_val[j];
            // mesh.setEdgeFlag(ind,5);
        }
        else if(marks[ind].first==1){
            marks[ind].first=3;
            sol.m_val[ind] = -(int)mesh.getEdgeIncidence(ind, j)/h1.m_val[ind] * (d0T*h1*sol).m_val[j];
            // mesh.setEdgeFlag(ind, 3);   
        }
        updateComplete(ind, marks, updates);
        // savePicture();
    }
    Text res;
    Text freq;
    Text freq2;
    for(uint i=0; i<mesh.getNodeSize(); i++){
        Vector3 pos = mesh.getNodePosition3(i);
        if(abs(pos.z-tmax)<1e-5){
            Vector3 interp = getField(sol.m_val, i);
            res << pos.x << " " << pos.y << " " << interp.x << " " << interp.y << " " <<  interp.z << "\n";
        }
        if(abs(pos.x-PI/2) + abs(pos.y-PI/2)<0.1){
            Vector3 interp = getField(sol.m_val, i);
            freq << interp.z << " " << pos.z << "\n";
        }
        if(pos.x + abs(pos.y-PI/2)<0.1){
            Vector3 interp = getField(sol.m_val, i);
            freq2 << interp.x << " " << pos.z << "\n";
        }
    }
    freq.save("freq.txt");
    freq2.save("freq2.txt");
    res.save("tsmsh2.txt");
}