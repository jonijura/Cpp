#include "../../GFD/Types/Random.hpp"
#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Discrete/Dec.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include <filesystem>
#include <iostream>
#include <queue>
using namespace std;
using namespace gfd;

PartMesh mesh(0,1,2);

/**
 * interpolate discrete 1-forms on mesh vertices
*/
Vector2 getField(const Buffer<double> &val, const uint node) {
	const Buffer<uint> par = mesh.getNodeEdges(node);
	Matrix2 A(0,0,0,0);
	Vector2 b(0,0);
	for(uint i=0; i<par.size(); i++) {
		const Vector2 v = mesh.getEdgeVector2(par[i]);
		A += v.outerProduct();
		b += v * val[par[i]];
	}
	return A.inverse() * b;
}

void savePicture( string s = "kuva.bmp"){
    Picture pc(500,500);
    pc.fillColor(Vector4(1,1,1,0));
    MeshDrawer md;
    md.initPicture(&pc);
    md.initPosition(Vector4(PI/2,PI/2,1,0), Vector4(PI/2,PI/2,0,0), 
                    Vector4(0.3,0,0,0), Vector4(0,0.3,0,0));
    md.drawDualEdges(mesh, Vector3(1,0,0));
    md.drawPrimalEdges(mesh, Vector3(0,0,1));
    pc.save(s, true);
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
                }
            }
        }
    }

    Buffer<uint> nodes = mesh.getEdgeNodes(edg);
    for(auto a : nodes){
        if(mesh.getNodePosition2(a).y>PI-1e-7){//avoid enforcing von neumann at the end of the mesh
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
                }
            }
        }
    }
}

int main() {
    BuilderMesh bm, bm2;
    // bm.createTriangleGrid(Vector2(0,0), Vector2(PI,PI), 0.2, true);
    // bm.createTriangleGrid(Vector2(0,0), Vector2(2*PI,PI/2-.1), 0.3, true);
    // bm2.createTriangleGrid(Vector2(0,PI/2+.1), Vector2(2*PI,PI), 0.2, true);
    // bm.insertMesh(bm2);
    int n = 8;
    double r = 0.7;
    double a = PI/((1-pow(r,n))/(1-r));
    double c = 0.9;
    double sz = a;
    vector<double> vals(n+1); vals[0]=0;
    vector<double> ts(n+1); ts[0]=c * sz;
    for(uint i=1; i<vals.size(); i++){
        vals[i]=vals[i-1]+sz;
        sz = sz*r;
        ts[i]=c*sz;
    }
    auto cmp = [](pair<double, uint> left, pair<double, uint> right) { return left.first < right.first; };
    priority_queue<pair<double, uint>, vector<pair<double, uint>>, decltype(cmp)> pq(cmp);
    for(uint i=0; i<vals.size(); i++)
        pq.push({0,i});
    while(pq.size()){
        auto a = pq.top();
        // cout << a.first << " " << a.second << endl;
        pq.pop();
        bm.insertNode(Vector4(vals[a.second], a.first,0,0),0,0,false);
        if(a.first!=PI){
            double next = min(a.first+ts[a.second], PI);
            pq.push({next,a.second});
        }
    }
    // bm.fillBoundaryFlags(1);
    bm.setMetric(SymMatrix4(1,0,-1,0,0,0,0,0,0,0));
    // bm.transform(Matrix4(0,1,0,0, 1,0,0,0, 0,0,1,0, 0,0,0,1));
    mesh.swap(bm);
    
    savePicture();
    //status.first: 1,2,3 = node, face, solved
    //second: facenum/nodenum
    Buffer<pair<uint, uint>> marks(mesh.getEdgeSize(), {0,0});
    queue<uint> updates;
    Column<double> sol(mesh.getEdgeSize(), 0.0);
    //edges at the bottom row
    Buffer<uint> initialValues;
    for(uint i = 0; i<mesh.getEdgeSize(); i++){
        if(mesh.getEdgePosition2(i).y == 0){
            initialValues.push_back(i);
            marks[i].first = 3;
            auto n = mesh.getEdgeNodes(i);
            sol.m_val[i] = cos(mesh.getNodePosition1(n[0]))-cos(mesh.getNodePosition1(n[1]));
        }
    }
    for(auto a : initialValues)
        updateComplete(a, marks, updates);
    
    Dec dec(mesh, 0, mesh.getDimension());
    Derivative d1, d0T;
    d1 = dec.integrateDerivative(fg_prim1,d1);
    d0T = dec.integrateDerivative(fg_dual1,d0T);//wont work for 3 dimensions?
    //integrate hodge is too difficult to use, this accounts for negative hodge-elements
    Buffer<double> h1val(mesh.getEdgeSize());
    for(uint i = 0 ; i<mesh.getEdgeSize(); i++)
        h1val[i] = mesh.getEdgeHodge(i);
    Diagonal<double> h1(h1val, 0.0);
    //timestepping
    while(updates.size()!=0){
        uint ind = updates.front();
        uint j = marks[ind].second;
        updates.pop();
        if(marks[ind].first==2){
            marks[ind].first=3;
            sol.m_val[ind] = -(int)mesh.getFaceIncidence(j, ind)*(d1*sol).m_val[j];
        }
        else if(marks[ind].first==1){
            marks[ind].first=3;
            sol.m_val[ind] = -(int)mesh.getEdgeIncidence(ind, j)/h1.m_val[ind] * (d0T*h1*sol).m_val[j];
        }
        updateComplete(ind, marks, updates);
    }

    Text res;
    for(uint i=0; i<mesh.getNodeSize(); i++){
        Vector2 pos = mesh.getNodePosition2(i);
        res << pos.x << " " << pos.y << " " << getField(sol.m_val, i).x << " " << getField(sol.m_val, i).y << "\n";
    }
    res.save("tsmsh.txt");
}