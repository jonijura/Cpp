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

void savePicture(){
    Picture pc(500,500);
    MeshDrawer md;
    md.initPicture(&pc);
    Vector4 campos = Vector4(.5, -1, .5,0); 
    md.initPosition(campos, Vector4(.5, .5, .5,0), 
                    Vector4(.5,0,0,0), Vector4(0,0,.5,0));
    md.drawPrimalEdges(mesh,Vector3(1,0,0), UintSet(1));
    md.drawPrimalEdges(mesh,Vector3(0,1,0), UintSet(3));
    md.drawPrimalEdges(mesh,Vector3(0,0,1), UintSet(4));
    md.drawPrimalEdges(mesh,Vector3(0,1,1), UintSet(5));
    md.drawPrimalEdges(mesh,Vector3(1,0,1), UintSet(6));
    md.drawPrimalEdges(mesh,Vector3(1,1,1));//draw all edges without flags last
    pc.save("ts3.bmp", true);
}

/**
 * update the list and detail of edges that can be solved via dF=0 or d*F=0
*/
void addSolvableEdges(uint edg, Buffer<pair<uint, uint>> &marks, queue<uint> &updates){
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
                    mesh.setEdgeFlag(b,5);
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
                    mesh.setEdgeFlag(b,6);
                }
            }
        }
    }
}


int main() {
    BuilderMesh bm;
    bm.createTriangleGrid(Vector2(0,0), Vector2(0.45,1), 0.3, true);
    bm.stretchLinear(Vector4(0,0,0.6,0),4);
    BuilderMesh bm2;
    bm2.createTriangleGrid(Vector2(0.55,0), Vector2(1,1), 0.2, true);
    bm2.stretchLinear(Vector4(0,0,0.6,0),6);
    bm.insertMesh(bm2);
    bm.fillBoundaryFlags(1);
    mesh.swap(bm);
    savePicture();
    Buffer<pair<uint, uint>> marks(mesh.getEdgeSize(), {0,0});
    Buffer<uint> initialValues;
    for(uint i = 0; i<mesh.getEdgeSize(); i++){
        auto b = mesh.getEdgeNodes(i);
        if((mesh.getNodeFlag(b[0]) + mesh.getNodeFlag(b[1])>=1) &&
            (mesh.getEdgeFlag(i)!=1) &&
            (mesh.getNodePosition(b[0]).z==0 || mesh.getNodePosition(b[1]).z==0) ){
            marks[i].first = 3;
            mesh.setEdgeFlag(i,4);
            initialValues.push_back(i);
        }
    }
    savePicture();
    for(uint i = 0; i<mesh.getEdgeSize(); i++){
        auto a = mesh.getEdgePosition3(i);
        if(a.z == 0 || a.x==0 || a.y==0 || a.x==1 || a.y==1){
            marks[i].first = 3;
            // mesh.setEdgeFlag(i,3);
        }
    }
    queue<uint> updates;
    for(auto a : initialValues){
        addSolvableEdges(a, marks, updates);
    }
    savePicture();
    while(updates.size()!=0){
        uint ind = updates.front();
        updates.pop();
        if(marks[ind].first==2){
            marks[ind].first=3;
            mesh.setEdgeFlag(ind, 4);
        }
        else if(marks[ind].first==1){
            marks[ind].first=3;
            mesh.setEdgeFlag(ind, 4);   
        }
        addSolvableEdges(ind, marks, updates);
    }
    savePicture();
}