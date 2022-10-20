#include "../../GFD/Types/Random.hpp"
#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include <filesystem>
#include <iostream>
using namespace std;
using namespace gfd;

void saveVerticesAndFaces(Mesh &mesh){
    Text vert;
    for (uint i = 0; i < mesh.getNodeSize(); i++)
    {
        vert << mesh.getNodePosition2(i).x << " " << mesh.getNodePosition2(i).y << "\n";
    }
    Text tria;
    for (uint i=0; i< mesh.getFaceSize(); i++)
    {
        auto nodes = mesh.getFaceNodes(i);
        for(uint i=0; i<nodes.size(); i++)
            tria << nodes[i] << " ";
        tria << "\n";
    }
    tria.save("build\\tria.txt");
    vert.save("build\\vert.txt");
}

void savePicture(Mesh &mesh, string &s){
    Picture pc(500,500);
    MeshDrawer md;
    md.initPicture(&pc);
    md.initPosition(Vector4(0.5,0.5,1,0), Vector4(0.5,0.5,0,0), Vector4(1,0,0,0), Vector4(0,1,0,0));
    md.drawPrimalEdges(mesh);
    md.drawDualEdges(mesh, Vector3(1,0,0));
    pc.save(s, true);
}

void refine(BuilderMesh &mesh){
    int near = 0;
    uint n = mesh.getFaceSize();
    Buffer<Vector4> refiningNodes;
    for(uint i=0; i<n; i++){
        Vector4 pos = mesh.getFacePosition(i);
        if(pos.x>=0.25 && pos.y>= 0.25 &&pos.x<=0.75 && pos.y<=0.75)
            refiningNodes.push_back(pos);
    }
    for(uint i=0; i< refiningNodes.size(); i++){
        near = mesh.insertNode(refiningNodes[i],0.0,near,true);
    }
    return;
    //collapse edges
    Buffer<uint> edgesToRemove;
    for(int i=mesh.getEdgeSize()-1; i>=0 ; i--){
        if(mesh.getEdgeHodge(i)<0.1){
            auto faces = mesh.getEdgeFaces(i);
            if(faces.size()<2)continue;
            auto edges = mesh.getFaceEdges(faces[0]).getUnion(mesh.getFaceEdges(faces[1]));
            edges.eraseFirst(i);
            mesh.addFace(edges);
            edgesToRemove.push_back(i);
        }
    }
    for(uint i=0; i<edgesToRemove.size(); i++)
        mesh.removeEdge(edgesToRemove[i]);
}

int main() {
    Random rnd(14);
    BuilderMesh mesh(2);
    mesh.createTriangleGrid(Vector2(0,0), Vector2(1,1), 0.1, true);
    refine(mesh);
    // BuilderMesh inside(2);
    // inside.createTriangleGrid(Vector2(0.25,0.25), Vector2(0.75,0.75), 0.05);
    // BuilderMesh inside2(2);
    // inside2.createTriangleGrid(Vector2(0.375,0.375), Vector2(0.625,0.625), 0.025);
    // inside.insertMesh(inside2);
    // mesh.insertMesh(inside);
    saveVerticesAndFaces(mesh);
    string s = "build\\kuva2.bmp";
    savePicture(mesh, s);
    mesh.fillBoundaryFlags(5);
    // mesh.improveNodeByHodge(0);
    // mesh.optimizeNodes(UintSet(0));
    string s2 = "build\\kuva3.bmp";
    savePicture(mesh, s2); 
}