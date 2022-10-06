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

void savePicture(Mesh &mesh, string s){
    Picture pc(500,500);
    MeshDrawer md;
    md.initPicture(&pc);
    md.initPosition(Vector4(0.5,0.5,1,0), Vector4(0.5,0.5,0,0), Vector4(1,0,0,0), Vector4(0,1,0,0));
    md.drawPrimalEdges(mesh);
    pc.save(s, true);
}

int main() {
    Random rnd(14);
    BuilderMesh mesh(2);
    mesh.createTriangleGrid(Vector2(0,0), Vector2(1,1), 0.1, true);
    BuilderMesh inside(2);
    inside.createTriangleGrid(Vector2(0.25,0.25), Vector2(0.75,0.75), 0.05);
    BuilderMesh inside2(2);
    inside2.createTriangleGrid(Vector2(0.375,0.375), Vector2(0.625,0.625), 0.025);
    inside.insertMesh(inside2);
    mesh.insertMesh(inside);
    saveVerticesAndFaces(mesh);
    savePicture(mesh, "build\\kuva1.bmp");
    mesh.fillBodyFlags(1);
    mesh.improveNodeByHodge(0);
    savePicture(mesh, "build\\kuva3.bmp");

}