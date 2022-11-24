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
    Vector4 campos = Vector4(.4, -1, .6,0); 
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


int main() {
    BuilderMesh bm;
    bm.createA15Grid(Vector3(0,0,0), Vector3(1,1,1), 1);
    // bm.createGrid(Vector4(0,0,0,0), Vector4(1,1,1,0), 0.5*Vector4(1,1,1,1));
    // bm.createGrid(Vector4(0,0,0,0), Vector4(1,1,1,0), Vector4(0.5,0.5,0.5,1));
    bm.fillBoundaryFlags(1);
    mesh.swap(bm);
    savePicture();

}