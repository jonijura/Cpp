#include <iostream>
#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Mesh/PartMesh.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include "../../GFD/Discrete/Dec.hpp"

using namespace std;
using namespace gfd;

const int N=10.0;
const int ts = 50;
const double tp = 2.5*PI/(sqrt(2.0));
const double dt = tp/double(ts);
const Vector4 refp(PI/2,PI/2,0,0);

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
    md.initPosition(Vector4(0.5,0.5,1,0), Vector4(0.5,0.5,0,0), Vector4(0.9,0,0,0), Vector4(0,0.9,0,0));
    md.drawPrimalEdges(mesh);
    md.drawDualEdges(mesh, Vector3(1,0,0));
    pc.save(s, true);
}

void makeMsh(PartMesh& mesh){
    BuilderMesh bm;
    bm.createTriangleGrid(Vector2(0.0,0.0), Vector2(PI,PI),PI/N, true);
    bm.fillBoundaryFlags(1);
    //flag nodes on left and right boundary
    for(uint i=0; i<bm.getNodeSize(); i++)
        if(bm.getNodePosition2(i).x<1e-8 || bm.getNodePosition2(i).x>PI-1e-8)
            bm.setNodeFlag(i,2);
    // bm.addNode(Vector4(0,0,0,0));
    // bm.addNode(Vector4(1,0,0,0));
    // bm.addNode(Vector4(1,1,0,0));
    // bm.addNode(Vector4(0,1,0,0));
    // Buffer<uint> edges;
    // edges.push_back(bm.addEdge(0,1));
    // edges.push_back(bm.addEdge(1,2));
    // edges.push_back(bm.addEdge(2,3));
    // edges.push_back(bm.addEdge(3,0));
    // bm.addFace(edges);

    // string s = "build\\recmsh.bmp";
    // savePicture(bm, s);
    mesh.swap(bm);
}

int main() {
    PartMesh mesh(0,1,2);
    makeMsh(mesh);
    Dec dec(mesh, 0, mesh.getDimension());
    Derivative d0;
    dec.integrateDerivative(fg_prim0, d0);
    Diagonal<double> h1,h0i;
    // dec.integrateHodge(HodgeUnit2, 0, fg_prim1, h1);
    // dec.integrateHodge(HodgeUnit2, 0, fg_dual2, h0i);
    Buffer<double> h1v(mesh.getEdgeSize());
    for(uint i=0; i<mesh.getEdgeSize(); i++)
        h1v[i]=mesh.getEdgeHodge(i);
    Buffer<double> h0iv(mesh.getNodeSize());
    for(uint i=0; i<mesh.getNodeSize(); i++)
        h0iv[i]=mesh.getNodeFlag(i)>0 ? 0 : 1.0/mesh.getNodeHodge(i);
    h1.setFull(h1v);
    h0i.setFull(h0iv);
    Sparse<double> A,B;
    A.setScale(dt,h1*d0);
    B.setScale(dt,-h0i*transpose(d0));
    
    // cout << h1.m_height << endl;
    // cout << d0.m_height << " " << d0.m_width << endl;
    // cout << h0i.m_height << endl;
    // cout << mesh.getNodeSize() << " " << mesh.getEdgeSize() << " " << mesh.getFaceSize() << endl;

    Buffer<double> fv(mesh.getNodeSize());
    Buffer<double> gv(mesh.getEdgeSize());
    for(uint i=0; i<mesh.getNodeSize(); i++){
        Vector2 pos = mesh.getNodePosition2(i);
        fv[i]=sin(pos.x)*sin(pos.y);
    }
    Column<double> f(fv,0);
    Column<double> g(gv,0);
    g.setScale(0.5,A*f);

    Buffer<double> hist(ts, 0);
    uint tp = mesh.findNode(refp, .1);

    for(int i=0; i<ts; i++){
        f+=B*g;
        g+=A*f;
        hist[i]=f.m_val[tp];
    }

    ofstream myfile;
    myfile.open("build\\solutions.txt");
    for(uint i=0; i<ts; i++){
        myfile << hist[i] << "\n";
    }
    myfile.close();
    saveVerticesAndFaces(mesh);
}