#include "../../GFD/Types/Random.hpp"
#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Discrete/Dec.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include <filesystem>
#include <iostream>
#include <queue>
using namespace std;
using namespace gfd;

// Solving pde system
//     de=0
//     d*e=0
// explicitly on a 2+1-dimensional space-time mesh. e is located at the
// mixed space-time edges and represented as a discrete one form.
// This system corresponds to the wave equation and has a solution
//     e = 1/sqrt(2)cos(x)sin(y)sin(sqrt(2)t)dx 
//       + 1/sqrt(2)sin(x)cos(y)sin(sqrt(2)t)dy
//       + sin(x)sin(y)cos(sqrt(2)t)dt
//       = d(1/sqrt(2)sin(x)sin(y)sin(sqrt(2)t))
// as a boundary condition e is set to 0 on all boundary edges except t=T_MAX
// this is the same as integrating the solution over these edges. As an initial value
// and to achieve an unique solution, e is also initialized by integration for all edges
// with one node at t=0.
// 
// The timestepping uses discrete versions
// of de (when the value of all but one boundary edge of a face is known) and
// d*e (when the value of all but one edge linked to a vertice is known / dual face is missing one edge)
// to advance the solution one edge at a time.
//
// there are different mesh options to modify and choose from at the start of the main function

PartMesh mesh(0,1,3);
// BuilderMesh tmp;
const uint BOUNDARYFLAG = 1;
const double T_MAX = PI;
const double XY_MAX = PI;
const double SOLUTION_REFERENCE_TIME = PI;


Vector3 interpolateOneForm(const Buffer<double> &values_on_edges, const uint node, Buffer<double> &h1) {
	const Buffer<uint> edge_nodes = mesh.getNodeEdges(node);
	Matrix3 A(0,0,0,0,0,0,0,0,0);
	Vector3 b(0,0,0);
	for(uint i=0; i<edge_nodes.size(); i++) {
		const Vector3 v = mesh.getEdgeVector3(edge_nodes[i]);
		A += h1[edge_nodes[i]]*v.outerProduct();
		b += h1[edge_nodes[i]]*v * values_on_edges[edge_nodes[i]];
	}
	return A.inverse() * b;
}

void savePicture(){
    Picture pc(500,500);
    pc.fillColor(Vector4(1,1,1,0));
    MeshDrawer md;
    md.initPicture(&pc);
    Vector4 camera_position = Vector4(PI/2.0, -3, PI/2+2,0);
    Vector4 lookDir = Vector4(PI/2.0, PI/2.0, PI/2.0-1.5, 0);
    md.initPosition(camera_position, lookDir, 
                    Vector4(.2,0,0,0), Vector4(0,.02,.2,0));
    md.drawPrimalEdges(mesh,Vector3(1,0,0), UintSet(BOUNDARYFLAG));
    md.drawPrimalEdges(mesh,Vector3(0,1,0), UintSet(3));
    md.drawPrimalEdges(mesh,Vector3(0,0,1), UintSet(4));
    md.drawPrimalEdges(mesh,Vector3(0,1,1), UintSet(5));
    md.drawPrimalEdges(mesh,Vector3(1,0,1), UintSet(6));

    //leikkauskuvan piirto
    // for(int i=tmp.getNodeSize(); i>=0; i--){
    //     auto a = tmp.getNodePosition3(i);
    //     if(a.y<a.z-.5)
    //         tmp.removeNode(i);
    // }
    // tmp.fillBoundaryFlags(4);
    // md.drawBoundaryFaces(tmp, Vector3(0,1,0), UintSet(4));
    // md.drawPrimalEdges(tmp,Vector3(0,0,1), UintSet(4), 0.01);

    pc.save("2+1mesh.bmp", true);
} 

/**
 * update the list and detail of edges that can be solved via dF=0 or d*F=0
*/
void addSolvableEdges(uint updatedEdge, Buffer<pair<uint, uint>> &marks, queue<uint> &solvableEdges){
    Buffer<uint> faces = mesh.getEdgeFaces(updatedEdge);
    for(auto face : faces){
        Buffer<uint> faceEdges = mesh.getFaceEdges(face);
        uint missing = 0;
        for(auto edge : faceEdges)
            if(marks[edge].first!=3)missing++;
        if(missing==1){
            for(auto edge : faceEdges){
                if(marks[edge].first==0){
                    solvableEdges.push(edge);
                    marks[edge]={2,face};
                    mesh.setEdgeFlag(edge,4);
                }
            }
        }
    }

    Buffer<uint> nodes = mesh.getEdgeNodes(updatedEdge);
    for(auto node : nodes){
        //using d*F=0 on the boundary assumes homogeneous von neumann boundary conditions.
        if(mesh.getNodeFlag(node)==BOUNDARYFLAG){
            continue;
        }
        Buffer<uint> nodeEdges = mesh.getNodeEdges(node);
        uint missing = 0;
        for(auto edge : nodeEdges)
            if(marks[edge].first!=3)missing++;
        if(missing==1){
            for(auto edge : nodeEdges){
                if(marks[edge].first==0){
                    solvableEdges.push(edge);
                    marks[edge]={1,node};
                    mesh.setEdgeFlag(edge,5);
                }
            }
        }
    }
}

void createReqularMesh(double dxy, double CFL_limit=0.5){
    BuilderMesh bm;
    bm.createGrid(Vector4(0,0,0,0), Vector4(PI,PI,T_MAX,0), Vector4(dxy,dxy,dxy*CFL_limit,1));
    bm.setMetric(SymMatrix4(1,0,1,0,0,-1,0,0,0,0));
    bm.fillBoundaryFlags(BOUNDARYFLAG);
    mesh.swap(bm);
}

void createPartlyRefinedMesh(double dxy_coarse=0.3, double dxy_fine=0.2, double CFL_limit=0.5){
    BuilderMesh bm ,bm2;
    double meshSpacing = (dxy_coarse+dxy_fine)/2;
    bm.createTriangleGrid(Vector2(0,0), Vector2(XY_MAX/2-meshSpacing/2,XY_MAX), dxy_coarse, true);
    bm.stretchLinear(Vector4(0,0,T_MAX,0),ceil(T_MAX/(CFL_limit*dxy_coarse)));
    bm2.createTriangleGrid(Vector2(XY_MAX/2+meshSpacing/2,0), Vector2(XY_MAX,XY_MAX), dxy_fine, true);
    bm2.stretchLinear(Vector4(0,0,T_MAX,0),ceil(T_MAX/(CFL_limit*dxy_fine)));
    bm.insertMesh(bm2);
    bm.setMetric(SymMatrix4(1,0,1,0,0,-1,0,0,0,0));
    bm.fillBoundaryFlags(1);
    Text statsBound;
    bm.writeStatistics(statsBound, UintSet(1));
    statsBound.save("statsBoundary.txt");
    // tmp.swap(bm);
    mesh.swap(bm);
    cout << "mesh done\n";
}

void createPartlyRefinedMesh2(double dxy = 0.5){
    BuilderMesh bm;
    bm.createTriangleGrid(Vector2(0,0), Vector2(XY_MAX, XY_MAX), dxy, true);
    //refine by inserting nodes at circumcenters
    int near = 0;
    Buffer<Vector4> refiningNodes;
    for(uint i=0; i<bm.getFaceSize(); i++){
        Vector4 pos = bm.getFacePosition(i);
        if(pos.x>=0.25*XY_MAX && pos.y>= 0.25*XY_MAX &&pos.x<=0.75*XY_MAX && pos.y<=0.75*XY_MAX)
            refiningNodes.push_back(pos);
    }
    for(uint i=0; i< refiningNodes.size(); i++)
        near = bm.insertNode(refiningNodes[i],0.0,near,true);


    bm.setMetric(SymMatrix4(1,0,1,0,0,-1,0,0,0,0));
    bm.fillBoundaryFlags(1);

    mesh.swap(bm);
    savePicture();
    //build spacetime mesh with tent pitcher
    cout << "mesh not implemented\n";
}

int main() {
    // createReqularMesh(0.3, 1);
    createPartlyRefinedMesh(0.4,0.2);//the results in the article draft used createPartlyRefinedMesh(0.4,0.2);
    // createPartlyRefinedMesh2();

    //status.first: 1,2,3 = node, face, solved
    //second: facenum/nodenum
    Buffer<pair<uint, uint>> marks(mesh.getEdgeSize(), {0,0});
    queue<uint> solvableEdges;
    Column<double> solution(mesh.getEdgeSize(), 0.0);
    //edges at the bottom row
    Buffer<uint> initialValues;
    for(uint edge = 0; edge<mesh.getEdgeSize(); edge++){
        auto edge_cc_position = mesh.getEdgePosition3(edge);
        auto edge_nodes = mesh.getEdgeNodes(edge);
        //lock value at boundary nodes to 0
        if(edge_cc_position.x == 0 || edge_cc_position.y==0 || edge_cc_position.z==0 ||
             edge_cc_position.x==XY_MAX || edge_cc_position.y==XY_MAX){
            initialValues.push_back(edge);
            marks[edge].first = 3;
            // mesh.setEdgeFlag(i,3);
        }
        //Set initial values on the first row of non boundary edges
        else if((mesh.getEdgeFlag(edge)!=1) &&
                (mesh.getNodePosition(edge_nodes[0]).z==0 || mesh.getNodePosition(edge_nodes[1]).z==0)){
            marks[edge].first = 3;
            // mesh.setEdgeFlag(i,4); 
            auto edge_end = mesh.getNodePosition3(edge_nodes[0]);
            auto edge_start = mesh.getNodePosition3(edge_nodes[1]);
            solution.m_val[edge] = (int)mesh.getEdgeIncidence(edge,edge_nodes[0])*sin(edge_end.x)*sin(edge_end.y)*sin(sqrt(2.0)*edge_end.z) / sqrt(2.0)
                        +(int)mesh.getEdgeIncidence(edge,edge_nodes[1])*sin(edge_start.x)*sin(edge_start.y)*sin(sqrt(2.0)*edge_start.z) / sqrt(2.0);
        }
    }
    for(auto edge : initialValues)
        addSolvableEdges(edge, marks, solvableEdges);
    savePicture();
    if(solvableEdges.size() == 0)throw std::invalid_argument( "could not find any edges to solve, check boundary and initial conditions" );
    
    Dec dec(mesh, 0, mesh.getDimension());
    Derivative d1, d0T;
    d1 = dec.integrateDerivative(fg_prim1,d1);
    d0T = dec.integrateDerivative(fg_prim0,d0T);
    d0T = d0T.setTranspose(d0T);
    //integrate hodge is too difficult to use, this seems to account for negative hodge-elements
    Buffer<double> h1_values(mesh.getEdgeSize());
    for(uint i = 0 ; i<mesh.getEdgeSize(); i++){
        h1_values[i] = mesh.getEdgeHodge(i);
    }
    Diagonal<double> h1(h1_values, 0.0);

    uint countNegative = 0;
    for(double value : h1_values)
        if(value<0)countNegative++;
    if(countNegative == 0)throw std::invalid_argument( "no negative hodge elements!" );

    //timestepping
    while(solvableEdges.size()!=0){
        uint ind = solvableEdges.front();
        uint j = marks[ind].second;
        solvableEdges.pop();
        if(marks[ind].first==2){
            solution.m_val[ind] = -(int)mesh.getFaceIncidence(j, ind)*(d1*solution).m_val[j];
            mesh.setEdgeFlag(ind,5);
        }
        else if(marks[ind].first==1){
            solution.m_val[ind] = -(int)mesh.getEdgeIncidence(ind, j)/h1.m_val[ind] * (d0T*h1*solution).m_val[j];
            mesh.setEdgeFlag(ind, 3);   
        }
        marks[ind].first=3;
        addSolvableEdges(ind, marks, solvableEdges);
    }

    Text interpolated_result;
    for(uint i=0; i<mesh.getNodeSize(); i++){
        Vector3 pos = mesh.getNodePosition3(i);
        if(abs(pos.z-SOLUTION_REFERENCE_TIME)<1e-5){
            Vector3 interp = interpolateOneForm(solution.m_val, i, h1.m_val);
            interpolated_result << pos.x << " " << pos.y << " " << interp.x << " " << interp.y << " " <<  interp.z << "\n";
        }
    }
    interpolated_result.save("2+1interpolated.txt");

    Text no_interpolation;
    for(uint i=0; i<mesh.getEdgeSize(); i++){
        auto nodes = mesh.getEdgeNodes(i);
        Vector3 pos1 = mesh.getNodePosition3(nodes[0]);
        Vector3 pos2 = mesh.getNodePosition3(nodes[1]);
        no_interpolation << pos1.x << " " << pos1.y << " " << pos1.z << " "
                         << pos2.x << " " << pos2.y << " " << pos2.z << " "
                         << solution.m_val[i] << "\n";
    }
    no_interpolation.save("2+1nointerpolation.txt");

    //find out the interpolation error, for this initialize exact solution and interpolate it:
    Text interpolationError;
    for(uint edge=0; edge<mesh.getEdgeSize(); edge++){
        auto edge_nodes = mesh.getEdgeNodes(edge);
        auto edge_end = mesh.getNodePosition3(edge_nodes[0]);
        auto edge_start = mesh.getNodePosition3(edge_nodes[1]);
        solution.m_val[edge] = (int)mesh.getEdgeIncidence(edge,edge_nodes[0])*sin(edge_end.x)*sin(edge_end.y)*sin(sqrt(2.0)*edge_end.z) / sqrt(2.0)
            +(int)mesh.getEdgeIncidence(edge,edge_nodes[1])*sin(edge_start.x)*sin(edge_start.y)*sin(sqrt(2.0)*edge_start.z) / sqrt(2.0);
    }
    for(uint i=0; i<mesh.getNodeSize(); i++){
        Vector3 pos = mesh.getNodePosition3(i);
        if(abs(pos.z-SOLUTION_REFERENCE_TIME)<1e-5){
            Vector3 interp = interpolateOneForm(solution.m_val, i, h1.m_val);
            interpolationError << pos.x << " " << pos.y << " " << interp.x << " " << interp.y << " " <<  interp.z << "\n";
        }
    }
    interpolationError.save("2+1exact_solution_interpolated.txt");

    Text stats;
    mesh.writeStatistics(stats);
    stats.save("stats.txt");
    cout << initialValues.size();
}