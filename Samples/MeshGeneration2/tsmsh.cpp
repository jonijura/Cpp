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
// explicitly on a 1+1-dimensional space-time mesh. e is located at the
// mixed space-time edges and represented as a discrete one form.
// This system corresponds to the wave equation and has a solution
//     e = sin(x)cos(t)dx + cos(x)sin(t)dt = d(-cos(x)cos(t))
// e is initialized at t=0, x=[0,X_MAX], by integrating the smooth solution over corresponding edges.
// The analytical value of this integral is simply cos(x2)cos(t2)-cos(x1)cos(t1) due to Stokes
// 
// The solution is chosen such that it conforms with the natural von neumann boundary e|_bound = 0
// condition, which comes from incoplete dual forms at the boundary. No additional
// consideration at the boundary is needed. The timestepping uses discrete versions
// of de (when the value of all but one boundary edge of a face is known) and
// d*e (when the value of all but one edge linked to a vertice is known / dual face is missing one edge)
// to advance the solution one edge at a time.
//
// there are three mesh options to modify and choose from at the start of the main function


PartMesh mesh(0,1,2);
const double T_MAX = PI;
const double X_MAX = PI;
/**
 * interpolate discrete 1-forms on mesh vertices
*/
Vector2 interpolateOneform(const Buffer<double> &values_on_edges, const uint node) {
	const Buffer<uint> node_edges = mesh.getNodeEdges(node);
	Matrix2 LHS(0,0,0,0);
	Vector2 RHS(0,0);
	for(uint i=0; i<node_edges.size(); i++) {
		const Vector2 v = mesh.getEdgeVector2(node_edges[i]);
		LHS += v.outerProduct();
		RHS += v * values_on_edges[node_edges[i]];
	}
	return LHS.inverse() * RHS;
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
void addSolvableEdges(uint updatedEdge, Buffer<pair<uint, uint>> &marks, queue<uint> &solvableEdges){
    Buffer<uint> faces = mesh.getEdgeFaces(updatedEdge);
    for(auto face : faces){
        Buffer<uint> faceEdges = mesh.getFaceEdges(face);
        uint unsolvedEdgesCount = 0;
        for(auto edge : faceEdges)
            if(marks[edge].first!=3)unsolvedEdgesCount++;
        if(unsolvedEdgesCount==1){
            for(auto edge : faceEdges){
                if(marks[edge].first==0){
                    solvableEdges.push(edge);
                    marks[edge]={2,face};
                }
            }
        }
    }

    Buffer<uint> nodes = mesh.getEdgeNodes(updatedEdge);
    for(auto node : nodes){
        //avoid enforcing von neumann boundary condition at the end of the mesh
        if(mesh.getNodePosition2(node).y>T_MAX-1e-7){
            continue;
        }
        Buffer<uint> nodeEdges = mesh.getNodeEdges(node);
        uint unsolvedEdgesCount = 0;
        for(auto edge : nodeEdges)
            if(marks[edge].first!=3)unsolvedEdgesCount++;
        if(unsolvedEdgesCount==1){
            for(auto edge : nodeEdges){
                if(marks[edge].first==0){
                    solvableEdges.push(edge);
                    marks[edge]={1,node};
                }
            }
        }
    }
}

void createReqularMesh(double dx = 0.2){
    BuilderMesh bm;
    bm.createTriangleGrid(Vector2(0,0), Vector2(T_MAX,X_MAX), dx, true);
    bm.setMetric(SymMatrix4(1,0,-1,0,0,0,0,0,0,0));
    bm.transform(Matrix4(0,1,0,0, 1,0,0,0, 0,0,1,0, 0,0,0,1));
    mesh.swap(bm);
}

void createOnceRefinedMesh( double coarsedx = 0.3, double finedx = 0.2){
    BuilderMesh bm, bm2;
    double meshSpacing = (coarsedx+finedx)/2;

    bm.createTriangleGrid(Vector2(0,0), Vector2(T_MAX, X_MAX/2-meshSpacing/2), coarsedx, true);
    bm2.createTriangleGrid(Vector2(0,X_MAX/2+meshSpacing/2), Vector2(T_MAX, X_MAX), finedx, true);
    bm.insertMesh(bm2);
    bm.setMetric(SymMatrix4(1,0,-1,0,0,0,0,0,0,0));
    bm.transform(Matrix4(0,1,0,0, 1,0,0,0, 0,0,1,0, 0,0,0,1));
    mesh.swap(bm);
}

void createContinouslyRefinedMesh(  int discretizationLevel = 12, 
                                    double refiningFactor = 0.8, 
                                    double CFL_limit = 0.9){

    BuilderMesh bm;
    double lengthModifier = X_MAX/((1-pow(refiningFactor,discretizationLevel))/(1-refiningFactor));
    cout << to_string(lengthModifier) << "\n";
    double currentLength = 1.0 * lengthModifier;
    vector<double> position_x(discretizationLevel+1); 
    position_x[0]=0;
    vector<double> timestep_size(discretizationLevel+1); 
    timestep_size[0]=T_MAX/ceil(T_MAX/(CFL_limit*currentLength));
    for(uint i=1; i<position_x.size(); i++){
        position_x[i]=position_x[i-1]+currentLength;
        currentLength = currentLength*refiningFactor;
        timestep_size[i]=T_MAX/ceil(T_MAX/(CFL_limit*currentLength));
    }
    auto comparator = [](pair<double, uint> left, pair<double, uint> right) { return left.first < right.first; };
    priority_queue<pair<double, uint>, vector<pair<double, uint>>, decltype(comparator)> nodePositions(comparator);
    for(uint i=0; i<position_x.size(); i++)
        nodePositions.push({0,i});
    while(nodePositions.size()){
        auto nextTimeValue = nodePositions.top();
        nodePositions.pop();
        bm.insertNode(Vector4(position_x[nextTimeValue.second], nextTimeValue.first,0,0),0,0,false);
        if(nextTimeValue.first!=T_MAX){
            double next = min(nextTimeValue.first+timestep_size[nextTimeValue.second], T_MAX);
            nodePositions.push({next,nextTimeValue.second});
        }
    }
    bm.setMetric(SymMatrix4(1,0,-1,0,0,0,0,0,0,0));
    mesh.swap(bm);
}

int main() {
    // choose mesh and refining parameters
    // createReqularMesh(0.2);
    // createOnceRefinedMesh(0.3, 0.2); //this mesh was used in some of the 1+1 simulations

    //doing some convergence tests, check that the required folder structure exists
    for(uint i=3; i<4; i++){   
        double discretization_level = i*i;
        double finest_coarsest_ratio = 1.0/3;
        double refinement_factor = exp(log(finest_coarsest_ratio)/discretization_level);
    createContinouslyRefinedMesh(discretization_level, refinement_factor, 0.9);
    // createReqularMesh(PI/discretization_level);

    savePicture("1+1tsmsh.bmp");
    //marks.first: 1,2,3 = node rule, face rule, solved
    //marks.second: node/face index to be used for solving
    Buffer<pair<uint, uint>> marks(mesh.getEdgeSize(), {0,0});
    queue<uint> solvableEdges;
    Column<double> solution(mesh.getEdgeSize(), 0.0);
    //edges at the bottom row
    Buffer<uint> initialValues;
    for(uint i = 0; i<mesh.getEdgeSize(); i++){
        if(mesh.getEdgePosition2(i).y == 0){
            initialValues.push_back(i);
            marks[i].first = 3;
            auto nodes = mesh.getEdgeNodes(i);
            solution.m_val[i] = cos(mesh.getNodePosition1(nodes[0]))-cos(mesh.getNodePosition1(nodes[1]));
        }
    }
    for(auto edge : initialValues)
        addSolvableEdges(edge, marks, solvableEdges);
    
    Dec dec(mesh, 0, mesh.getDimension());
    Derivative d1, d0T;
    d1 = dec.integrateDerivative(fg_prim1, d1);
    d0T = dec.integrateDerivative(fg_dual1, d0T);//wont work for 3 dimensions?
    //construction of hodge, negative sign for timelike edges comes automatically
    Buffer<double> h1_values(mesh.getEdgeSize());
    for(uint i = 0 ; i<mesh.getEdgeSize(); i++)
        h1_values[i] = mesh.getEdgeHodge(i);
    Diagonal<double> h1(h1_values, 0.0);
    //timestepping, replacing matrix multiplications should speed up the program a lot
    while(solvableEdges.size()!=0){
        uint i = solvableEdges.front();
        uint j = marks[i].second;
        solvableEdges.pop();
        if(marks[i].first==2){
            solution.m_val[i] = -(int)mesh.getFaceIncidence(j, i)*(d1*solution).m_val[j];
        }
        else if(marks[i].first==1){
            solution.m_val[i] = -(int)mesh.getEdgeIncidence(i, j)/h1.m_val[i] * (d0T*h1*solution).m_val[j];
        }
        marks[i].first=3;
        addSolvableEdges(i, marks, solvableEdges);
    }

    Text resultsInterpolated;
    for(uint i=0; i<mesh.getNodeSize(); i++){
        Vector2 pos = mesh.getNodePosition2(i); 
        resultsInterpolated << pos.x << " " << pos.y << " " << interpolateOneform(solution.m_val, i).x << " " << interpolateOneform(solution.m_val, i).y << "\n";
    }
    resultsInterpolated.save("conv\\" + to_string((int)round(discretization_level)) + "_interpolated.txt");
    Text resultsForm;
    for(uint i=0; i<mesh.getEdgeSize(); i++){
        auto nodes = mesh.getEdgeNodes(i);
        Vector2 pos = mesh.getNodePosition2(nodes[0]),  pos2 = mesh.getNodePosition2(nodes[1]);
        resultsForm << pos.x << " " << pos.y << " " <<  pos2.x << " " << pos2.y << " "  << solution.m_val[i] << "\n";
    }
    resultsForm.save("conv\\" + to_string((int)round(discretization_level)) + "_1+1form.txt");
    }
}