/**
 * Ohjelma, jonka avulla voi selvittää läpäisykertoimen fotonikiteen läpi TM aallolle
 */

#include "../../GFD/Mesh/BuilderMesh.hpp"	//for building and modifying a mesh
#include "../../GFD/Mesh/PartMesh.hpp"	//for working with parts of a mesh?
#include "../../GFD/Output/MeshDrawer.hpp"	//for drawing a mesh
#include "../../GFD/Types/Types.hpp"	//constants, shorter names for types, type couple
#include "../../GFD/Types/MpiEasy.hpp"	//for parallel computing
#include "../../GFD/Discrete/Dec.hpp"	//for dec related calculations
#include <iostream>
#include <chrono>	//For recording time

using namespace std;
using namespace gfd;

// modify these constants:
const uint space_steps = 16; // number of space elements per unit length
const uint time_steps = 64; // number of time elements per unit time
const uint T = 200; //simulation time
const uint grid_type = 1; // 2d grid type: 0 = squares, 1 = triangles, 2 = snubsquare
const double fm = 0.5;
const double lev = 3.0/(PI*fm);
const double epsilon = 8.9;
const double r = 0.2;

bool createMesh(PartMesh &mesh) {
	if(mesh.getNumberOfParts() > 1) {
		if(mesh.getPart() == 0) cout << "MPI mesh generation is not implemented." << endl;
		return false;
	}
	BuilderMesh bmesh(mesh.getDimension());	//cpp class instance cmp java: BuildeMesh bmesh = new BuilderMesh(mesh.getDimension());
	const double hx = 1.0 / double(space_steps);
	if(grid_type == 0) bmesh.createGrid(Vector4(0,0,0,0), Vector4(10,1,0,0), hx);
	else if(grid_type == 1) bmesh.createTriangleGrid(Vector2(0,0), Vector2(10,1), hx, true);
	else if(grid_type == 2) bmesh.createSnubSquareGrid(Vector2(0,0), Vector2(10,1), hx);
	bmesh.fillBoundaryFlags(1);	//Requires 2 parameters, second parameter set to default defined in BuilderMesh.hpp
	
	if(mesh.getPart() == 0) mesh.swap(bmesh);	//set the mesh to the mesh build using bmesh
	return true;
}

int main() {
	initMPI();
	auto starttime = chrono::system_clock::now();

	// create mesh
	PartMesh mesh(0, 1, 2);
	if(!createMesh(mesh)) return 0;

	// initialize matrices and vectors
	Dec dec(mesh, 0,  mesh.getDimension());
    Derivative d0; 
    d0=dec.integrateDerivative(fg_prim0, d0);
    Diagonal<double> h1;
    h1=dec.integrateHodge(HodgeUnit2, 0, fg_prim1, h1);//mistä tulee HodgeUnit2??
    Diagonal<double> h0i;
    h0i=dec.integrateHodge(HodgeUnit1, 0, fg_dual0, h0i);

    const double dt = 1.0 / double(time_steps);
    const double dx = 1.0 / double(space_steps);

    Buffer<double> holes(mesh.getNodeSize());
    Buffer<double> s_m(mesh.getNodeSize());
    Buffer<double> source;
	double xn,yn,xm,ym;
    double eps = 1e-5;
    uint tp,tpr;
	for(double i=0;i<mesh.getNodeSize();i++){
        s_m[i]=0;
        xn=mesh.getNodePosition2(i).x;
        yn=mesh.getNodePosition2(i).y;
        if(xn<eps){
            source.push_back(i);
            //cout << xn << ", " << yn << endl;
            s_m[i]=1.0;
            if(yn<eps | yn > 1-eps) s_m[i]=.5;
        }
        if(xn>10-eps){
            s_m[i]=1.0;
            if(yn<eps | yn > 1-eps) s_m[i]=.5;
        }
        if(abs(xn-0.5)+abs(yn-0.5)<0.05){
            tp=i;
            cout << "Tarkastelupiste 1: " << xn << ", " << yn << endl;
        }
        if(abs(xn-9.5)+abs(yn-0.5)<0.05){
            tpr=i;
            cout << "Tarkastelupiste 2: " << xn << ", " << yn << endl;
        } 
		xm=fmod(xn,1.0)-0.5;
		ym=fmod(yn,1.0)-0.5;
		if((xn>2)&(xn<8)&(xm*xm+ym*ym)<(r*r)) holes[i]=1.0/epsilon;
		else holes[i]=1;
	}
	Diagonal<double> P2(holes,0.0);

    Column<double> h(0.0);
	dec.integrateZeroForm(fg_dual1, h);
	Column<double> e(0.0);
	dec.integrateZeroForm(fg_prim0, e);
    Column<double> h2(0.0);
	dec.integrateZeroForm(fg_dual1, h2);
	Column<double> e2(0.0);
	dec.integrateZeroForm(fg_prim0, e2);
    Sparse<double> delta;
	delta.setScale(dt, P2*h0i*transpose(d0));
    Sparse<double> delta2;
	delta2.setScale(dt, h0i*transpose(d0));
	Sparse<double> gamma;
	gamma.setScale(dt, h1*d0);
    Diagonal<double> P1(s_m,0.0);
    P1.setScale(0.055555*dt,h0i*P1);
	double solt[T*time_steps];
    double solt2[T*time_steps];

    for(double i=0; i<T*time_steps; i++) {
		e+=delta*h-P1*e;
        e2+=delta2*h2-P1*e2;
        for(uint j=0; j<source.size(); j++){
            e.m_val[source[j]]+=exp(-(i*dt-3*lev)*(i*dt-3*lev)/(lev*lev))*sin(PIx2*fm*(i*dt-3*lev));
            e2.m_val[source[j]]+=exp(-(i*dt-3*lev)*(i*dt-3*lev)/(lev*lev))*sin(PIx2*fm*(i*dt-3*lev));
        }
        h-=gamma*e;
        h2-=gamma*e2;
		solt[int(i+0.1)]=e2.m_val[tp];
        solt2[int(i+0.1)]=e.m_val[tpr];
    }

	ofstream myfile;
	myfile.open ("solutions.txt");
	for(uint i=0;i<T*time_steps;i++){
		myfile << solt2[i] << " " << solt[i] << endl;
	}

    double dualsum=0;
    double holearea=0;
    for(uint i=0;i<mesh.getNodeSize();i++){
		dualsum+=1.0/h0i.getValue(i);
        if(holes[i]!=1) holearea+=1.0/h0i.getValue(i);
	}
    cout << "sum of dual areas: " << dualsum << " Should be: " << 10 << endl;
	cout << "area of holes: " << holearea/6.0 << " Should be: " << PI*r*r << endl;

	myfile.close();
    gfd::finalizeMPI();
}
