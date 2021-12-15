/*Yritetään toistaa pyDEC:illä tuotettuja simulaatiotuloksia

MPI ei toimi, tekee saman homman useammassa prosessissa. muuta kääntäjäasetuksia?
Kuinka voin debugatessa tarkastella bufferin lukuja ilman segventation faulttia?

koodipohjan ongelmia:
index out of bounds ei kaada ohjelmaa tai heitä varoitusta
väärädimensioisien matriisien kertominen ei kaada ohjelmaa tai heitä varoitusta
verkkojen tallennusmuoto ei ole järkevä, pitäisi olla ainakin numeropohjainen vaihtoehto

kommentit puuttuu lähes kaikkialta
testifunktiot puuttuu kaikkialta?
paljon turhaa toistoa
useita hyödyllisiä funktioita puuttuu, esimerkiksi funktioita ulkoisen verkon lukemiseen
*/



#include "../../GFD/Output/MeshDrawer.hpp"	//for drawing a mesh
#include "../../GFD/Mesh/BuilderMesh.hpp"	//for building and modifying a mesh
#include "../../GFD/Mesh/PartMesh.hpp"	//for working with parts of a mesh?
#include "../../GFD/Types/Types.hpp"	//constants, shorter names for types, type couple
#include "../../GFD/Types/MpiEasy.hpp"	//for parallel computing
#include "../../GFD/Discrete/Dec.hpp"	//for dec related calculations
#include <iostream>
#include <fstream>
#include <chrono>

using namespace std;
using namespace gfd;

const uint space_steps = 15; // number of space elements per unit length
const uint time_steps = 4*space_steps; // number of time elements per unit time
const uint T = 60; //simulation time
const double fm = 0.5;
const double lev = 20/(PI*fm);
const double epsilon = 8.9;
const double r = 0.2;
const double korjaus = 0;

double getDouble(const std::string &str, const double def) {
	char *rest;
	const double val = strtod(str.c_str(), &rest);
	if(rest == NULL) return def;
	return val;
}
Buffer<double> getDoubles(const std::string &str) {
	Buffer<double> res;
	const string space = " \t";
	size_t pos = 0;
	pos = str.find_first_not_of(space, pos);
	while(pos != string::npos) {
		const size_t until = str.find_first_of(space, pos);
		res.push_back(getDouble(str.substr(pos, until), 0.0));
		pos = str.find_first_not_of(space, until);
	}
	return res;
}
uint getUint(const std::string &str, const uint def) {
	const double val = getDouble(str, double(def));
	if(val < 0.0) return 0;
	return uint(val + 0.5);
}
Buffer<uint> getUints(const std::string &str) {
	Buffer<uint> res;
	const string space = " \t";
	size_t pos = 0;
	pos = str.find_first_not_of(space, pos);
	while(pos != string::npos) {
		const size_t until = str.find_first_of(space, pos);
		res.push_back(getUint(str.substr(pos, until), 0));
		pos = str.find_first_not_of(space, until);
	}
	return res;
}

void line(){
	std::cout << "---------------"<<endl;
}
void loadMesh(PartMesh &mesh, string file){
    Text input;
    if(input.load(file)) std::cout << "mesh file loaded succesfully." << endl;
    while(input.hasRow()) {
		if(input.getRow().substr(0, 6).compare("$Nodes") == 0){
            std::cout << "reading nodes." << endl;
            break;
        } 
	}
    const uint numVert=getUint(input.getRow(),0);
    std::cout << numVert << endl;
    mesh.resizeNodeBuffer(numVert);
	for(uint i=0; i<numVert; i++) {
		const Buffer<double> pos = getDoubles(input.getRow());
		mesh.addNode(Vector4(pos[0], pos[1], 0.0, 0.0));
	}
    if(input.getRow().substr(0, 9).compare("$EndNodes") == 0) std::cout << "Nodes read succesfully." << endl;
	if(input.getRow().substr(0, 6).compare("$Faces") == 0){
        std::cout << "reading faces and edges." << endl;
    }
    const uint numFace=getUint(input.getRow(),0);
	std::cout << numFace << endl;
    for(uint i=0; i<numFace; i++) {
		const Buffer<uint> face = getUints(input.getRow());
		Buffer<uint> edge(face.size());
		for(uint i=0; i<face.size(); i++){
			edge[i] = mesh.addEdge(face[i%face.size()]-1, face[(i+1)%face.size()]-1);
		}
		mesh.addFace(edge);
	}
    if(input.getRow().substr(0, 9).compare("$EndFaces") == 0) std::cout << "Faces and edges read succesfully." << endl;
}

int main(int argc, const char* argv[]) {
	initMPI();
    PartMesh mesh(0, 1, 2);
    if(argc==1){
        cout << "give input mesh" << endl;
		return 0;
    }
    else{
        loadMesh(mesh,argv[1]);
    }

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
    uint source;
	double xn,yn,xm,ym;
    double eps = 1e-5;
    uint tp=1,tpr;
	for(double i=0;i<mesh.getNodeSize();i++){
        s_m[i]=0;
        xn=mesh.getNodePosition2(i).x;
        yn=mesh.getNodePosition2(i).y;
        if((xn<eps) | (xn>15-eps) | (yn<eps) | (yn>7-eps)) s_m[i]=1.0;
        if(abs(xn)+abs(yn-3.5)<0.05){
            source=i;
            cout << "Lähdepiste: " << xn << ", " << yn << endl;
        }
        if(abs(xn-9)+abs(yn-3.5)<0.05){
            tpr=i;
            cout << "Tarkastelupiste: " << xn << ", " << yn << endl;
        } 
		xm=fmod(xn,1.0)-0.5;
		ym=fmod(yn,1.0)-0.5;
		if(((yn<3)|(yn>4))&(xm*xm+ym*ym)<(r*r+korjaus)) holes[i]=1.0/epsilon;
		else holes[i]=1;
	}
	Diagonal<double> P2(holes,0.0);

    Column<double> h(0.0);
	dec.integrateZeroForm(fg_dual1, h);
	Column<double> e(0.0);
	dec.integrateZeroForm(fg_prim0, e);
    Sparse<double> delta;
	delta.setScale(dt, P2*h0i*transpose(d0));
	Sparse<double> gamma;
	gamma.setScale(dt, h1*d0);
    Diagonal<double> P1(s_m,0.0);
    P1.setScale(dx*dt,h0i*P1);
	double solt[4*T*time_steps];

    for(double i=0; i<T*time_steps; i++) {
		e+=delta*h-P1*e;
        e.m_val[source]+=exp(-(i*dt-3*lev)*(i*dt-3*lev)/(lev*lev))*sin(PIx2*fm*(i*dt-3*lev));
        h-=gamma*e;
		solt[int(4*i+0.1)]=e.m_val[tp];
        solt[int(4*i+1.1)]=e.m_val[tpr];
		if(int(i)%time_steps==0)cout << '*';
    }
	cout << endl;
	ofstream myfile;
	
	myfile.open ("solutions.txt");
	uint n_n=mesh.getNodeSize();
	for(uint i=0;i<n_n;i++){
		myfile << e.m_val[i] << endl;
	}
	myfile.close();
	
    double dualsum=0;
    double holearea=0;
    for(uint i=0;i<mesh.getNodeSize();i++){
		dualsum+=1.0/h0i.getValue(i);
        if(holes[i]!=1) holearea+=1.0/h0i.getValue(i);
	}
    cout << "sum of dual areas: " << dualsum << " Should be: " << 10 << endl;
	cout << "area of holes: " << holearea/6.0 << " Should be: " << PI*r*r << " difference: " << abs(PI*r*r-holearea/6.0) << endl;

	gfd::finalizeMPI();
}