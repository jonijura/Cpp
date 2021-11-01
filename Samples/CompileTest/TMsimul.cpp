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

const uint space_steps = 5; // number of space elements per unit length
const uint time_steps = 1000; // number of time steps
const double dt = 0.003; //time step size

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
	cout << "---------------"<<endl;
}
void loadMesh(PartMesh &mesh, string file){
    Text input;
    if(input.load(file)) cout << "mesh file loaded succesfully." << endl;
    while(input.hasRow()) {
		if(input.getRow().substr(0, 6).compare("$Nodes") == 0){
            cout << "reading nodes." << endl;
            break;
        } 
	}
    const uint numVert=getUint(input.getRow(),0);
    cout << numVert << endl;
    mesh.resizeNodeBuffer(numVert);
	for(uint i=0; i<numVert; i++) {
		const Buffer<double> pos = getDoubles(input.getRow());
		mesh.addNode(Vector4(pos[0], (pos.size() > 1 ? pos[1] : 0.0), 
													(pos.size() > 2 ? pos[2] : 0.0), 0.0));
	}
    if(input.getRow().substr(0, 9).compare("$EndNodes") == 0) cout << "Nodes read succesfully." << endl;
	if(input.getRow().substr(0, 6).compare("$Faces") == 0){
        cout << "reading faces and edges." << endl;
    }
    const uint numFace=getUint(input.getRow(),0);
    for(uint i=0; i<numFace; i++) {
		const Buffer<uint> face = getUints(input.getRow());
		Buffer<uint> edge(face.size());
		for(uint i=0; i<face.size(); i++){
			edge[i] = mesh.addEdge(face[i%face.size()]-1, face[(i+1)%face.size()]-1);
		}
		mesh.addFace(edge);
	}
    cout << numFace << endl;
    //read faces and edges
    if(input.getRow().substr(0, 9).compare("$EndFaces") == 0) cout << "Faces and edges read succesfully." << endl;
}

int main(int argc, const char* argv[]) {
	initMPI();

	auto starttime = chrono::system_clock::now();

    PartMesh mesh(0, 1, 2);
    if(argc==1){
        BuilderMesh bmesh(2);
        const double dx = 1.0 / double(space_steps);
        bmesh.createTriangleGrid(Vector2(0,0), Vector2(1,1), dx);
        mesh.swap(bmesh);
    }
    else{
        loadMesh(mesh,argv[1]);
    }
    Dec dec(mesh, 0, mesh.getDimension());

    Derivative d0; 
    d0=dec.integrateDerivative(fg_prim0, d0);
    Diagonal<double> h1;
    h1=dec.integrateHodge(HodgeUnit2, 0, fg_prim1, h1);//mistä tulee HodgeUnit2??
    Diagonal<double> h0i;
    h0i=dec.integrateHodge(HodgeUnit1, 0, fg_dual0, h0i);
    Column<double> h(0.0);
	dec.integrateZeroForm(fg_dual1, h);
	Column<double> e(0.0);
	dec.integrateZeroForm(fg_prim0, e);
    Sparse<double> delta;
	delta.setScale(dt, h0i*transpose(d0));
	Sparse<double> gamma;
	gamma.setScale(dt, h1*d0);
	double sol44[time_steps];
	double sol70[time_steps];

	cout << "Mesh and derivates generated in " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;
	starttime = chrono::system_clock::now();

    for(double i=0; i<time_steps; i++) {
		e+=delta*h;
        e.m_val[1700]=sin(PIx2*i/300.0);//60 158
        h-=gamma*e;
		sol44[int(i+0.1)]=e.m_val[7505];//36 44
		sol70[int(i+0.1)]=e.m_val[4505];//24 70
    }
	
	cout << "Timeintegration performed in " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;
	starttime = chrono::system_clock::now();

	ofstream myfile;
	myfile.open ("solutions.txt");
	for(uint i=0;i<time_steps;i++){
		myfile << sol44[i] << " " << sol70[i] << endl;
	}
	myfile.close();
    finalizeMPI();
}