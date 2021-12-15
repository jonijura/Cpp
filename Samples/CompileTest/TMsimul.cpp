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
#include <iomanip>
#include <iostream>
#include <fstream>
#include <chrono>

using namespace std;
using namespace gfd;

const uint space_steps = 21; // number of space elements per unit length
const uint time_steps = 4*space_steps; // number of time elements per unit time
const uint T = 200; //simulation time
const uint grid_type = 1; // 2d grid type: 0 = squares, 1 = triangles, 2 = snubsquare
const double fm = 0.5;
const double lev = 3.0/(PI*fm);
const double epsilon = 5.9;
const double r = 0.3;
const double korjaus = 0;
const bool saveresult = true;

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
void loadMesh(PartMesh &mesh, string file, Buffer<double> &boundaryNodes){
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
	//Buffer<uint> flags(numFace*4);//record edges that are linked to only one face -> boundary edges
	//flags.fill(0);
    for(uint i=0; i<numFace; i++) {
		const Buffer<uint> face = getUints(input.getRow());
		Buffer<uint> edge(face.size());
		for(uint i=0; i<face.size(); i++){
			edge[i] = mesh.addEdge(face[i%face.size()]-1, face[(i+1)%face.size()]-1);
			//flags[edge[i]]++;
		}
		mesh.addFace(edge);
	}
	/*
	flags.resize(mesh.getEdgeSize());
	for(uint i=0; i<flags.size();i++){
		flags[i]=flags[i]%2;//flags[i]=1 iff edge number i is on the boundary
	}
	
	//for(uint i=0;i<flags.size();i++) cout << flags[i] << " ";
	Buffer<double> flagsV(mesh.getNodeSize());
	flagsV.fill(0.0);
	for(uint i=0; i<flags.size();i++){
		if(flags[i]==1){//flags[i]=1 iff node number i is on the boundary
			flagsV[mesh.getEdgeNodes(i)[0]]=1.0;
			flagsV[mesh.getEdgeNodes(i)[1]]=1.0;
			//mesh.setNodeFlag(mesh.getEdgeNodes(i)[0],1);
			//mesh.setNodeFlag(mesh.getEdgeNodes(i)[1],1);
		}
	}
	//for(uint i=0;i<flagsV.size();i++) cout << flagsV[i] << " ";
	boundaryNodes=flagsV;
	*/
    if(input.getRow().substr(0, 9).compare("$EndFaces") == 0) std::cout << "Faces and edges read succesfully." << endl;
}

int main(int argc, const char* argv[]) {
	initMPI();

	auto starttime = chrono::system_clock::now();

    PartMesh mesh(0, 1, 2);
	Buffer<double> boundaryNodes;
    if(argc==1){
        cout << "give input mesh" << endl;
		return 0;
    }
    else{
        loadMesh(mesh,argv[1],boundaryNodes);
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
        if(abs(xn-1)+abs(yn-0.5)<0.8/space_steps){
            tp=i;
            cout << "Tarkastelupiste 1: " << xn << ", " << yn << endl;
        }
        if(abs(xn-9)+abs(yn-0.5)<0.8/space_steps){
            tpr=i;
            cout << "Tarkastelupiste 2: " << xn << ", " << yn << endl;
        } 
		xm=fmod(xn,1.0)-0.5;
		ym=fmod(yn,1.0)-0.5;
		if((xn>2)&(xn<8)&(xm*xm+ym*ym)<(r*r+korjaus)) holes[i]=1.0/epsilon;
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
    P1.setScale(dx*dt,h0i*P1);
	double solt[4*T*time_steps];

    for(double i=0; i<T*time_steps; i++) {
		e+=delta*h-P1*e;
        e2+=delta2*h2-P1*e2;
        for(uint j=0; j<source.size(); j++){
            e.m_val[source[j]]+=exp(-(i*dt-3*lev)*(i*dt-3*lev)/(lev*lev))*sin(PIx2*fm*(i*dt-3*lev));
            e2.m_val[source[j]]+=exp(-(i*dt-3*lev)*(i*dt-3*lev)/(lev*lev))*sin(PIx2*fm*(i*dt-3*lev));
        }
        h-=gamma*e;
        h2-=gamma*e2;
		solt[int(4*i+0.1)]=e.m_val[tp];
        solt[int(4*i+1.1)]=e.m_val[tpr];
		solt[int(4*i+2.1)]=e2.m_val[tp];
        solt[int(4*i+3.1)]=e2.m_val[tpr];
		if(int(i)%time_steps==0)cout << '*';
    }
	cout << endl;
	ofstream myfile;
	myfile << setprecision(14);
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

	if(saveresult){
		srand((unsigned)time(0)*1942302);
		int name = rand()*1000000+rand();
		string fileloc = "C:\\MyTemp\\Tutkimusavustaja2021\\Simulointitulokset\\" + to_string(name) + ".txt";
		myfile.open (fileloc);
		myfile << setprecision(14);
		for(uint i=0;i<T*time_steps;i++){
			myfile << solt[4*i] << "\t" << solt[4*i+1] << "\t" << solt[4*i+2] << "\t" << solt[4*i+3] << endl;
		}
		myfile.close();
		
		uint e_n = mesh.getEdgeSize();
		double min = dx;
		for(uint i=0; i<e_n; i++)
		{
			const Vector4 v = mesh.getEdgeVector(i);
			double len = v.len();
			if(len < min) min = len;
		}

		fileloc = "C:\\MyTemp\\Tutkimusavustaja2021\\Simulointitulokset\\simulation_details.txt";
		myfile.open(fileloc, std::ios_base::app);
		myfile << setprecision(14);
		myfile << r << "\t" << 
			sqrt(holearea/(PI*6)) << "\t" << 
			epsilon << "\t" <<
			min << "\t" << 
			dx << "\t" <<
			dt << "\t" <<
			0 << "\t" <<
			T*time_steps << "\t" <<
			20011111100 << "\t" <<
			fm << "\t" <<
			1/(PI*lev) << "\t" <<
			mesh.getNodePosition2(tp).y << "\t" <<
			mesh.getNodePosition2(tp).x << "\t" <<
			mesh.getNodePosition2(tpr).y << "\t" <<
			mesh.getNodePosition2(tpr).x << "\t" <<
			name << endl;
		myfile.close();
	}

	gfd::finalizeMPI();
}