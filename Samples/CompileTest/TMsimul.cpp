/*Yritetään toistaa pyDEC:illä tuotettuja simulaatiotuloksia
TODO: generoidun verkon sijaan käytä samaa verkkoa kuin pyDECissä
TODO: sinimuotoinen aalto ja tuloksien talteenotto
*/




#include "../../GFD/Mesh/BuilderMesh.hpp"	//for building and modifying a mesh
#include "../../GFD/Mesh/PartMesh.hpp"	//for working with parts of a mesh?
#include "../../GFD/Types/Types.hpp"	//constants, shorter names for types, type couple
#include "../../GFD/Types/MpiEasy.hpp"	//for parallel computing
#include "../../GFD/Discrete/Dec.hpp"	//for dec related calculations
#include <iostream>

using namespace std;
using namespace gfd;

const uint space_steps = 64; // number of space elements per unit length
const uint time_steps = 128; // number of time elements per unit time

int main() {
	initMPI();

    PartMesh mesh(0, 1, 2);
    BuilderMesh bmesh(2);
    const double dx = 1.0 / double(space_steps);
    bmesh.createTriangleGrid(Vector2(0,0), Vector2(1,1), dx);
    mesh.swap(bmesh);
    Dec dec(mesh, 0, bmesh.getDimension());

    Derivative d0;
    d0=dec.integrateDerivative(fg_prim0, d0);
    Diagonal<double> h1;
    h1=dec.integrateHodge(HodgeUnit1, 0, fg_prim1, h1);
    Diagonal<double> h0i;
    h0i=dec.integrateHodge(HodgeUnit2, 0, fg_dual2, h0i);
    Column<double> h(0.0);
	dec.integrateZeroForm(fg_dual1, h);
	Column<double> e(0.0);
	dec.integrateZeroForm(fg_prim0, e);
    Sparse<double> delta1;
	const double dt = 1.0 / double(time_steps);
	delta1.setScale(dt, h0i*transpose(d0));
    Sparse<double> delta2;
    delta1.setScale(dt, h1*d0);

    for(uint i=0; i<time_steps; i++) {
        e+=delta1*h;
        h-=delta2*e;
    }

    finalizeMPI();
}