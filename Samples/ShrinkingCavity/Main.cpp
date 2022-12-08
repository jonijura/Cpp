/**
 * Program for computing (2+1)-dimensional wave in a rotating cavity.
 * Idea of the solution is following:
 * - create (2+1)-dimensional space-time mesh.
 * - consider discrete forms of grade 2 (values located at faces).
 * - initialize wave at t=0.
 * - solve field finding gradient zero dF=0 (using Minkowski-metric).
 * - in this occasion, dF=0 leads to explicit and causal system.
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include "../../GFD/Types/Types.hpp"
#include "../../GFD/Types/Random.hpp"
#include "../../GFD/Types/Buffer.hpp"
#include <iostream>
#include <chrono>

using namespace std;
using namespace gfd;

TwoVector3 getField(const Mesh& mesh, const Buffer<double>& val, const uint node) {
	const Buffer<uint> par = mesh.getNodeFaces(node);
	SymTwoMatrix3 A(0, 0, 0, 0, 0, 0);
	TwoVector3 b(0, 0, 0);
	for (uint i = 0; i < par.size(); i++) {
		const TwoVector3 v = mesh.getFaceVector(par[i]).toTwoVector3();
		A += v.outerProduct();
		b += v * val[par[i]];
	}
	return A.inverse() * b;
}

Buffer<uint> findEdges(const uint n0, const uint n1, const uint n2, BuilderMesh& mesh) {
	Buffer<uint> e(3);
	e[0] = mesh.findEdge(n0, n1);
	e[1] = mesh.findEdge(n1, n2);
	e[2] = mesh.findEdge(n2, n0);
	return e;
}
Buffer<uint> findFaces(const uint n0, const uint n1, const uint n2, const uint n3, BuilderMesh& mesh) {
	Buffer<uint> f(4);
	f[0] = mesh.findFace(findEdges(n0, n1, n2, mesh));
	f[1] = mesh.findFace(findEdges(n0, n1, n3, mesh));
	f[2] = mesh.findFace(findEdges(n0, n2, n3, mesh));
	f[3] = mesh.findFace(findEdges(n1, n2, n3, mesh));
	return f;
}
Buffer<uint> findBodies(const uint n0, const uint n1, const uint n2, const uint n3, const uint n4, BuilderMesh& mesh) {
	Buffer<uint> b(5);
	b[0] = mesh.findBody(findFaces(n0, n1, n2, n3, mesh));
	b[1] = mesh.findBody(findFaces(n0, n1, n2, n4, mesh));
	b[2] = mesh.findBody(findFaces(n0, n1, n3, n4, mesh));
	b[3] = mesh.findBody(findFaces(n0, n2, n3, n4, mesh));
	b[4] = mesh.findBody(findFaces(n1, n2, n3, n4, mesh));
	return b;
}
uint createSpaceTimeMesh(const double h, const double dtime, const uint steps, const double slope, const bool async, BuilderMesh& mesh) {
	uint i, j, k, l, m;
	const double N = 1.0 / h + 1e-5;
	const double n = 0.4 / h + 1e-5;
	const double SQRT3_4 = sqrt(0.75);
	BuilderMesh basemesh(2);
	basemesh.createTriangleGrid(Vector2(-1, -SQRT3_4), Vector2(1, SQRT3_4), h, false);
	const Vector2 e0 = Vector2(0.0, SQRT3_4) / (h * 0.75);
	const Vector2 e1 = Vector2(0.75, 0.5 * SQRT3_4) / (h * 0.75);
	const Vector2 e2 = Vector2(0.75, -0.5 * SQRT3_4) / (h * 0.75);
	for (i = basemesh.getNodeSize(); i-- > 0; ) {
		const Vector2 p = basemesh.getNodePosition2(i);
		if (e0.dot(p) > N || e0.dot(p) < -N) basemesh.removeNode(i);
		else if (e1.dot(p) > N || e1.dot(p) < -0.3 / h - 0.001) basemesh.removeNode(i);
		else if (e2.dot(p) > N || e2.dot(p) < -N) basemesh.removeNode(i);
		else if (e0.dot(p) < n && e1.dot(p) < n && e2.dot(p) < n) basemesh.removeNode(i);
	}
	const uint nodes = basemesh.getNodeSize();
	for (i = 0; i < nodes; i++) {
		const Vector2 x = basemesh.getNodePosition2(i);
		basemesh.setNodeFlag(i, uint((3.0 + x.x + 1.5 * x.y / sqrt(0.75)) / h + 0.5) % 3);
	}

	// create space-time mesh
	mesh.setMetric(SymMatrix4(1, 0, 1, 0, 0, -1, 0, 0, 0, 0));
	for (i = 0; i <= steps; i++) {
		for (j = 0; j < nodes; j++) {
			const Vector2 x = basemesh.getNodePosition2(j);
			const double z = dtime * (async ? i - 1.0 + double(basemesh.getNodeFlag(j)) / 3.0 : i);
			const Vector2 b(-0.0 * z, 0);
			const double cosi = cos(slope * z);
			const double sini = sin(slope * z);
			const Matrix2 a(cosi, sini, -sini, cosi);
			mesh.addNode(Vector4(a * x + b, z, 0.0));
		}
	}
	for (j = 0; j < basemesh.getEdgeSize(); j++) {
		const Buffer<uint>& n = basemesh.getEdgeNodes(j);
		mesh.addEdge(n[0], n[1]);
	}
	for (j = 0; j < basemesh.getFaceSize(); j++) mesh.addFace(basemesh.getFaceEdges(j));
	for (j = 0; j < basemesh.getBodySize(); j++) mesh.addBody(basemesh.getBodyFaces(j));
	Buffer<uint> nn(nodes);
	for (k = 0; k < nodes; k++) nn[k] = k;
	for (i = 1; i <= steps; i++) {
		for (j = 0; j < 3; j++) {
			for (k = 0; k < nodes; k++) {
				if (basemesh.getNodeFlag(k) != j) continue;
				nn[k] += nodes;
				mesh.addEdge(nn[k], nn[k] - nodes);
				const Buffer<uint>& e = basemesh.getNodeEdges(k);
				for (l = 0; l < e.size(); l++) {
					Buffer<uint> en = basemesh.getEdgeNodes(e[l]);
					for (m = 0; m < en.size(); m++) en[m] = nn[en[m]];
					mesh.addEdge(en[0], en[1]);
					mesh.addFace(findEdges(en[0], en[1], nn[k] - nodes, mesh));
				}
				const Buffer<uint> f = basemesh.getNodeFaces(k);
				for (l = 0; l < f.size(); l++) {
					Buffer<uint> fn = basemesh.getFaceNodes(f[l]);
					for (m = 0; m < fn.size(); m++) fn[m] = nn[fn[m]];
					mesh.setFaceFlag(mesh.addFace(findEdges(fn[0], fn[1], fn[2], mesh)), 1);
					mesh.addBody(findFaces(fn[0], fn[1], fn[2], nn[k] - nodes, mesh));
				}
				const Buffer<uint> b = basemesh.getNodeBodies(k);
				for (l = 0; l < b.size(); l++) {
					Buffer<uint> bn = basemesh.getBodyNodes(b[l]);
					for (m = 0; m < bn.size(); m++) bn[m] = nn[bn[m]];
					mesh.addBody(findFaces(bn[0], bn[1], bn[2], bn[3], mesh));
					mesh.addQuad(findBodies(bn[0], bn[1], bn[2], bn[3], nn[k] - nodes, mesh));
				}
			}
		}
	}
	return basemesh.getFaceSize() + (mesh.getFaceSize() - basemesh.getFaceSize()) / steps; // return number of faces in a block
}

void drawMeshBar(const double h, const double dtime, const uint steps, const double slope, const uint picwidth = 400) {
	BuilderMesh mesh(3);
	createSpaceTimeMesh(h, dtime, steps, slope, false, mesh);
	for (uint i = 0; i <= steps; i++) {
		cout << "Drawing mesh " << i << "..." << endl;
		const double height = dtime * i;
		Buffer<Vector3> col(mesh.getFaceSize());
		for (uint j = 0; j < col.size(); j++) {
			if (mesh.getFaceAverage(j).z < height) col[j] = Vector3(1, 0, 0);
			else col[j] = Vector3(1, 1, 1);
		}

		MeshDrawer drawer;
		const Vector3 vo(0, 0, 3.7);
		const Vector3 vp(0, -100, 100);
		const Vector3 vx = TwoVector3(Vector3(0, 0, 1), vp).dual().unit() / 6.7;
		const Vector3 vy = TwoVector3(vp, vx).dual().unit() / 6.7;
		drawer.initPosition(Vector4(vp + vo, 0), Vector4(vo, 0), Vector4(vx, 0), Vector4(vy, 0));
		Picture pic(picwidth, picwidth);
		drawer.initPicture(&pic);
		drawer.drawBoundaryFaces(mesh, col);
		Text path;
		path << "mesh" << i << ".bmp";
		pic.save(path.str(), false);
	}
}

void drawWave(const double h, const double dtime, const uint microsteps, const uint steps, const double slope, const uint picwidth = 400) {
	BuilderMesh mesh(3);
	const uint fblock = createSpaceTimeMesh(h, dtime, microsteps + 1, slope, true, mesh);
	const uint nblock = microsteps * (mesh.getNodeSize() / (microsteps + 2));

	Buffer<bool> edge_flags(mesh.getEdgeSize(), false); // ************************************************* flags for edges that are used to update values
	Buffer<bool> body_flags(mesh.getBodySize(), false); // ************************************************* flags for bodies that are used to update values

	// initialize values and operator
	uint i, j, k;
	const Vector2 p0(-0.7 * sin(PI / 6.0), 0.7 * cos(PI / 6.0));
	Buffer<double> val(mesh.getFaceSize(), 0.0);
	Buffer< Buffer< pair<uint, double> > > buf(val.size());
	for (i = 0; i < val.size(); i++) {
		if (i < fblock) { // copy value from the end
			const TwoVector3 tv = mesh.getFaceVector3(i);
			const Vector2 p = mesh.getFaceAverage3(i).toVector2() - p0;
			const double plen = 50.0 * p.len();
			const Buffer<uint> n = mesh.getFaceNodes(i);
			const uint ii = mesh.findFace(findEdges(n[0] + nblock, n[1] + nblock, n[2] + nblock, mesh));
			if (plen < PI) val[ii] = 3.0 * (1.0 + cos(plen)) * tv.xy;
			buf[i].push_back(pair<uint, double>(ii, 1.0));
		}
		else if (mesh.getFaceFlag(i) == 1) { // update by a body
			const Buffer<uint>& ele = mesh.getFaceBodies(i);
			uint body = ele[0];
			for (j = 1; j < ele.size(); j++) {
				if (ele[j] < body) body = ele[j];
			}
			//const uint body = mesh.getFaceBodies(i)[0];
			const double sign = -mesh.getBodyIncidence(body, i);
			const Buffer<uint>& f = mesh.getBodyFaces(body);
			buf[i].resize(f.size() - 1);
			for (j = 0, k = 0; j < f.size(); j++) {
				if (f[j] == i) continue;
				buf[i][k++] = pair<uint, double>(f[j], sign * mesh.getBodyIncidence(body, f[j]));
			}
			body_flags[body] = true; // ************************************************* mark the body to indicate that the condition dF=0 is imposed on it
		}
		else { // update by an edge
			const Buffer<uint>& ele = mesh.getFaceEdges(i);
			uint edge = ele[0];
			for (j = 1; j < ele.size(); j++) {
				if (ele[j] < edge) edge = ele[j];
			}
			//const uint edge = mesh.getFaceEdges(i)[0];
			const double sign = -mesh.getFaceIncidence(i, edge) / mesh.getFaceHodge(i);
			const Buffer<uint>& f = mesh.getEdgeFaces(edge);
			buf[i].resize(f.size() - 1);
			for (j = 0, k = 0; j < f.size(); j++) {
				if (f[j] == i) continue;
				buf[i][k++] = pair<uint, double>(f[j], sign * mesh.getFaceIncidence(f[j], edge) * mesh.getFaceHodge(f[j]));
			}
			edge_flags[edge] = true; // ************************************************* mark the edge to indicate the condition *d*F=0 is imposed on it
		}
	}

	// iterate and draw fields
	const uint node0 = mesh.findNode(Vector4(p0, 0, 0), 0.5 * h * h, 0, false);
	const uint node1 = mesh.findNode(Vector4(0.7 * cos(PI / 6.0), -0.7 * sin(PI / 6.0), 0, 0), 0.5 * h * h, 0, false);
	const double z = dtime * microsteps;
	const double cosi = cos(slope * z);
	const double sini = sin(slope * z);
	const Matrix4 A(cosi, sini, 0, 0, -sini, cosi, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0);
	Picture picture(picwidth, picwidth);
	picture.fillColor(Vector4(-1, -1, -1, 0));
	for (uint iter = 0; iter <= steps; iter++) {
		cout << "Iteration " << iter << "..." << endl;

		// update values
		for (i = 0; i < buf.size(); i++) {
			val[i] = 0.0;
			for (j = 0; j < buf[i].size(); j++) val[i] += buf[i][j].second * val[buf[i][j].first];
		}

		// *************************************************
		if (iter == 1) { // check which conditions are satisfied by the solution after the first iteration
			
			// test which bodies satisfy the condition dF=0
			uint a = 0, b = 0, c = 0, d = 0;
			for (i = 0; i < mesh.getBodySize(); i++) {
				const Buffer<uint>& faces = mesh.getBodyFaces(i);
				double sum = 0;
				for (j = 0; j < faces.size(); j++)
					sum += mesh.getBodyIncidence(i, faces[j]) * val[faces[j]];
				if (fabs(sum) > 1e-9) {
					a++;
					if (body_flags[i] == true)
						b++;
				}
				else {
					c++;
					if (body_flags[i] == true)
						d++;
				}
			}
			cout << "Condition dF=0 satisfied on " << c << " bodies, " << d << " of which have body flag set\n";
			cout << "Condition dF=0 violated on " << a << " bodies, " << b << " of which have body flag set\n";

			// test which edges satisfy the condition *d*F=0
			Buffer<bool> satisfied(mesh.getEdgeSize(), true);
			a = 0, b = 0, c = 0, d = 0;
			for (i = 0; i < mesh.getEdgeSize(); i++) {
				bool boundary = false;
				const Buffer<uint>& fac = mesh.getEdgeFaces(i);
				for (j = 0; j < fac.size(); j++) {
					if (mesh.getFaceBodies(fac[j]).size() == 1)
						boundary = true;
				}
				if (boundary) continue; //the condition is violated throughout the boundary (comment this line to show results) 
				const Buffer<uint>& faces = mesh.getEdgeFaces(i);
				double sum = 0;
				for (j = 0; j < faces.size(); j++)
					sum += mesh.getFaceIncidence(faces[j], i) * mesh.getFaceHodge(faces[j]) * val[faces[j]];
				sum /= mesh.getEdgeHodge(i);
				if (abs(sum) > 1e-3) {
					a++;
					if (edge_flags[i] == true) {
						b++;
						cout << "Condition *d*F=0 violated on edge " << i << " although it was imposed. Value (should be zero): " << sum << endl;
					}
					else {
						cout << "Condition *d*F=0 violated on edge " << i << " (it was not imposed). Value (should be zero): " << sum << endl;
					}
					satisfied[i] = false;
				}
				else {
					c++;
					if (edge_flags[i] == true)
						d++;
				}
			}
			cout << "Condition *d*F=0 satisfied on " << c << " edges, " << d << " of which have edge flag set\n";
			cout << "Condition *d*F=0 violated on " << a << " edges, " << b << " of which have edge flag set\n";
			
			//print the edges on which the condition was violated and nearby mesh elements for Mathematica to visualise
			std::cout << "edges = {";
			for (uint i = 0; i < mesh.getEdgeSize(); ++i) {
				if (satisfied[i]) continue;
				Buffer<uint> nodes = mesh.getEdgeNodes(i);
				Vector3 p0 = mesh.getNodePosition3(nodes[0]);
				Vector3 p1 = mesh.getNodePosition3(nodes[1]);
				std::cout << "{{" << p0.x << ", " << p0.y << ", " << p0.z << "}, {";
				std::cout << p1.x << ", " << p1.y << ", " << p1.z << "}},\n";
			}
			std::cout << "};\n";
			std::cout << "elements = {";
			for (uint b = 0; b < mesh.getBodySize(); ++b) {
				Buffer<uint> nodes = mesh.getBodyNodes(b);
				if (nodes.size() == 4 && (mesh.getBodyPosition3(b).toVector2() - p0).lensq() < 0.01) {
					Vector3 p0 = mesh.getNodePosition3(nodes[0]);
					Vector3 p1 = mesh.getNodePosition3(nodes[1]);
					Vector3 p2 = mesh.getNodePosition3(nodes[2]);
					Vector3 p3 = mesh.getNodePosition3(nodes[3]);
					std::cout << "{{" << p0.x << ", " << p0.y << ", " << p0.z << "}, {";
					std::cout << p1.x << ", " << p1.y << ", " << p1.z << "}, {";
					std::cout << p2.x << ", " << p2.y << ", " << p2.z << "}, {";
					std::cout << p3.x << ", " << p3.y << ", " << p3.z << "}},\n";
				}
			}
			std::cout << "};\n";
		}

		// use the following command (with the lists of edges and elements printed above) to visualise with Mathematica
		// Show[Graphics3D[{{Cyan, EdgeForm[Thickness[0.002]], Opacity[0.3], Tetrahedron[elements]}, {Red, Thickness[0.02], Line[edges]}}], {Axes -> True, Boxed -> False}]
		// *************************************************

		// draw picture
		uint node = node0;
		const TwoVector3 tv = getField(mesh, val, node);
		Vector4 color(tv.xy, tv.xz, tv.yz, 1);
		for (j = 0; j < picture.getHeight(); j++) {
			for (i = 0; i < picture.getWidth(); i++) {
				if (iter > 0 && picture.getColor(i, j).t == 0.0 && i >= 5 && picture.getColor(i - 5, j).t == 0.0 &&
					i + 5 < picture.getWidth() && picture.getColor(i + 5, j).t == 0.0 &&
					j >= 5 && picture.getColor(i, j - 5).t == 0.0 &&
					j + 5 < picture.getHeight() && picture.getColor(i, j + 5).t == 0.0) continue;

				const Vector4 p(2.1 * i / double(picture.getWidth()) - 1.05, 2.1 * j / double(picture.getHeight()) - 1.05, 0, 0);
				uint newnode = mesh.findNode(p, 0.5 * h * h, node, false);
				if (newnode == NONE) {
					const double sq0 = (mesh.getNodePosition(node0) - p).lensq();
					const double sq1 = (mesh.getNodePosition(node1) - p).lensq();
					if (sq0 < sq1) {
						if (sq0 < (mesh.getNodePosition(node) - p).lensq()) newnode = mesh.findNode(p, 0.5 * h * h, node0, false);
					}
					else if (sq1 < (mesh.getNodePosition(node) - p).lensq()) newnode = mesh.findNode(p, 0.5 * h * h, node1, false);
					if (newnode == NONE) {
						picture.setColor(i, j, Vector4(-1, -1, -1, 0));
						continue;
					}
				}
				if (newnode != node) {
					node = newnode;
					const TwoVector3 tv = getField(mesh, val, node);
					color = Vector4(tv.xy, tv.xz, tv.yz, 1);
				}
				picture.setColor(i, j, color);
			}
		}
		Text path;
		path << "field" << iter << ".bmp";
		picture.save(path.str(), true);

		// update mesh node positions for the next iteration
		for (i = 0; i < mesh.getNodeSize(); i++) {
			mesh.setNodePosition(i, A * mesh.getNodePosition(i));
		}
	}

	// draw and save statistics
	/*MeshDrawer drawer;
	const Vector3 vo(0, 0, 5);
	const Vector3 vp(0, 20, 10);
	const Vector3 vx = TwoVector3(Vector3(0, 0, 1), vp).dual().unit() / 3.0;
	const Vector3 vy = TwoVector3(vp, vx).dual().unit() / 3.0;
	drawer.initPosition(Vector4(vp + vo, 0), Vector4(vo, 0), Vector4(vx, 0), Vector4(vy, 0));
	Picture pic(picwidth, picwidth);
	drawer.initPicture(&pic);
	drawer.initSvg(picwidth, picwidth);
	drawer.drawBoundaryFaces(mesh, Vector3(1, 1, 1));
	pic.save("mesh.bmp", false);
	drawer.saveSvg("mesh.svg");

	// save mesh statistics
	Text text;
	mesh.writeStatistics(text);
	text.save("stat.txt");
	cout << "Saved statistics" << endl;*/
}

int main() {
	auto starttime = chrono::system_clock::now();

	drawMeshBar(1.0 / 4.0, 1.0 / 10.0, 80, PI / 12.0, 400);
	drawWave(1.0 / 40.0, 1.0 / 100.0, 10, 80, PI / 12.0, 400);

	cout << "Elapsed time: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;
	return 0;
}
