#include "../../GFD/Types/Random.hpp"
#include "../../GFD/Mesh/BuilderMesh.hpp"
#include <iostream>
using namespace std;
using namespace gfd;

int main() {
    Random rnd(14);
    BuilderMesh mesh(3);
    cout << "hello to worlds" << endl;
    cout << rnd.getUint() << endl;
    cout << rnd.getUint() << endl;
    cout << rnd.getUint() << endl;
    cout << rnd.getUint() << endl;
    cout << rnd.getUint() << endl;
    cout << rnd.getUint() << endl;
    cout << rnd.getUint() << endl;
}