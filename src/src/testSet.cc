#include <iostream>
#include "sets.hpp"

using namespace std;

void printSet(StdSet<int> set) {
  cout << "Printing set:" << endl;
  for (auto element : set) {
    cout << element << endl;
  }
}

void testCreationFromVector() {
  cout << "Testing creation from vector..." << endl;
  vector<int> v = {1, 2, 3, 4};
  StdSet<int> set_a = StdSet<int>(v.begin(), v.end());
  if (set_a.size() == 4) {
    cout << "Success!" << endl;
    return;
  }
  cout << "Fail!" << endl;
}

void testCreationFromSet() {
  cout << "Testing creation from set..." << endl;
  std::set<int> set_a;
  set_a.insert(1);
  set_a.insert(2);
  StdSet<int> set_b = StdSet<int>(set_a.begin(), set_a.end());
  set_a.erase(2);
  if (set_b.size() == 2) {
    cout << "Success!" << endl;
    return;
  }
  cout << "Fail!" << endl;
}

void testInsertion() {
  cout << "Testing set insertion..." << endl;
  StdSet<int> set_a;
  set_a.insert(1);
  set_a.insert(1);
  set_a.insert(2);
  set_a.insert(1);
  set_a.insert(10);
  if (set_a.size() == 3) {
    cout << "Success!" << endl;
    return;
  }
  cout << "Fail!" << endl;
}

void testVectorInsertion() {
  cout << "Testing set insertion with vectors..." << endl;
  StdSet<int> set_a;
  vector<int> v = {1, 2, 3, 4};
  set_a.insert(v);
  if (set_a.size() == 4) {
    cout << "Success!" << endl;
    return;
  }
  cout << "Fail!" << endl;
}

void testIntersection() {
  cout << "Testing set intersection..." << endl;
  StdSet<int> set_a;
  StdSet<int> set_b;
  set_a.insert(1);
  set_a.insert(2);
  set_a.insert(3);
  set_b.insert(1);
  set_b.insert(2);
  StdSet<int> intersection = set_a.intersect(set_b);
  if (intersection.size() == 2) {
    cout << "Success!" << endl;
    return;
  }
  cout << "Fail!" << endl;
}

void testEmpty() {
  cout << "Testing empty..." << endl;
  StdSet<int> set_a;
  if (!set_a.empty()) {
    cout << "Fail!" << endl;
  }
  cout << "Success!" << endl;
  set_a.insert(1);
  if (set_a.empty()) {
    cout << "Fail!" << endl;
    return;
  }
  cout << "Success!" << endl;
}

void testContains() {
  cout << "Testing contains..." << endl;
  StdSet<int> set_a;
  set_a.insert(1);
  if (set_a.contains(1)) {
    cout << "Success!" << endl;
    return;
  }
  cout << "Fail!" << endl;
}

void testErase() {
  cout << "Testing erase..." << endl;
  StdSet<int> set_a;
  int e = 1;
  set_a.insert(e);
  set_a.erase(e);
  if (set_a.contains(e)) {
    cout << "Fail!" << endl;
    return;
  }
  cout << "Success!" << endl;
}

void testIteratorErase() {
  cout << "Testing iterator erase..." << endl;
  StdSet<int> set_a;
  set_a.insert(1);
  auto it = set_a.begin();
  if (set_a.erase(it) == set_a.end()) {
    cout << "Success!" << endl;
    return;
  }
  cout << "Fail!" << endl;
}

int main(){
    testCreationFromVector();
    testCreationFromSet();
    testInsertion();
    testVectorInsertion();
    testIntersection();
    testEmpty();
    testContains();
    testErase();
    testIteratorErase();
    return 0;
}
