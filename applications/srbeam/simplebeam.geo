// Gmsh project created on Thu Nov 15 12:27:57 2018

lc1 = 0.005;
lc2 = 0.005;

Point(0) = {0,   0,   0, lc1};
Point(1) = {0.1, 0,   0, lc1};
Point(2) = {0,   0.1, 0, lc1};
Point(3) = {0.1, 0.1, 0, lc1};

Point(4) = {0,   0,   1, lc2};
Point(5) = {0.1, 0,   1, lc2};
Point(6) = {0,   0.1, 1, lc2};
Point(7) = {0.1, 0.1, 1, lc2};



//+
Line(1) = {4, 5};
//+
Line(2) = {5, 7};
//+
Line(3) = {7, 6};
//+
Line(4) = {6, 4};
//+
Line(5) = {5, 1};
//+
Line(6) = {1, 3};
//+
Line(7) = {3, 2};
//+
Line(8) = {2, 0};
//+
Line(9) = {0, 1};
//+
Line(10) = {3, 7};
//+
Line(11) = {6, 2};
//+
Line(12) = {0, 4};
//+
Line Loop(1) = {9, 6, 7, 8};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {11, 8, 12, -4};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {3, 11, -7, 10};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {6, 10, -2, 5};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {5, -9, 12, 1};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {1, 2, 3, 4};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {3, 6, 5, 4, 1, 2};
//+
Volume(1) = {1};
