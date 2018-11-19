// Gmsh project created on Thu Nov 15 13:47:20 2018

lc1 = 0.005;
lc2 = 0.005;

r = 0.08*0.2;

Point(0) = {0, 0, 0, lc1};
Point(1) = {0.1, 0, 0, lc1};
Point(2) = {0, 0.1, 0, lc1};
Point(3) = {0.1, 0.1, 0, lc1};

Point(4)  = {0.05,   0.05,   0, lc1};
Point(5)  = {0.05+r, 0.05,   0, lc2};
Point(6)  = {0.05,   0.05+r, 0, lc2};
Point(7)  = {0.05-r, 0.05,   0, lc2};
Point(8)  = {0.05,   0.05-r, 0, lc2};
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 1};
//+
Line(3) = {1, 0};
//+
Line(4) = {0, 2};
//+
Circle(5) = {5, 4, 6};
//+
Circle(6) = {6, 4, 7};
//+
Circle(7) = {7, 4, 8};
//+
Circle(8) = {8, 4, 5};
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Line Loop(2) = {6, 7, 8, 5};
//+
Plane Surface(1) = {1, 2};
//+
Plane Surface(2) = {2};
//+
Extrude {0, 0, 1} {
  Surface{1}; Layers{50};
}
//+
Extrude {0, 0, 1} {
  Surface{2}; Layers{50};
}
//+
Physical Volume(1) = {1};
Physical Volume(2) = {2};
