// Gmsh project created on Wed Nov  7 09:52:10 2018

SetFactory("Built-in");

L   = 1;
W   = 0.1;
R   = 0.01;

lc1 = 0.01*L;
lc2 = 0.01*W;

Point(0) = {0, 0, 0, lc1};
Point(1) = {W, 0, 0, lc1};
Point(2) = {0, W, 0, lc1};
Point(3) = {W, W, 0, lc1};

Point(4) = {W/2, W/2,   0, lc1};
Point(5) = {W/2+R, W/2, 0, lc2};
Point(6) = {W/2, W/2+R, 0, lc2};
Point(7) = {W/2-R, W/2, 0, lc2};
Point(8) = {W/2, W/2-R, 0, lc2};

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

/*
Point(9)  = {0, 0, L, lc1};
Point(10) = {W, 0, L, lc1};
Point(11) = {0, W, L, lc1};
Point(12) = {W, W, L, lc1};

Point(13) = {W/2, W/2,   L, lc2};
Point(14) = {W/2+R, W/2, L, lc2};
Point(15) = {W/2, W/2+R, L, lc2};
Point(16) = {W/2-R, W/2, L, lc2};
Point(17) = {W/2, W/2-R, L, lc2};

//+
Line(9) = {12, 10};
//+
Line(10) = {10, 9};
//+
Line(11) = {9, 11};
//+
Line(12) = {11, 12};
//+
Circle(13) = {17, 13, 14};
//+
Circle(14) = {14, 13, 15};
//+
Circle(15) = {15, 13, 16};
//+
Circle(16) = {16, 13, 17};
//+
Line(17) = {12, 3};
//+
Line(18) = {10, 1};
//+
Line(19) = {9, 0};
//+
Line(20) = {11, 2};
//+
Line(21) = {17, 8};
//+
Line(22) = {14, 5};
//+
Line(23) = {15, 6};
//+
Line(24) = {16, 7};
//+
Line Loop(1) = {12, 17, -1, -20};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {9, 18, -2, -17};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {11, 20, -4, -19};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {10, 19, -3, -18};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {12, 9, 10, 11};
//+
Line Loop(6) = {14, 15, 16, 13};
//+
Plane Surface(5) = {5, 6};
//+
Line Loop(7) = {13, 14, 15, 16};
//+
Plane Surface(6) = {7};
//+
Line Loop(8) = {3, 4, 1, 2};
//+
Line Loop(9) = {7, 8, 5, 6};
//+
Plane Surface(7) = {8, 9};
//+
Line Loop(10) = {8, 5, 6, 7};
//+
Plane Surface(8) = {10};
*/

//+
Line Loop(1) = {4, 1, 2, 3};
//+
Line Loop(2) = {6, 7, 8, 5};
//+
Plane Surface(1) = {1, 2};
//+
Plane Surface(2) = {2};
//+
Extrude {0, 0, L} {
  Surface{1}; Surface{2}; Layers{1/lc1};
}
