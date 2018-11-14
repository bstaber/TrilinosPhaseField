// Gmsh project created on Sun Sep 23 18:09:47 2018
SetFactory("OpenCASCADE");

lc1 = 0.5;
lc2 = 0.025;

H = 0.1;

a = 2;
b = 5;
c = 10;
d = 3;

Point(1) = {0, 0, 0, lc1};
Point(2) = {b, 0, 0, lc1};
Point(3) = {0, a, 0, lc1};
Point(4) = {b, a, 0, lc1};

Point(5) = {0, 0, c, lc1};
Point(6) = {b, 0, c, lc1};
Point(7) = {0, a, c, lc1};
Point(8) = {b, a, c, lc1};

Point(9)  = {d, 0, b-H, lc2};
Point(10) = {b, 0, b-H, lc1};
Point(11) = {d, a, b-H, lc2};
Point(12) = {b, a, b-H, lc1};

Point(13) = {d, 0, b+H, lc2};
Point(14) = {b, 0, b+H, lc1};
Point(15) = {d, a, b+H, lc2};
Point(16) = {b, a, b+H, lc1};

Point(17) = {0, 0, b-H, lc2};
Point(18) = {0, a, b-H, lc2};
Point(19) = {0, 0, b+H, lc2};
Point(29) = {0, a, b+H, lc2};

//+
Line(1) = {2, 4};
//+
Line(2) = {2, 1};
//+
Line(3) = {1, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 12};
//+
Line(6) = {12, 10};
//+
Line(7) = {10, 2};
//+
Line(8) = {10, 9};
//+
Line(9) = {9, 13};
//+
Line(10) = {13, 14};
//+
Line(11) = {14, 16};
//+
Line(12) = {16, 15};
//+
Line(13) = {15, 11};
//+
Line(14) = {11, 12};
//+
Line(15) = {9, 11};
//+
Line(16) = {13, 15};
//+
Line(17) = {16, 8};
//+
Line(18) = {14, 6};
//+
Line(19) = {6, 8};
//+
Line(20) = {8, 7};
//+
Line(21) = {7, 5};
//+
Line(22) = {5, 6};
//+
Line(23) = {7, 29};
//+
Line(24) = {29, 19};
//+
Line(25) = {19, 5};
//+
Line(26) = {19, 17};
//+
Line(27) = {17, 18};
//+
Line(28) = {18, 29};
//+
Line(29) = {17, 1};
//+
Line(30) = {18, 3};
//+
Line(31) = {17, 9};
//+
Line(32) = {19, 13};
//+
Line(33) = {18, 11};
//+
Line(34) = {29, 15};
