
lx=1;
ly=1.;
lz=1.;

lc=lx/5;

Point(1) = {0. , 0., 0, lc};
Point(2) = {lx , 0., 0, lc};
Point(3) = {lx , 0 , lz*0.5, lc};
Point(6) = {0., 0. , lz*0.5, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 6};
Line(6) = {6, 1};


Point(11) = {0. , ly, 0, lc};
Point(12) = {lx , ly, 0, lc};
Point(13) = {lx , ly, lz*0.5, lc};
Point(16) = {0. , ly , lz*0.5, lc};

Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 16};
Line(16) = {16, 11};

Line(17) = {13, 3};
Line(18) = {2, 12};
Line(21) = {6, 16};
Line(22) = {1, 11};

Line Loop(23) = {17, -2, 18, 12};
Plane Surface(24) = {23};
//Line Loop(25) = {19, 14, -20, -4};
//Plane Surface(26) = {25};
//Line Loop(27) = {5, 21, -15, -20};
//Plane Surface(28) = {27};
Line Loop(29) = {3, 21, -13, 17};
Plane Surface(30) = {29};
Line Loop(31) = {1, 2, 3, 6};
Plane Surface(32) = {31};
Line Loop(33) = {11, -18, -1, 22};
Plane Surface(34) = {33};
Line Loop(35) = {21, 16, -22, -6};
Plane Surface(36) = {35};
Line Loop(37) = {12, 13, 16, 11};
Plane Surface(38) = {37};
Coherence;