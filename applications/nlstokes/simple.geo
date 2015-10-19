lc = .01;

Point(1) = {0., 0., 0., lc};
Point(2) = {.5, -.1, 0., lc}; 
Point(3) = {1., -.15, 0., lc};
Point(4) = {1.2, -.15, 0., lc};
Point(5) = {.3, .0, 0., lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,1};

Line Loop(6) = {1,2,3,4,5};
Plane Surface(7) = {6};

Physical Line("Bottom") = {1,2,3};
Physical Line("Top") = {4,5};
Physical Surface("Ice") = {7};
