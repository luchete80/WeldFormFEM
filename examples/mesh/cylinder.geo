//GMSH project
lc = 1.0e-3;
ld = 0.5e-3;
thck = 3.0e-4;
R2 = 0.025;
rn = 5.0e-4;
h1 = 0.0125;
h2 = 0.0175;
hend = 0.01;
lcal = 0.025; //calibrated length
lth1 = 0.03; //Total length of narrow section

Point(1)  = {0.0,0.0   ,0.0,lc}; 
Point(2)  = {0.001,0.0   ,0.0,lc}; 
Point(3)  = {0.0,  h1-rn,0.0,lc};
Point(4)  = {0.001,h1-rn,0.0,lc};
Point(5)  = {0.0,h1,0.0,lc};
Point(6)  = {rn, h1+rn,0.0,lc};
Point(7)  = {0.0 ,  h1+rn,0.0,lc};
Point(8)  = {0.001 ,h1+rn,0.0,lc};

Point(9)  = {0.005, 0.0, 0.0, lc};
Point(10) = {0.005, h1+rn, 0.0, lc};
Point(11) = {0.005, h1-rn, 0.0, lc};

Point(12) = {lcal, 0.0, 0.0, lc};
Point(13) = {lcal, h1+rn, 0.0, lc};

Point(14) = {lth1, 0.0,   0.0, lc};
Point(15) = {lth1, h1 + rn,     0.0, lc};
Point(16) = {lth1, h1 + rn + R2,0.0, lc};



a2 = Acos ((h1 + rn + R2-h2)/R2); //Angle of the curve of probe clam area
test = (h1 + rn + R2-h2);

Printf("angle: %f ", a2);
Printf("test: %f",test);

Point(17) = {lth1+ R2*Sin(a2) , 0.0, 0.0, lc};
Point(18) = {lth1+ R2*Sin(a2) , h2,0.0};

Point(19) = {lth1+ R2*Sin(a2)+hend, 0.0,0.0,lc};
Point(20) = {lth1+ R2*Sin(a2)+hend, h2,0.0,lc};

// From https://chi-tech.github.io/d4/db9/_gmsh_example_01.html
//The "Transfinite Line" command is used. This command specifies opposing faces of the four sided
// shapes we will apply a structured mesh to and also how many grid points will be used when drawing mesh lines between the two defined faces. 
//"Using Progression 1" indicates that the structured mesh should have uniform grading. 
//Note that in the outer region of the geometry, "Using Progression 0.85" is used. 
//This causes the meshing to be graded. Also note the input "Transfinite Line {-77, 79}" where the negative sign ensures that both lines have the same orientation which is required when using a grading. Test the input without this and see that the grading is flipped on one side.

//After the transfinite lines are specified, each surface is specified to be mesh with a structured mesh by the syntax "Transfinite Surface {i}". Lastly, "Recombine Surface "*";" is included to specify that the mesh be made up of quadrilaterials. If this input is not included, then gmsh will create a mesh of structured triangles.

Line(1)  = {1,2};
Line(2)  = {2,4};
Line(3)  = {4,3};
Line(4)  = {3,1};

Transfinite Line {1,3} = 10;
Transfinite Line {2,4} = 30;

Line(5)   = {3,5};
Circle(6) = {5, 7, 6};
Line(7)   = {6,8};
Line(8)   = {8,4};

Transfinite Line {5,6,7,8} = 10;


Line(9)   = {2,9};
Line(10)  = {9,11};
Line(11)  = {11,10};
Line(12)  = {10,8};
Line(13)  = {11,4};

Line(14) = {9,12};
Line(15) = {12,13};
Line(16) = {13,10};

//calibrated zone

Circle(17) = {15, 16, 18}; //Clamp Zone
Line(18) = {12,14};
Line(19) = {14,15};
Line(20) = {15,13};
Line(21) = {14,17};
Line(22) = {17,18};

Line(23) = {17,19};
Line(24) = {18,20};
Line(25) = {19,20};

Transfinite Line {8} = 8;

Transfinite Line {11} = 3;


Transfinite Line {10} = 30; //As  2and 4, Top next to znotch

Transfinite Line {12,13,9} = 10; //bottom

//Transfinite Line {12,13} = 10; //bottom

Curve Loop(1) = {1,2,3,4};
Surface(1) = {1};

Transfinite Surface{1};

Curve Loop(2) = {5,6,7,8,3};
Plane Surface(2) = {2};

//Transfinite Surface{2};

//Curve Loop(3) = {9,10,11,12,8,-2};
//Plane Surface(3) = {3};
Curve Loop(3) = {-2,9,10,13};
Surface(3) = {3};
Transfinite Surface{3};

Curve Loop(4) = {-13,11,12,8};
Surface(4) = {4};


Curve Loop(5) = {-11,-10,14,15,16};
Plane Surface(5) = {5};

Transfinite Line {14,16} = 40;
Transfinite Line {15,19} = 10;

Transfinite Line {18,20} = 4;
Curve Loop(6) = {-15,18,19,20};
Plane Surface(6) = {6};

Transfinite Line {15,19,22,25} = 10;
Transfinite Line {17,21} = 10;

Transfinite Line {23,24} = 10; //Horizontal, clamp area

Transfinite Line {23,24} = 10; //Vertical, clamp area


Curve Loop(7) = {-19,21,22,-17};
Plane Surface(7) = {7};

Curve Loop(8) = {23,25,-24,-22};
Plane Surface(8) = {8};

Transfinite Surface{6,7,8};

Recombine Surface "*";


Extrude {0,0,thck} {Surface{1}; Layers{2};Recombine ;}
Extrude {0,0,thck} {Surface{2}; Layers{2};Recombine ;}
Extrude {0,0,thck} {Surface{3}; Layers{2};Recombine ;}
Extrude {0,0,thck} {Surface{4}; Layers{2};Recombine ;}
Extrude {0,0,thck} {Surface{5}; Layers{2};Recombine ;}
Extrude {0,0,thck} {Surface{6}; Layers{2};Recombine ;}
Extrude {0,0,thck} {Surface{7}; Layers{2};Recombine ;}
Extrude {0,0,thck} {Surface{8}; Layers{2};Recombine ;}

// /*//+*/
// Show "*";

/*//+*/
Physical Surface("left",    1) = {46,57};
Physical Surface("bottom",  2) = {34,87,136,158,180,198};
Physical Surface("back",    3) = {1,2,3,4,5,6,7,8};
Physical Surface("right",   4) = {202};

Physical Volume("vol1", 1) = {1};
Physical Volume("vol2", 2) = {2};
Physical Volume("vol3", 3) = {3};
Physical Volume("vol4", 4) = {4};
Physical Volume("vol5", 5) = {5};
Physical Volume("vol6", 6) = {6};
Physical Volume("vol7", 7) = {7};
Physical Volume("vol8", 8) = {8};
//+

