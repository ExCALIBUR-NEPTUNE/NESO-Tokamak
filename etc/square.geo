//Create construction points for the corners of the mesh
//Point(i) = (x, y, z, scale)
Point(1) = {0.5, -0.5, 0, 1.0};
Point(2) = {0.5, 0.5, 0, 1.0};
Point(3) = {1.5, 0.5, 0, 1.0};
Point(4) = {1.5, -0.5, 0, 1.0};
//Point(5) = {2.5, 0, 0, 1.0};

//Create construction lines for the edges of the mesh
//Line(i) = (Start_i, End_i) 
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//Line(5) = {5, 3};
//Line(6) = {4, 5};

//Add physical lines from the construction lines
Physical Line(5) = {1};
Physical Line(6) = {2};
Physical Line(7) = {3};
Physical Line(8) = {4};
//Physical Line(9) = {5};
//Physical Line(10) = {6};

//Close the lines with a curve loop
Curve Loop(1) = {1,2,3,4};
//Curve Loop(2) = {3,6,5};
//Create a construction surface out of the curve loop
Plane Surface(1) = {1};
//Plane Surface(2) = {2};
//Create a physical surface out of the construction surface
Physical Surface(1) = {1};
//Physical Surface(2) = {2};

//A transfinite line means that 65 points are created along lines
// 1 and 3, which will form the corners of 64 mesh elements.
//Progression 1 means they are uniformly spaced
Transfinite Line {1, 3} = 9 Using Progression 1;
Transfinite Line {2, 4} = 9 Using Progression 1;
//Set the surface
Transfinite Surface {1};

//The default is tris, this line is necessary to produce quads
Recombine Surface {1};