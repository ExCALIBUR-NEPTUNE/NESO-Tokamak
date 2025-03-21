// geo file for meshing with Gmsh meshing software created by FreeCAD

// enable multi-core processing
General.NumThreads = 6;
// open brep geometry
Merge "Face001_Geometry.brep";

// Characteristic Length
// no boundary layer settings for this mesh
// min, max Characteristic Length
Mesh.CharacteristicLengthMax = 50.0;
Mesh.CharacteristicLengthMin = 0.0;
Mesh.MeshSizeFromCurvature = 12; // number of elements per 2*pi radians, 0 to deactivate

// optimize the mesh
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 0;
// High-order meshes optimization (0=none, 1=optimization, 2=elastic+optimization, 3=elastic, 4=fast curving)
Mesh.HighOrderOptimize = 0;

// mesh order
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0; // Second order nodes are created by linear interpolation instead by curvilinear

// mesh algorithm, only a few algorithms are usable with 3D boundary layer generation
// 2D mesh algorithm (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=BAMG, 8=DelQuad, 9=Packing Parallelograms)
Mesh.Algorithm = 2;
// 3D mesh algorithm (1=Delaunay, 2=New Delaunay, 4=Frontal, 7=MMG3D, 9=R-tree, 10=HTX)
Mesh.Algorithm3D = 1;

// meshing
Geometry.Tolerance = 1e-06; // set geometrical tolerance (also used for merging nodes)
Mesh  2;
Coherence Mesh; // Remove duplicate vertices

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79};
Physical Surface(1) = {1};
Physical Curve(100) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79};

// save
// Ignore Physical definitions and save all elements;
// Mesh.SaveAll = 1;
// Save "Face001_Mesh.unv";



// **********************************************************************
// Gmsh documentation:
// https://gmsh.info/doc/texinfo/gmsh.html#Mesh
//
// We do not check if something went wrong, like negative jacobians etc. You can run Gmsh manually yourself: 
//
// to see full Gmsh log, run in bash:
// C:/Program Files/FreeCAD 0.21/bin/gmsh.exe - C:\Users\jedge\AppData\Local\Temp\fcfem_9mex2f_3\shape2mesh.geo
//
// to run Gmsh and keep file in Gmsh GUI (with log), run in bash:
// C:/Program Files/FreeCAD 0.21/bin/gmsh.exe C:\Users\jedge\AppData\Local\Temp\fcfem_9mex2f_3\shape2mesh.geo
