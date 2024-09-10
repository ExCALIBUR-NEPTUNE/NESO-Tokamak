SetFactory("OpenCASCADE");
Merge "scaled_reactor.step";
//+
Physical Volume("fluid", 45) = {1};
//+
Physical Surface("outer_wall", 46) = {3, 2, 1, 17, 16, 15, 12, 10, 5, 6, 7, 8, 9, 14, 4, 11, 13};
//+
Wire(19) = {43};
Extrude { Curve{8}; Curve{44}; } Using Wire {19}

//+
Wire(22) = {43};
Extrude { Curve{8}; Curve{44}; } Using Wire {22}

//+
Wire(25) = {43};
Extrude { Curve{8}; Curve{44}; } Using Wire {25}

