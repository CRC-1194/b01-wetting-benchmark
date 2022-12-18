// Define sphere using two circle arcs of 90 degree 
// Use center at (4.0 4.0 0.1) as in test cases of Francois
// Parameters to Determine center and radius
x = 0.0025;
y = 0.0025;
z = 0.00068578;
r = 0.001;
// Parameters to control mesh size
h = 0.00001;
// Define support points
// Center point
Point(1) = {x, y, z, h};
// Support points in plane parallel to x-y-plane through center
Point(2) = {x-r, y, z, h};
Point(3) = {x, y-r, z, h};
Point(4) = {x+r, y, z, h};
Point(5) = {x, y+r, z, h};
// Support points in plane parallel to x-z-plane through center
Point(6) = {x, y, z-r, h};
Point(7) = {x, y, z+r, h};
// Define circle arcs in plane parallel to x-y-plane
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
// Define circle arcs in plane parallel to x-z-plane
Circle(5) = {6, 1, 4};
Circle(6) = {4, 1, 7};
Circle(7) = {7, 1, 2};
Circle(8) = {2, 1, 6};
// Define circle arcs in plane parallel to y-z-plane
Circle(9) = {6, 1, 3};
Circle(10) = {3, 1, 7};
Circle(11) = {7, 1, 5};
Circle(12) = {5, 1, 6};
// Define Sphere segments' boundaries:
// Begin in positive sector, iterate counter clockwise, segments with x >= 0
Line Loop(1) = {3, -11, -6};
Line Loop(2) = {3, 12, 5};
Line Loop(3) = {-2, -9, 5};
Line Loop(4) = {-2, 10, -6};
// Begin in sector (x < 0, y > 0, z > 0) and iterate counter clockwise
Line Loop(5) = {-4, -11, 7};
Line Loop(6) = {-4, 12, -8};
Line Loop(7) = {1, -9, -8};
Line Loop(8) = {1, 10, 7};
// Define sphere surface patches
Ruled Surface(1) = {1} In Sphere {1};
Ruled Surface(2) = {2} In Sphere {1};
Ruled Surface(3) = {3} In Sphere {1};
Ruled Surface(4) = {4} In Sphere {1};
Ruled Surface(5) = {5} In Sphere {1};
Ruled Surface(6) = {6} In Sphere {1};
Ruled Surface(7) = {7} In Sphere {1};
Ruled Surface(8) = {8} In Sphere {1};
// Test to orient normal vectors consistently
// They should now point outwards
Reverse Surface {2};
Reverse Surface {4};
Reverse Surface {5};
Reverse Surface {7};
// Combine to sphere surface
Physical Surface(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Mesh.RemeshAlgorithm = 1; // Choose automatically
Mesh.RemeshParametrization = 0; // harmonic
Mesh.Algorithm = 6; // Frontal

