
Include "parameters.dat";

lc_host      = lambda0/(parmesh);
lc_rod	     = lambda0/(parmesh);


Point(1)  = {-dx/2.,-dy/2, 0. , lc_host};
Point(2)  = {-dx/2.,dy/2, 0. , lc_host};
Point(3)  = {dx/2.,dy/2, 0. , lc_host};
Point(4)  = {dx/2.,-dy/2, 0. , lc_host};


Point(5)  = {0.,0., 0. , lc_rod};

Line(11) = {1, 2};
Line(12) = {2, 3};
Line(13) = {3, 4};
Line(14) = {4, 1};
/* Periodic Line {11} = {13};
Periodic Line {12} = {14}; */

Line Loop(20) = {11, 12, 13, 14};


Plane Surface(30) = {20};

Physical Line(101) = {11};		        // Bloch_LeftX-
Physical Line(103) = {13};		        // Bloch_RightX+
Physical Line(102) = {12};		        // Bloch_TopY+
Physical Line(104) = {14};		        // Bloch_BotY-



Physical Surface(1000) = {30};        // host


Physical Point(10000) = {5};        // Printpoint

Coherence;
Coherence;
Coherence;
Coherence;
Coherence;
Coherence;
