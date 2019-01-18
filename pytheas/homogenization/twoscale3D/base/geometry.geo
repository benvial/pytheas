
Include "parameters.dat";
lc_host      = l_carac/(parmesh*Sqrt[Fabs[eps_host_re]]);
lc_incl	     = l_carac/(parmesh_incl*Sqrt[Fabs[eps_incl_re]]);


Point(1) = {-dx/2, -dy/2, -dz/2, lc_host};
Point(2) = {dx/2, -dy/2, -dz/2, lc_host};
Point(3) = {dx/2, dy/2, -dz/2, lc_host};
Point(4) = {-dx/2, dy/2, -dz/2, lc_host};
Point(5) = {-dx/2, -dy/2, dz/2, lc_host};
Point(6) = {dx/2, -dy/2, dz/2, lc_host};
Point(7) = {dx/2, dy/2, dz/2, lc_host};
Point(8) = {-dx/2, dy/2, dz/2, lc_host};
Line (1) = {5, 6};
Line (2) = {6, 2};
Line (3) = {2, 1};
Line (4) = {1, 5};
Line (5) = {8, 7};
Line (6) = {7, 3};
Line (7) = {3, 4};
Line (8) = {4, 8};
Line (9) = {5, 8};
Line (10) = {6, 7};
Line (11) = {2, 3};
Line (12) = {1, 4};
Line Loop (1000014) = {5, -10, -1, 9};
Plane Surface (14) = {1000014};
Line Loop (1000016) = {8, 5, 6, 7};
Plane Surface (16) = {1000016};
Line Loop (1000018) = {4, 9, -8, -12};
Plane Surface (18) = {1000018};
Line Loop (1000020) = {2, 11, -6, -10};
Plane Surface (20) = {1000020};
Line Loop (1000022) = {3, 12, -7, -11};
Plane Surface (22) = {1000022};
Line Loop (1000024) = {4, 1, 2, 3};
Plane Surface (24) = {1000024};
Surface Loop (102) = {14, 16, 18, 24, 20, 22};


/* #################################### */


Periodic Line {1} = {5};
Periodic Line {3} = {7};
Periodic Line {9} = {10};
Periodic Line {12} = {11};
Periodic Line {4} = {2} ;
Periodic Line {4} = {8} ;
Periodic Line {1} = {5} ;
Periodic Line {3} = {1} ;
Periodic Line {8} = {6} ;
Periodic Line {7} = {5} ;
Periodic Line {12} = {9} ;
Periodic Line {11} = {10} ;



Periodic Surface 18 {9,8,12,4} = 20 {10,-6,11,-2};
Periodic Surface 24 {1,2,3,4} = 16 {5,6,7,8};
Periodic Surface  22 {11,-7,12,-3} = 14 {10,5,9,1};


/*---------------------------------------------------*/

If (inclusion_flag)
  Point(10) = {0,0,0,lc_incl};
  Point(11) = {ax,0,0,lc_incl};
  Point(12) = {0,ay,0,lc_incl};
  Point(13) = {0,0,az,lc_incl};
  Point(14) = {-ax,0,0,lc_incl};
  Point(15) = {0,-ay,0,lc_incl};
  Point(16) = {0,0,-az,lc_incl};


  Ellipsis(21) = {11,10,11,12};
  Ellipsis(22) = {12,10,13,13};
  Ellipsis(23) = {13,10,11,11};
  Ellipsis(24) = {16,10,11,11};
  Ellipsis(25) = {12,10,13,16};
  Ellipsis(26) = {14,10,11,12};
  Ellipsis(27) = {16,10,11,14};
  Ellipsis(28) = {14,10,11,15};
  Ellipsis(29) = {13,10,11,14};
  Ellipsis(30) = {15,10,13,13};
  Ellipsis(31) = {11,10,11,15};
  Ellipsis(32) = {15,10,13,16};


  Line Loop(1) = {21,22,23};   Surface(1) = {1};
  Line Loop(2) = {21,24,25};   Surface(2) = {2};
  Line Loop(3) = {25,27,26};   Surface(3) = {3};
  Line Loop(4) = {26,22,29};   Surface(4) = {4};
  Line Loop(5) = {31,30,23}; Surface(5) = {5};
  Line Loop(6) = {31,32,24}; Surface(6) = {6};
  Line Loop(7) = {32,27,28};  Surface(7) = {7};
  Line Loop(8) = {28,30,29};  Surface(8) = {8};

  Surface Loop(101) = {1,2,3,4,5,6,7,8};

EndIf

/*---------------------------------------------------*/


If (inclusion_flag)

  Volume (201) = {101};/* ellipsoïd */
  Volume (202) = {101,102};/* host */
Else
  Volume (202) = {102};
EndIf

If (inclusion_flag)
  Physical Volume(1001)={201};
EndIf

Physical Volume(1002)={202};

/* Physical Surface(10)={14:24:2};
Physical Surface(500)={1,2,3,4,5,6,7,8}; */


Physical Surface(501)={24};
Physical Surface(502)={16};
Physical Surface(503)={14};
Physical Surface(504)={22};
Physical Surface(505)={18};
Physical Surface(506)={20};

Physical Point(10000) = {1};



Coherence;
Coherence;
Coherence;
Coherence;


/*---------------------------------------------------*/

/*

Mesh.Algorithm   = 1; // // 1=MeshAdapt, 5=Delaunay, 6=Frontal
Mesh.Algorithm3D = 1; // // 1=Delaunay, 4=Frontal*/
