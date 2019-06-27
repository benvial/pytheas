
Include "parameters.dat";

l_carac = 0.5*(dx+dy);

lc_host      = l_carac/(parmesh*Sqrt[Fabs[eps_host_re]]);
lc_incl	     = l_carac/(parmesh_incl*Sqrt[Fabs[eps_incl_re]]);


Point(1)  = {-dx/2.,-dy/2, 0. , lc_host};
Point(2)  = {-dx/2.,dy/2, 0. , lc_host};
Point(3)  = {dx/2.,dy/2, 0. , lc_host};
Point(4)  = {dx/2.,-dy/2, 0. , lc_host};


Point(5)  = {0.,0., 0. , lc_incl};

Line(11) = {1, 2};

/* Periodic Line {11} = {13};
Periodic Line {12} = {14}; */


If (quad_mesh_flag)
/* Periodic Line {11} = {12}; */
/* Periodic Line {13} = {14}; */
  out[] = Extrude{dx,0,0}{Line{11};Layers{dx/lc_incl};Recombine;};
  tag_des =  out[1];
  Line Loop(tag_des) = {11, -12, -13, 14};

Else
  tag_des = 20;
  Line(12) = {2, 3};
  Line(13) = {3, 4};
  Line(14) = {4, 1};
  Line Loop(tag_des) = {11, 12, 13, 14};
EndIf





tag_host = 100;
If (inclusion_flag)
  tag_incl =101;
  Include "inclusion.geo";
  Line Loop(30) = {1000};
  Plane Surface(tag_incl) = {30};
  Plane Surface(tag_host) = {tag_des, 30};
Else
  Plane Surface(tag_host) = {tag_des};
EndIf

If (quad_mesh_flag)
  Physical Line(101) = {11};		        // Bloch_LeftX-
  Physical Line(103) = {12};		        // Bloch_RightX+
  Physical Line(102) = {14};		        // Bloch_TopY+
  Physical Line(104) = {13};		        // Bloch_BotY-
Else
  Physical Line(101) = {11};		        // Bloch_LeftX-
  Physical Line(103) = {13};		        // Bloch_RightX+
  Physical Line(102) = {12};		        // Bloch_TopY+
  Physical Line(104) = {14};		        // Bloch_BotY-
EndIf

Physical Surface(1000) = {tag_host};        // host
If (inclusion_flag)
  Physical Surface(2000) = {tag_incl};        // inclusion
  Physical Line(100) = {1000};		        // surface inclusion
EndIf


Physical Point(10000) = {5};        // Printpoint


Coherence;Coherence;Coherence;Coherence;

/* Geometry.Tolerance=1e-3; */

/*
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
Coherence; */


/*

##################################

#####################
Include "parameters.dat";

l_carac = 0.5*(dx+dy);

lc_host      = l_carac/(parmesh*Sqrt[Fabs[eps_host_re]]);
lc_incl	     = l_carac/(parmesh_incl*Sqrt[Fabs[eps_incl_re]]);


Point(1)  = {-dx/2.,-dy/2, 0. , lc_host};
Point(2)  = {-dx/2.,dy/2, 0. , lc_host};
Point(3)  = {dx/2.,dy/2, 0. , lc_host};
Point(4)  = {dx/2.,-dy/2, 0. , lc_host};

Point(5)  = {0.,0., 0. , lc_incl};

Line(11) = {1, 2};
Line(12) = {2, 3};
Line(13) = {3, 4};
Line(14) = {4, 1};


Line Loop(20) = {11, 12, 13, 14};


tag_host = 100;
If (inclusion_flag)
  tag_incl =101;
  Include "inclusion.geo";
  Line Loop(30) = {1000};
  Plane Surface(tag_incl) = {30};
  Plane Surface(tag_host) = {20, 30};
Else

  Plane Surface(tag_host) = {20};
EndIf


Physical Line(101) = {11};		        // Bloch_LeftX-
Physical Line(103) = {13};		        // Bloch_RightX+
Physical Line(102) = {12};		        // Bloch_TopY+
Physical Line(104) = {14};		        // Bloch_BotY-



Physical Surface(1000) = {tag_host};        // host
If (inclusion_flag)
  Physical Surface(2000) = {tag_incl};        // inclusion
EndIf


Physical Point(10000) = {5};        // Printpoint

Coherence;Coherence;Coherence;Coherence;Coherence; */
