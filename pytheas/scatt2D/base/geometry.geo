Include "parameters.dat";
/*Geometry.Tolerance = 1e-4 ;*/
lc_des          = lambda_mesh/(parmesh_des*Sqrt[Fabs[eps_des_re]]);
lc_host         = lambda_mesh/(parmesh*Sqrt[Fabs[eps_host_re]]);
lc_pml          = lambda_mesh/(parmesh_pml*Sqrt[Fabs[eps_host_re]]);
lc_incl         = lambda_mesh/(parmesh_incl*Sqrt[Fabs[eps_incl_re]]);


Point(1) = {domX_L, domY_B, 0, lc_host};
Point(2) = {domX_R, domY_B, 0, lc_host};
Point(3) = {domX_R, domY_T, 0, lc_host};
Point(4) = {domX_L, domY_T, 0, lc_host};
Point(5) = {domX_L-h_pml, domY_B-h_pml, 0, lc_pml};
Point(6) = {domX_R+h_pml, domY_B-h_pml, 0, lc_pml};
Point(7) = {domX_R+h_pml, domY_T+h_pml, 0, lc_pml};
Point(8) = {domX_L-h_pml, domY_T+h_pml, 0, lc_pml};

Point(9) = {domX_L, domY_B-h_pml, 0, lc_pml};
Point(10)= {domX_L-h_pml, domY_B, 0, lc_pml};
Point(11)= {domX_R, domY_B-h_pml, 0, lc_pml};
Point(12)= {domX_R+h_pml, domY_B, 0, lc_pml};
Point(13)= {domX_R, domY_T+h_pml, 0, lc_pml};
Point(14)= {domX_R+h_pml, domY_T, 0, lc_pml};
Point(15)= {domX_L, domY_T+h_pml, 0, lc_pml};
Point(16)= {domX_L-h_pml, domY_T, 0, lc_pml};

Point(17)= {-hx_des/2,-hy_des/2,0.0,lc_des};
Point(18)= {-hx_des/2, hy_des/2,0.0,lc_des};
/* Point(19)= { hx_des/2, hy_des/2,0.0,lc_des};
Point(20)= { hx_des/2,-hy_des/2,0.0,lc_des}; */


/* Point(51)= {x_target, y_target, 0.0,lc_host}; */

z_target=0;
RadiusSourceIntX = r_target;
RadiusSourceIntY = r_target;
lcSourceInt=lc_host/1;
/* lcSourceExt=lc_host; */
PS = newp; Point(PS) = {x_target, y_target, z_target, lcSourceInt};

PSourceInt1 = newp; Point(PSourceInt1) = {x_target + RadiusSourceIntX, y_target ,z_target, lcSourceInt};
PSourceInt2 = newp; Point(PSourceInt2) = {x_target, y_target + RadiusSourceIntY, z_target, lcSourceInt};
PSourceInt3 = newp; Point(PSourceInt3) = {x_target - RadiusSourceIntX, y_target, z_target, lcSourceInt};
PSourceInt4 = newp; Point(PSourceInt4) = {x_target, y_target - RadiusSourceIntY, z_target, lcSourceInt};




Curve(1) = {1, 4};
Curve(2) = {4, 3};
Curve(3) = {3, 2};
Curve(4) = {2, 1};
Curve(5) = {9, 5};
Curve(6) = {5, 10};
Curve(7) = {10, 16};
Curve(8) = {16, 8};
Curve(9) = {8, 15};
Curve(10) = {15, 13};
Curve(11) = {13, 7};
Curve(12) = {7, 14};
Curve(13) = {14, 12};
Curve(14) = {12, 6};
Curve(15) = {6, 11};
Curve(16) = {11, 9};
Curve(17) = {9, 1};
Curve(18) = {1, 10};
Curve(19) = {16, 4};
Curve(20) = {4, 15};
Curve(21) = {13, 3};
Curve(22) = {3, 14};
Curve(23) = {11, 2};
Curve(24) = {2, 12};
/* Curve(25) = {17, 20};
Curve(26) = {20, 19};
Curve(27) = {19, 18}; */
Curve(28) = {18, 17};

Curve Loop(30) = {6, -18, -17, 5};
Plane Surface(30) = {30};
Curve Loop(32) = {17, -4, -23, 16};
Plane Surface(32) = {32};
Curve Loop(34) = {23, 24, 14, 15};
Plane Surface(34) = {34};
Curve Loop(36) = {24, -13, -22, 3};
Plane Surface(36) = {36};
Curve Loop(38) = {22, -12, -11, 21};
Plane Surface(38) = {38};
Curve Loop(40) = {21, -2, 20, 10};
Plane Surface(40) = {40};
Curve Loop(42) = {9, -20, -19, 8};
Plane Surface(42) = {42};
Curve Loop(44) = {7, 19, -1, 18};
Plane Surface(44) = {44};



If (quad_mesh_flag)
  out[] = Extrude{hx_des,0,0}{Curve{28};Layers{hx_des/lc_des};Recombine;};
Else
  /* out[] = Extrude{hx_des,0,0}{Curve{28};Layers{hx_des/lc_des};}; */
  out[] = Extrude{hx_des,0,0}{Curve{28};};
EndIf

tag_des =  out[1];
Curve Loop(tag_des) = {-28, 46, 45, -47}; //box

If (inclusion_flag)
loopholes={};
loopholes[0] = tag_des;
  For k In {1:nb_incl}
  tag_incl=1000*(k);
  /* name = Sprintf("ellipse", k); */
    Include Sprintf("ellipse%g.geo", k-1);
    Curve Loop(tag_incl) = {tag_incl};
    Plane Surface(tag_incl) = {tag_incl};
    loopholes[k]=tag_incl;
    /* Plane Surface(tag_des) = {tag_des, tag_incl}; */
  EndFor

  Plane Surface(tag_des+1) = loopholes[];

Else
  Plane Surface(tag_des+1) = {tag_des};
  /*Transfinite Surface(tag_des) = {tag_des};*/
EndIf

Curve Loop(47) = {1,2,3,4}; //host


Physical Surface(1) = {30,34,38,42}; //PML corners
Physical Surface(2) = {32,40};       //PML top bot
Physical Surface(3) = {44,36};       //PML left right
Physical Surface(4) = {47};  //host
Physical Surface(5) = {tag_des+1};  //box

If (inclusion_flag)
  /* Physical Surface(6) = {tag_incl};  // inclusion */
  physloop={};
  For k In {0:nb_incl-1:1}
  	physloop[k]=1000*(k+1);
  EndFor
  Physical Surface(6) = physloop[];    //holes
  /* Physical Curve(100) = {1000};  // incl Curve */
EndIf



If (target_flag)
  LSourceInt1 = newreg; Ellipse(LSourceInt1) = {PSourceInt1, PS, PSourceInt1, PSourceInt2};
  LSourceInt2 = newreg; Ellipse(LSourceInt2) = {PSourceInt2, PS, PSourceInt1, PSourceInt3};
  LSourceInt3 = newreg; Ellipse(LSourceInt3) = {PSourceInt3, PS, PSourceInt1, PSourceInt4};
  LSourceInt4 = newreg; Ellipse(LSourceInt4) = {PSourceInt4, PS, PSourceInt1, PSourceInt1};

  LLSourceInt = newreg; Curve Loop(LLSourceInt) = {LSourceInt1, LSourceInt2, LSourceInt3, LSourceInt4};

  SSourceExt = newreg; Plane Surface(SSourceExt) = {LLSourceInt};


  Plane Surface(47) = {47, tag_des, LLSourceInt};

  Physical Surface(7) = {SSourceExt};
Else
    Plane Surface(47) = {47, tag_des};
EndIf

Physical Curve(100) = {25};  //box bot
Physical Curve(200) = {26};  //box right
Physical Curve(300) = {27};  //box top
Physical Curve(400) = {28};  //box left

Physical Point(10000) = {1};        // PrintPoint
Physical Point(20000) = {PS};        // TargetPoint

Coherence;Coherence;Coherence;Coherence;Coherence;
/* Coherence; */
/* Physical Curve(200) = {46};  //box Curve */
/* Physical Point(10000) = {1};        // Printpoint
Physical Point(20000) = {25};        // Sourcepoint */


/*Mesh.RemeshAlgorithm = 1; // automatic
Mesh.RemeshParametrization = 7; // conformal finite element
Mesh.Algorithm = 6; // Frontal*/
