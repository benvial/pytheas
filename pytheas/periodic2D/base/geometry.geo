Include "parameters.dat";
/* Geometry.Tolerance = 1e-7 ; */
lc_des          = lambda_mesh/(parmesh_des*Sqrt[Fabs[eps_des_re]]);
lc_layer1       = lambda_mesh/(parmesh*Sqrt[Fabs[eps_layer1_re]]);
lc_layer2       = lambda_mesh/(parmesh*Sqrt[Fabs[eps_layer2_re]]);
lc_sub          = lambda_mesh/(parmesh*Sqrt[Fabs[eps_sub_re]]);
lc_sup          = lambda_mesh/(parmesh*Sqrt[Fabs[eps_sup_re]]);
lc_pmlbot       = lambda_mesh/(parmesh_pml*Sqrt[Fabs[eps_sub_re]]);
lc_pmltop       = lambda_mesh/(parmesh_pml*Sqrt[Fabs[eps_sup_re]]);
lc_incl       = lambda_mesh/(parmesh_incl*Sqrt[Fabs[eps_incl_re]]);

/*lc_des          = lambda_mesh/(parmesh_des);
lc_layer1       = lambda_mesh/(parmesh);
lc_layer2       = lambda_mesh/(parmesh);
lc_sub          = lambda_mesh/(parmesh);
lc_sup          = lambda_mesh/(parmesh);
lc_pmlbot       = lambda_mesh/(parmesh_pml);
lc_pmltop       = lambda_mesh/(parmesh_pml);
*/


/*If (inclusion_flag)
  Include "inclusion.geo";
EndIf*/

Point(1)  = {-d/2.,-h_sub-h_pmlbot, 0. , lc_pmlbot};
Point(2)  = {-d/2.,-h_sub        , 0. , lc_sub};
Point(3)  = {-d/2., 0.        , 0. , lc_layer1};
Point(4)  = {-d/2., h_layer1         , 0. , lc_des};
Point(5)  = {-d/2., h_des  + h_layer1         , 0. , lc_des};
Point(6)  = {-d/2., h_des  + h_layer1+ h_layer2          , 0. , lc_layer2};
Point(7)  = {-d/2., h_des  + h_layer1+ h_layer2 + h_sup         , 0. , lc_sup};
Point(8)  = {-d/2., h_des  + h_layer1+ h_layer2 + h_sup +h_pmltop    , 0. , lc_pmltop};


Curve(1) = {1, 2};
Curve(2) = {2, 3};
Curve(3) = {3, 4};
Curve(4) = {4, 5};
Curve(5) = {5, 6};
Curve(6) = {6, 7};
Curve(7) = {7, 8};


If (quad_mesh_flag)
  out[] = Extrude{d,0,0}{Curve{1};Curve{2};Curve{3};Layers{d/lc_des};};
  out[] = Extrude{d,0,0}{Curve{4};Layers{d/lc_des};Recombine;};
  out[] = Extrude{d,0,0}{Curve{5};Curve{6};Curve{7};Layers{d/lc_des};};
  /*out[] = Extrude{d,0,0}{Curve{1};Curve{2};Curve{3};};
  out[] = Extrude{d,0,0}{Curve{4};Recombine;};
  out[] = Extrude{d,0,0}{Curve{5};Curve{6};Curve{7};};*/
  Else
    /* If (extrude_mesh_flag) */
      /* out[] = Extrude{d,0,0}{Curve{1};Curve{2};Curve{3};Curve{4};Curve{5};Curve{6};Curve{7};Layers{d/lc_des};}; */
    /* Else */
      out[] = Extrude{d,0,0}{Curve{1};Curve{2};Curve{3};Curve{4};Curve{5};Curve{6};Curve{7};};
    /* EndIf */


    /* Periodic Curve {1} = {8};
    Periodic Curve {2} = {12};
    Periodic Curve {3} = {16};
    Periodic Curve {4} = {20};
    Periodic Curve {5} = {24};
    Periodic Curve {6} = {28};
    Periodic Curve {7} = {32}; */

EndIf

tag_des = 23;
If (inclusion_flag)
  tag_des = 730;
  Include "inclusion.geo";
  Curve Loop(721) = {4, 22, -20, -18};
  Curve Loop(722) = {1000};
  Plane Surface(tag_des) = {721, 722};
  Plane Surface(731) = {722};
EndIf



Physical Curve(101) = {1, 2, 3, 4, 5, 6, 7};          // Bloch_LeftX-
Physical Curve(102) = {8, 12, 16, 20, 24, 28, 32};    // Bloch_RightX+
Physical Curve(110) = {9, 34};                        // Dirichlet
Physical Curve(120) = {10,14,18,22,26,30};            // Continuity

Physical Surface(1000) = {11};        // PML_bot
Physical Surface(2000) = {15};        // sub
Physical Surface(3000) = {19};        // layer 1
Physical Surface(4000) = {tag_des};        // layer des
Physical Surface(5000) = {27};        // layer 2
Physical Surface(6000) = {31};        // sup
Physical Surface(7000) = {35};        // pmltop
If (inclusion_flag)
  Physical Surface(8000) = {731};        // inclusion
EndIf

Physical Point(10000) = {1};        // PrintPoint



Coherence;
Coherence;
Coherence;
Coherence;
Coherence;
Coherence;
Coherence Mesh;
/* Mesh.Algorithm = 6; // Frontal */
