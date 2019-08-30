Include "parameters.dat";


// #############################################################################
// #############################################################################

Group {
    // Domains
    host         = Region[1000];

    If (inclusion_flag)
      incl          = Region[2000];
      Omega          = Region[{host, incl}];
    Else
      Omega          = Region[{host}];
    EndIf


    // Boundaries
  	SurfBlochLeft   = Region[101];
  	SurfBlochRight  = Region[103];
  	SurfBlochTop   = Region[102];
  	SurfBlochBot  = Region[104];

    Borders =  Region[{SurfBlochLeft, SurfBlochRight, SurfBlochTop, SurfBlochBot}];

  	// Points
  	PrintPoint      = Region[10000];

    Omega1          = ElementsOf [{host}, Not Borders];
}
// #############################################################################

Function{
  I3[]=TensorDiag[1,1,1];
  If (inclusion_flag)
    /* epsilonr[host]         = Complex[eps_host_re,eps_host_im] * I3[]; */
    epsilonr[incl]           = Complex[eps_incl_re,eps_incl_im] * I3[];
    If (aniso)
        epsilonr_xx[]  = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}];
        epsilonr_yy[]  = Complex[ScalarField[XYZ[], 0, 1 ]{2}, ScalarField[XYZ[], 0, 1 ]{3}];
        epsilonr_zz[]  = Complex[ScalarField[XYZ[], 0, 1 ]{4}, ScalarField[XYZ[], 0, 1 ]{5}];

        epsilonr[host] =  TensorDiag[epsilonr_xx[],epsilonr_yy[],epsilonr_zz[]];
      Else
        If (interp)
          epsilonr[host]         = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}] * I3[];
        Else
          epsilonr[host]         = Complex[eps_host_re,eps_host_im] * I3[];
        EndIf
      EndIf
  Else
  If (aniso)
      epsilonr_xx[]  = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}];
      epsilonr_yy[]  = Complex[ScalarField[XYZ[], 0, 1 ]{2}, ScalarField[XYZ[], 0, 1 ]{3}];
      epsilonr_zz[]  = Complex[ScalarField[XYZ[], 0, 1 ]{4}, ScalarField[XYZ[], 0, 1 ]{5}];

      epsilonr[Omega] =  TensorDiag[epsilonr_xx[],epsilonr_yy[],epsilonr_zz[]];
    Else
      epsilonr[Omega]    = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}] * I3[];
    EndIf
  EndIf
  If (aniso)
    istore=6;
  Else
    istore=2;
  EndIf
    epsXX[]=CompXX[epsilonr[]];
    epsYY[]=CompYY[epsilonr[]];
    epsXY[]=CompXY[epsilonr[]];
    epsYX[]=CompYX[epsilonr[]];
    deteps[] = epsXX[]*epsYY[] - epsXY[]*epsYX[];
    /* xsi[Omega] = Tensor[epsXX[], epsXX[],epsXX[]];  */
    /* xsi[Omega] = TensorDiag[1/CompYY[epsilonr[]], 1/CompXX[epsilonr[]],1]; */

    /* xsi[Omega] = Transpose[epsilonr[]]/Det[epsilonr[]]; */

    xsi[Omega] = Transpose[epsilonr[]]/deteps[];

    /* xsi[Omega] = TensorDiag[1/CompYY[epsilonr[]], 1/CompXX[epsilonr[]],111]; */




		If (y_flag)
			E[] = Vector[0,1,0] ;
		Else
			E[] = Vector[1,0,0] ;
		EndIf

    source[]=xsi[] * E[];




  // Topology optimization
/* DefineFunction[sadj]; */

xsi_0_xx=1/7;
xsi_0_xy = 0;
coef_obj[] = 1/xsi_0_xx^2;
int_xsi[] = CompY[xsi[]* ( E[] + $1)];

objective[] = (int_xsi[$1]  - xsi_0_xx)^2 ;

adj_source_int[]  = - 2 *  coef_obj[] *  (CompYY[xsi[]]) * Conj[CompYY[xsi[]] + CompY[xsi[] * $1] - xsi_0_xx];

adj_source_int1[]  = (CompYY[xsi[]]) ;



db_deps[] = 0 * E[];
dA_deps[] = 1;//- TensorDiag[1/epsilonr_xx[]^2,1/epsilonr_yy[]^2,1/epsilonr_zz[]^2];
dEq_deps[] = db_deps[] - dA_deps[] * ($1);



}

// #############################################################################

Constraint {

		{Name Bloch;
		    Case {
              { Region SurfBlochRight; Type LinkCplx ; RegionRef SurfBlochLeft;
                Coefficient Complex[1.0,0.0]; Function Vector[$X-dx,$Y,$Z] ;
              }
              { Region SurfBlochTop; Type LinkCplx ; RegionRef SurfBlochBot;
                Coefficient Complex[1.0,0.0]; Function Vector[$X,$Y-dy,$Z] ;
              }
			 }
		}
}

// #############################################################################

Jacobian {
  { Name JVol ;
    Case {
      { Region All ; Jacobian Vol ; }
    }
  }
  { Name JSur ;
    Case {
      { Region All ; Jacobian Sur ; }
    }
  }
    { Name JLin ;
    Case {
      { Region All ; Jacobian Lin ; }
    }
  }
}

// #############################################################################

Integration {
  { Name Int_1 ;
    Case {
      { Type Gauss ;
        Case {
	  { GeoElement Point       ; NumberOfPoints  1 ; }
	  { GeoElement Line        ; NumberOfPoints  4 ; }
	  { GeoElement Triangle    ; NumberOfPoints  6 ; }
	  { GeoElement Quadrangle  ; NumberOfPoints  7 ; }
	}
      }
    }
  }
}
Integration {
  { Name GaussOnePoint ; Case {
      { Type Gauss ;
        Case {
          { GeoElement Point       ; NumberOfPoints  1; }
          { GeoElement Line        ; NumberOfPoints  4; }
          { GeoElement Triangle    ; NumberOfPoints  6; }
          { GeoElement Quadrangle  ; NumberOfPoints  7; }
          { GeoElement Prism       ; NumberOfPoints  1; }
          { GeoElement Tetrahedron ; NumberOfPoints  1; }
          { GeoElement Hexahedron  ; NumberOfPoints  1; }
          { GeoElement Pyramid     ; NumberOfPoints  1; }
        }
      }
    }
  }
}
// #############################################################################

FunctionSpace {
  { Name Hgrad; Type Form0;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_Node;    Support Region[Omega]; Entity NodesOf[Omega]; }
      { Name sn2; NameOfCoef un2; Function BF_Node_2E; Support Region[Omega]; Entity EdgesOf[Omega]; }
      /* { Name sn3; NameOfCoef un3; Function BF_Node_3E; Support Region[Omega]; Entity EdgesOf[Omega]; } */
   }
    Constraint {
      /*{ NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Dirichlet; }*/
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Bloch; }
      /*{ NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }*/
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Bloch; }
      /* { NameOfCoef un3; EntityType EdgesOf ; NameOfConstraint Bloch; } */
   }
  }

}

// #############################################################################

Formulation {
	    {Name electrostat_scalar; Type FemEquation;
    		Quantity {{ Name u; Type Local; NameOfSpace Hgrad;}}
		Equation {
		        Galerkin { [ xsi[] * Dof{d u} , {d u}];
         		In Omega; Jacobian JVol; Integration Int_1; }
            Galerkin { [ ($Source ? source[] : 0) , {d u}];
            In Omega; Jacobian JVol; Integration Int_1; }

            /* Galerkin { [ ($SourceAdj ? Complex[ScalarField [XYZ[], 0, 1 ]{istore}, ScalarField [XYZ[], 0, 1 ]{istore+1}]  : 0) , {d u}];
            In Omega ;Jacobian JVol; Integration Int_1; } */

            /* Galerkin { [ ($SourceAdj ? Complex[ScalarField [XYZ[], 0, 1 ]{istore}, ScalarField [XYZ[], 0, 1 ]{istore+1}] * Complex[ScalarField [XYZ[], 0, 1 ]{istore+2}, ScalarField [XYZ[], 0, 1 ]{istore+3}]  : 0) ,{d u}]; */
            /* In Omega ;Jacobian JVol; Integration Int_1; } */

            Galerkin { [ ($SourceAdj ? adj_source_int[{d u}] : 0) , {d u}];
            In Omega; Jacobian JVol; Integration Int_1; }


                }
            }

        }
// #############################################################################


Resolution {
  { Name electrostat_scalar;
    System {
      { Name S; NameOfFormulation electrostat_scalar; Type ComplexValue;}
    }
    Operation {
      Evaluate[$Source = 1, $SourceAdj = 0];
      Generate[S] ;  Solve[S] ;SaveSolution[S] ;

      If (adjoint_flag)
        /* PostOperation[postop_solution]; */
        PostOperation[postop_int_objective];
        PostOperation[postop_dEq_deps];
        PostOperation[postop_source_adj];
        Evaluate[$Source = 0, $SourceAdj = 1];
        GenerateRHSGroup[S, Omega]; SolveAgain[S] ; SaveSolution[S] ;
        /* Generate[S] ;  Solve[S] ;SaveSolution[S] ; */
        PostOperation[postop_adjoint];
      EndIf
    }
  }

}

// #############################################################################


PostProcessing {
    { Name postpro; NameOfFormulation electrostat_scalar; NameOfSystem S;
            Quantity {
              { Name epsilonr_xx; Value { Local { [CompXX[ epsilonr[]]] ; In Omega; Jacobian JVol; } } }
              { Name epsilonr_yy; Value { Local { [CompYY[ epsilonr[]]] ; In Omega; Jacobian JVol; } } }
              { Name solution; Value { Local { [ {u}] ; In Omega; Jacobian JVol; } } }
							{ Name vx; Value { Local {[CompX[{d u}]] ; In Omega; Jacobian JVol; } } }
							{ Name vy; Value { Local {[CompY[{d u}]] ; In Omega; Jacobian JVol; } } }
              { Name Ix; Value { Integral { [CompX[xsi[]*{d u}]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name Iy; Value { Integral { [CompY[xsi[]*{d u}]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
							{ Name I_inveps_xx; Value { Integral { [CompXX[xsi[]]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
          		{ Name I_inveps_yy; Value { Integral { [CompYY[xsi[]]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name I_inveps_xy; Value { Integral { [CompXY[xsi[]]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name I_inveps_yx; Value { Integral { [CompYX[xsi[]]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name V; Value { Integral { [1]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name u_adj   ; Value { Local { [ {d u} * ElementVol[] ] ; In Omega; Jacobian JVol; } } }
              /* { Name u_adj   ; Value { Local { [ {u}] ; In Omega; Jacobian JVol; } } } */
              { Name int_objective  ; Value { Integral { [objective[{d u}]] ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name dEq_deps_x   ; Value { Local { [CompX[dEq_deps[{d u}]]] ; In Omega; Jacobian JVol; } } }
              { Name dEq_deps_y   ; Value { Local { [CompY[dEq_deps[{d u}]]] ; In Omega; Jacobian JVol; } } }
              { Name dEq_deps  ; Value { Local { [dEq_deps[CompY[{d u}]]] ; In Omega; Jacobian JVol; } } }
              /* { Name sadj_var   ; Value { Local { [$sadj  ] ; In Omega; Jacobian JVol; } } } */
              /* { Name sadj_scalfield   ; Value { Local { [Complex[ScalarField[XYZ[], 0, 1 ]{istore}, ScalarField[XYZ[], 0, 1 ]{istore+1}]  ] ; In Omega; Jacobian JVol; } } } */
              { Name u_re   ; Value { Local { [ Re[{u}]] ; In Omega; Jacobian JVol; } } }
              { Name u_im  ; Value { Local { [ Im[{u}]] ; In Omega; Jacobian JVol; } } }

              /* { Name sadj_int_re  ; Value { Integral { [Re[adj_source_int[{d u}]]] ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name sadj_int_im  ; Value { Integral { [Im[adj_source_int[{d u}]]]; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name sadj_int_re1  ; Value { Integral { [Re[adj_source_int1[]]] ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name sadj_int_im1  ; Value { Integral { [Im[adj_source_int1[]]]; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } } */



              { Name sadj_int_re   ; Value { Local { [Re[adj_source_int[{d u}]]] ; In Omega; Jacobian JVol; } } }
              { Name sadj_int_im   ; Value { Local { [Im[adj_source_int[{d u}]]] ; In Omega; Jacobian JVol; } } }
              { Name sadj_int_re1   ; Value { Local { [Re[adj_source_int1[]]] ; In Omega; Jacobian JVol; } } }
              { Name sadj_int_im1   ; Value { Local { [Im[adj_source_int1[]]] ; In Omega; Jacobian JVol; } } }

              /* { Name sadj_re1_times_gradu   ; Value { Local { [Re[adj_source_int[{d u}]] * {d u}] ; In Omega; Jacobian JVol; } } } */


         }
    }

}

// #############################################################################

PostOperation {
{ Name postop; NameOfPostProcessing postpro ;
Operation {
	If (y_flag)
	 	If (save_solution)
			Print [ solution , OnElementsOf Omega, File "solutiony.pos", Name "uy"];
			Print [ vx , OnElementsOf Omega, File "vyx.pos", Name "vyx"];
			Print [ vy , OnElementsOf Omega, File "vyy.pos" , Name "vyy"];
		EndIf
		Print [ Ix[Omega], OnElementsOf PrintPoint, Format SimpleTable, File "Phiyx.txt"];
		Print [ Iy[Omega], OnElementsOf PrintPoint, Format SimpleTable, File "Phiyy.txt"];
	Else
	If (save_solution)
		Print [ solution , OnElementsOf Omega, File "solutionx.pos" , Name "ux"];
		Print [ vx , OnElementsOf Omega, File "vxx.pos" , Name "vxx"];
		Print [ vy , OnElementsOf Omega, File "vxy.pos" , Name "vxy"];
	EndIf
		Print [ Ix[Omega], OnElementsOf PrintPoint, Format SimpleTable, File "Phixx.txt"];
		Print [ Iy[Omega], OnElementsOf PrintPoint, Format SimpleTable, File "Phixy.txt"];

	EndIf
  Print [ I_inveps_xx[Omega], OnElementsOf PrintPoint, Format SimpleTable, File "I_inveps_xx.txt"];
  Print [ I_inveps_yy[Omega], OnElementsOf PrintPoint, Format SimpleTable, File "I_inveps_yy.txt"];
  Print [ I_inveps_yx[Omega], OnElementsOf PrintPoint, Format SimpleTable, File "I_inveps_yx.txt"];
  Print [ I_inveps_xy[Omega], OnElementsOf PrintPoint, Format SimpleTable, File "I_inveps_xy.txt"];
  Print [ V[Omega], OnElementsOf PrintPoint, Format SimpleTable, File "Vol.txt"];
	}
}
{ Name postop_source_adj; NameOfPostProcessing postpro ;
    Operation {
      /* Print[sadj_re1, OnElementsOf Omega, StoreInField  istore];
      Print[sadj_im1, OnElementsOf Omega, StoreInField  istore+1]; */
      /* Print[u_re, OnElementsOf Omega, StoreInField  istore];
      Print[u_im, OnElementsOf Omega, StoreInField  istore+1]; */
      Print[sadj_int_re, OnElementsOf Omega, StoreInField  istore];
      Print[sadj_int_im, OnElementsOf Omega, StoreInField  istore+1];
      /* Print[sadj_int_re1, OnElementsOf Omega, StoreInField  istore+2]; */
      /* Print[sadj_int_im1, OnElementsOf Omega, StoreInField  istore+3]; */
      /* Print[sadj_int_re, OnElementsOf Omega, File "sadj_int_re.pos", Name "sadj_int_re"]; */
      /* Print[sadj_int_im, OnElementsOf Omega, File "sadj_int_im.pos", Name "sadj_int_im"]; */

      /* Print[sadj_int_re1, OnElementsOf Omega, File "sadj_int_re1.pos", Name "sadj_int_re1"]; */
      /* Print[sadj_int_im1, OnElementsOf Omega, File "sadj_int_im1.pos", Name "sadj_int_im1"]; */


  }
}
{ Name postop_adjoint; NameOfPostProcessing postpro ;
    Operation {
    If (nodes_flag)
      Print[u_adj, OnElementsOf Omega ,  Format NodeTable, File "adjoint.txt" ];
    Else
      Print[u_adj, OnElementsOf Omega , Depth 0, Format SimpleTable, File "adjoint.txt" ];
    EndIf
    /* Print [ u_adj , OnElementsOf Omega, File "u_adj.pos", Name "adjoint"]; */

    /* Print [ sadj_var , OnElementsOf Omega, File "sadj_var.pos", Name "sadj_var"]; */

    /* Print [ sadj_scalfield , OnElementsOf Omega, File "sadj_scalfield.pos", Name "sadj_scalfield"]; */
  }
}
{ Name postop_solution; NameOfPostProcessing postpro ;
    Operation {
    If (nodes_flag)
      Print[solution, OnElementsOf Omega ,  Format NodeTable, File "u.txt" ];
    Else
      Print[solution, OnElementsOf Omega , Depth 0, Format SimpleTable, File "u.txt" ];
    EndIf
    /* Print [ solution , OnElementsOf Omega, File "u.pos", Name "solution"]; */

    /* Print [ sadj_var , OnElementsOf Omega, File "sadj_var.pos", Name "sadj_var"]; */

    /* Print [ sadj_scalfield , OnElementsOf Omega, File "sadj_scalfield.pos", Name "sadj_scalfield"]; */
  }
}

    { Name postop_int_objective; NameOfPostProcessing postpro ;
       Operation {
         /*Print [ u  , OnElementsOf Omega, File "u.pos" ];*/
       Print[ int_objective[Omega],  OnElementsOf PrintPoint, File "objective.txt" , Format SimpleTable ];
       }
    }

    { Name postop_dEq_deps; NameOfPostProcessing postpro ;
        Operation {
          If (nodes_flag)
            Print[dEq_deps_x, OnElementsOf Omega ,  Format NodeTable, File "dEq_deps_x.txt" ];
            Print[dEq_deps_y, OnElementsOf Omega ,  Format NodeTable, File "dEq_deps_y.txt" ];
          Else
            Print[dEq_deps_x, OnElementsOf Omega ,   Depth 0, Format SimpleTable, File "dEq_deps_x.txt" ];
            Print[dEq_deps_y, OnElementsOf Omega ,   Depth 0, Format SimpleTable, File "dEq_deps_y.txt" ];
          EndIf
      }
    }

    /* { Name postop_dEq_deps; NameOfPostProcessing postpro ;
        Operation {
          If (nodes_flag)
            Print[dEq_deps, OnElementsOf Omega ,  Format NodeTable, File "dEq_deps.txt" ];
          Else
            Print[dEq_deps, OnElementsOf Omega ,   Depth 0, Format SimpleTable, File "dEq_deps.txt" ];
          EndIf
      }
    } */


}

// #############################################################################

PostOperation {
    { Name postop_fields_pos; NameOfPostProcessing postpro ;
        Operation {

          Print [ epsilonr_xx , OnElementsOf Omega, File "epsilonr_xx.pos", Name "epsilonr_xx"];
          Print [ epsilonr_yy , OnElementsOf Omega, File "epsilonr_yy.pos", Name "epsilonr_yy"];
          If (y_flag)
            Print [ solution , OnElementsOf Omega, File "uy.pos", Name "uy"];
            Print [ vx , OnElementsOf Omega, File "vyx.pos", Name "vyx"];
            Print [ vy , OnElementsOf Omega, File "vyy.pos", Name "vyy"];
          Else
            Print [ solution , OnElementsOf Omega, File "ux.pos", Name "ux"];
            Print [ vx , OnElementsOf Omega, File "vxx.pos", Name "vxx"];
            Print [ vy , OnElementsOf Omega, File "vxy.pos", Name "vxy"];
          EndIf
        	}
    }
  { Name postop_fields_txt; NameOfPostProcessing postpro ;
        Operation {
              Print [ epsilonr_xx , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
              { Niy-1, Nix-1} ,Format SimpleTable, File "epsilonr_xx.txt" ];
              Print [ epsilonr_yy , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
              { Niy-1, Nix-1} ,Format SimpleTable, File "epsilonr_yy.txt" ];
              Print [ solution , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
              { Niy-1, Nix-1} ,Format SimpleTable, File "v.txt" ];
              Print [ vx , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
              { Niy-1, Nix-1} ,Format SimpleTable, File "vx.txt" ];
              Print [ vy , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
              { Niy-1, Nix-1} ,Format SimpleTable, File "vy.txt" ];


        	}
    }

    { Name postop_fields_circle; NameOfPostProcessing postpro ;
          Operation {
            For i In {0:N_pts_circ-1}
              t=i * 2*Pi/(N_pts_circ-1);
              xi = R_circ*Cos[t];
              yi = R_circ*Sin[t];
              Print [ solution , OnPoint {xi,yi,0} ,Format SimpleTable, File >> "sol_circ.txt" ];
              Print [ vx , OnPoint  {xi,yi,0}  ,Format SimpleTable, File >> "vx_circ.txt" ];
              Print [ vy , OnPoint  {xi,yi,0}  ,Format SimpleTable, File >> "vy_circ.txt" ];

            EndFor


            }
      }

}
