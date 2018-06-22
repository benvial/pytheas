Include "parameters.dat";


// #############################################################################
// #############################################################################

Group {
    // Domains
    cell         = Region[1000];
    Omega          = Region[{cell}];

    // Boundaries
  	SurfBlochLeft   = Region[101];
  	SurfBlochRight  = Region[103];
  	SurfBlochTop   = Region[102];
  	SurfBlochBot  = Region[104];

  	// Points
  	PrintPoint      = Region[10000];
}
// #############################################################################

Function{
    epsilonr[Omega]  = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}] * TensorDiag[1,1,1];
    contrast_epsilonr[Omega] = epsilonr[] - TensorDiag[1,1,1];
    ex[] = Vector[1,0,0] ;
    E_bias[] = E_static * ex[];
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

// #############################################################################

FunctionSpace {
  { Name Hgrad; Type Form0;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_Node;    Support Region[Omega]; Entity NodesOf[Omega]; }
      { Name sn2; NameOfCoef un2; Function BF_Node_2E; Support Region[Omega]; Entity EdgesOf[Omega]; }
     }
    Constraint {
      /*{ NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Dirichlet; }*/
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Bloch; }
      /*{ NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }*/
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Bloch; }
   }
  }

}

// #############################################################################

Formulation {
	    {Name electrostat_scalar; Type FemEquation;
    		Quantity {{ Name v; Type Local; NameOfSpace Hgrad;}}
		Equation {
          		Galerkin { [epsilonr[]*Dof{d v} , {d v}];
              In Omega; Jacobian JVol; Integration Int_1; }
              Galerkin { [ contrast_epsilonr[] * E_bias[] , {d v} ];
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
      Generate[S] ;  Solve[S] ;SaveSolution[S] ;
    }
  }

}

// #############################################################################


PostProcessing {
    { Name postpro; NameOfFormulation electrostat_scalar; NameOfSystem S;
            Quantity {
              { Name v; Value { Local { [{v}] ; In Omega; Jacobian JVol; } } }
              { Name vx; Value { Local {[-CompX[{d v}]] ; In Omega; Jacobian JVol; } } }
              { Name vy; Value { Local {[-CompY[{d v}]] ; In Omega; Jacobian JVol; } } }
              { Name abs_dv; Value { Local {[Sqrt[(CompX[{d v}])^2 + (CompY[{d v}])^2]] ; In Omega; Jacobian JVol; } } }

		     }
    }

}

// #############################################################################

PostOperation {
    { Name postop_fields_pos; NameOfPostProcessing postpro ;
        Operation {
        			Print [ v , OnElementsOf Omega, File "solution.pos", Name "v"];
              Print [ vx , OnElementsOf Omega, File "vx.pos", Name "vx"];
              Print [ vy , OnElementsOf Omega, File "vy.pos", Name "vy"];
              /*Print [ abs_dv , OnElementsOf Omega, File "abs_dv.pos", Name "abs_dv"];*/

        	}
    }
    { Name postop_fields_txt; NameOfPostProcessing postpro ;
        Operation {
              Print [ v , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
              { Niy-1, Nix-1} ,Format SimpleTable, File "v.txt" ];
              Print [ vx , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
              { Niy-1, Nix-1} ,Format SimpleTable, File "vx.txt" ];
              Print [ vy , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
              { Niy-1, Nix-1} ,Format SimpleTable, File "vy.txt" ];


        	}
    }

}
