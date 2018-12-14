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

    /* epsilonr_xx[Omega]  = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}]; */
    /* epsilonr_yy[Omega]  = Complex[ScalarField[XYZ[], 0, 1 ]{2}, ScalarField[XYZ[], 0, 1 ]{3}]; */
    /* epsilonr[Omega] =  TensorDiag[epsilonr_xx[],epsilonr_yy[],1]; */
    dephx[] = Complex[ Cos[kx*dx]       , Sin[kx*dx]      ];
  	dephy[] = Complex[ Cos[ky*dy]       , Sin[ky*dy]      ];
    shift = (2*Pi/lambda0search)^2;
}

// #############################################################################

Constraint {

		{Name Bloch;
		    Case {
                    { Region SurfBlochRight; Type LinkCplx ; RegionRef SurfBlochLeft;
                      Coefficient dephx[]; Function Vector[$X-dx,$Y,$Z] ;
                    }
                    { Region SurfBlochTop; Type LinkCplx ; RegionRef SurfBlochBot;
                      Coefficient dephy[]; Function Vector[$X,$Y-dy,$Z] ;
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
	    {Name wave_eq_scalar_modal_TM; Type FemEquation;
    		Quantity {{ Name h; Type Local; NameOfSpace Hgrad;}}
		Equation {
		Galerkin { [-1/epsilonr[] *Dof{d h} , {d h}];
                 		In Omega; Jacobian JVol; Integration Int_1; }
                 Galerkin { DtDtDof[ - Dof{h} , {h} ];
                  		In Omega; Jacobian JVol; Integration Int_1; }
                }
            }
	    {Name wave_eq_scalar_modal_TE; Type FemEquation;
    		Quantity {{ Name e; Type Local; NameOfSpace Hgrad;}}
		Equation {
		Galerkin { [-1 *Dof{d e} , {d e}];
                 		In Omega; Jacobian JVol; Integration Int_1; }
                 Galerkin { DtDtDof[ -epsilonr[]* Dof{e} , {e} ];
                  		In Omega; Jacobian JVol; Integration Int_1; }
                }
            }

}
// #############################################################################

Resolution {
  { Name TE;
    System {
      { Name MTE; NameOfFormulation wave_eq_scalar_modal_TE; Type ComplexValue;}
    }
    Operation {
    GenerateSeparate[MTE]; EigenSolve[MTE,neig,shift,0]; SaveSolutions[MTE];
    }
  }
  { Name TM;
    System {
      { Name MTM; NameOfFormulation wave_eq_scalar_modal_TM; Type ComplexValue;}
    }
    Operation {
    GenerateSeparate[MTM]; EigenSolve[MTM,neig,shift,0]; SaveSolutions[MTM];
    }
  }


}

// #############################################################################




PostProcessing {
    	{ Name postpro_TE; NameOfFormulation wave_eq_scalar_modal_TE; NameOfSystem MTE;
		Quantity {
			{ Name EV_TE;  Value { Local{ [$EigenvalueImag]; In PrintPoint; Jacobian JVol; } } }
			{ Name modes_TE; Value { Local { [ {e}]  ; In Omega; Jacobian JVol; } } }
		}
	}
	{ Name postpro_TM; NameOfFormulation wave_eq_scalar_modal_TM; NameOfSystem MTM;
		Quantity {
			{ Name EV_TM;  Value { Local{ [$EigenvalueImag]; In PrintPoint; Jacobian JVol; } } }
			{ Name modes_TM; Value { Local { [ {h}]  ; In Omega; Jacobian JVol; } } }
		}
	}

}

// #############################################################################

PostOperation {


      /* ----------- Modal analysis --------------------*/
      { Name postop_eigenvalues_TE; NameOfPostProcessing postpro_TE ;
          Operation {
            Print [EV_TE, OnElementsOf PrintPoint, Format TimeTable, File "EV_TE.txt"];
            }
      }
      { Name postop_eigenvectors_TE_pos; NameOfPostProcessing postpro_TE ;
          Operation {
              Print [ modes_TE   , OnElementsOf Omega, File "modes_TE.pos" , EigenvalueLegend];
          }
      }

      { Name postop_eigenvectors_TE_txt; NameOfPostProcessing postpro_TE ;
          Operation {
          Print [ modes_TE , OnPlane    { { -dx/2,-dy/2,0 } { -dx/2,dy/2,0 } { dx/2,-dy/2,0 } }
          { Niy-1, Nix-1} , Format TimeTable, File "modes_TE.txt" ];
          }
      }

      { Name postop_eigenvalues_TM; NameOfPostProcessing postpro_TM ;
          Operation {
            Print [EV_TM, OnElementsOf PrintPoint, Format TimeTable, File "EV_TM.txt"];
            }
      }
      { Name postop_eigenvectors_TM_pos; NameOfPostProcessing postpro_TM ;
          Operation {
              Print [ modes_TM   , OnElementsOf Omega, File "modes_TM.pos" , EigenvalueLegend];
          }
      }

      { Name postop_eigenvectors_TM_txt; NameOfPostProcessing postpro_TM ;
          Operation {
          Print [ modes_TM , OnPlane    { { -dx/2,-dy/2,0 } { -dx/2,dy/2,0 } { dx/2,-dy/2,0 } }
          { Niy-1, Nix-1} , Format TimeTable, File "modes_TM.txt" ];
          }
      }


}
