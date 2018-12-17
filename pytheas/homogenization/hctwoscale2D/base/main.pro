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
  	SurfInclusion  = Region[100];

  	// Points
  	PrintPoint      = Region[10000];
}
// #############################################################################

Function{

  ksearch = 2.*Pi/lambda0search;

  If (inclusion_flag)
    
    epsilonr[incl]           = Complex[eps_incl_re,eps_incl_im] * TensorDiag[1,1,1];
    epsilonr[host] = Complex[eps_host_re,eps_host_im] * TensorDiag[1,1,1];
    /* epsilonr_xx[]  = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}]; */
    /* epsilonr_yy[]  = Complex[ScalarField[XYZ[], 0, 1 ]{2}, ScalarField[XYZ[], 0, 1 ]{3}]; */
    /* epsilonr_zz[]  = Complex[ScalarField[XYZ[], 0, 1 ]{4}, ScalarField[XYZ[], 0, 1 ]{5}]; */
    /* epsilonr[incl] =  TensorDiag[epsilonr_xx[],epsilonr_yy[],epsilonr_zz[]]; */

  Else
    /* epsilonr[design]         = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}] * TensorDiag[1,1,1]; */
    epsilonr_xx[Omega]  = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}];
    epsilonr_yy[Omega]  = Complex[ScalarField[XYZ[], 0, 1 ]{2}, ScalarField[XYZ[], 0, 1 ]{3}];
    epsilonr_zz[Omega]  = Complex[ScalarField[XYZ[], 0, 1 ]{4}, ScalarField[XYZ[], 0, 1 ]{5}];

    epsilonr[Omega] =  TensorDiag[epsilonr_xx[],epsilonr_yy[],epsilonr_zz[]];
  EndIf

  xsi[] = TensorDiag[1/CompYY[epsilonr[]],1/CompXX[epsilonr[]],0];

		If (y_flag)
			E[] = Vector[0,1,0] ;
		Else
			E[] = Vector[1,0,0] ;
		EndIf
}

// #############################################################################

Constraint {


  {Name Dirichlet; Type Assign;
      Case {
          { Region SurfInclusion; Value 0.; }
      }
  }

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

FunctionSpace {
  { Name Hgrad_incl; Type Form0;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_Node;    Support Region[incl]; Entity NodesOf[incl]; }
      { Name sn2; NameOfCoef un2; Function BF_Node_2E; Support Region[incl]; Entity EdgesOf[incl]; }
     }
    Constraint {
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Dirichlet; }
      /* { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Bloch; } */
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
      /* { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Bloch; } */
   }
  }

}

// #############################################################################

Formulation {
	    {Name electrostat_scalar_host; Type FemEquation;
    		Quantity {{ Name u; Type Local; NameOfSpace Hgrad;}}
		/* Equation {
        		Galerkin { [Dof{d u} , {d u}];
                     		In host; Jacobian JVol; Integration Int_1; }
            Galerkin { [ Dof{d u} * Normal[] , {u} ];
                  		In SurfInclusion; Jacobian JSur; Integration Int_1; }
            Galerkin { [  -Normal[] * E[]  , {u} ];
                       In SurfInclusion; Jacobian JSur; Integration Int_1; }
                } */

    Equation {
        		Galerkin { [xsi[] * Dof{d u} , {d u}];
                     		In host; Jacobian JVol; Integration Int_1; }
            Galerkin { [xsi[] * E[] , {d u} ];
                  		In host; Jacobian JVol; Integration Int_1; }
                }
            }

        }
// #############################################################################


Formulation {
	    {Name spectral_problem_incl; Type FemEquation;
    		Quantity {{ Name phi; Type Local; NameOfSpace Hgrad_incl;}}

    Equation {
      Galerkin {  DtDtDof[ -epsilonr[]*Dof{phi} , {phi}];
      In incl; Jacobian JVol; Integration Int_1;  }
      Galerkin { [-Dof{d phi} , {d phi}];
      In incl; Jacobian JVol; Integration Int_1; }
                }
            }

        }
// #############################################################################



Resolution {
  { Name electrostat_scalar_host;
    System {
      { Name S; NameOfFormulation electrostat_scalar_host; Type ComplexValue;}
    }
    Operation {
      Generate[S] ;  Solve[S] ;SaveSolution[S] ;
    }
  }

  { Name spectral_problem_incl;
      System {
      { Name Smodal; NameOfFormulation spectral_problem_incl; Type ComplexValue;}
      }
      Operation {
      GenerateSeparate[Smodal]; EigenSolve[Smodal,neig,ksearch^2,0]; SaveSolutions[Smodal];
      }

  }

}

// #############################################################################


PostProcessing {
    { Name postpro; NameOfFormulation electrostat_scalar_host; NameOfSystem S;
            Quantity {
              { Name epsilonr_xx; Value { Local { [CompXX[ epsilonr[]]] ; In Omega; Jacobian JVol; } } }
              { Name epsilonr_yy; Value { Local { [CompYY[ epsilonr[]]] ; In Omega; Jacobian JVol; } } }
              { Name solution; Value { Local { [ {u}] ; In Omega; Jacobian JVol; } } }
							{ Name vx; Value { Local {[CompX[{d u}]] ; In Omega; Jacobian JVol; } } }
							{ Name vy; Value { Local {[CompY[{d u}]] ; In Omega; Jacobian JVol; } } }
              { Name Ix; Value { Integral { [CompX[1/epsilonr[]*{d u}]]  ; In host    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name Iy; Value { Integral { [CompY[1/epsilonr[]*{d u}]]  ; In host    ; Integration Int_1 ; Jacobian JVol ; } } }
							{ Name I_inveps_xx; Value { Integral { [CompXX[1/epsilonr[]]]  ; In host    ; Integration Int_1 ; Jacobian JVol ; } } }
          		{ Name I_inveps_yy; Value { Integral { [CompYY[1/epsilonr[]]]  ; In host    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name V; Value { Integral { [1]  ; In host    ; Integration Int_1 ; Jacobian JVol ; } } }

         }
    }

    { Name postpro_modal; NameOfFormulation spectral_problem_incl;
      Quantity {
        { Name EigenValues;  Value { Local{ [$EigenvalueImag]; In PrintPoint; Jacobian JVol; } } }
        { Name EigenVectors ;   Value { Local { [     {phi}  ] ; In incl; Jacobian JVol; } } }
        { Name Normalization ; Value { Integral { [ {phi}*Conj[{phi}] ] ; In incl    ; Integration Int_1 ; Jacobian JVol ; } } }
        { Name Coefs_mu ; Value { Integral { [ {phi} ] ; In incl    ; Integration Int_1 ; Jacobian JVol ; } } }
        { Name V_incl ; Value { Integral { [ 1 ] ; In incl    ; Integration Int_1 ; Jacobian JVol ; } } }

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
		Print [ I_inveps_xx[Omega], OnElementsOf PrintPoint, Format SimpleTable, File "I_inveps_xx.txt"];
		Print [ I_inveps_yy[Omega], OnElementsOf PrintPoint, Format SimpleTable, File "I_inveps_yy.txt"];
		Print [ V[Omega], OnElementsOf PrintPoint, Format SimpleTable,LastTimeStepOnly, File "Vol.txt"];
	EndIf
	}
}

{ Name postop_eigenvalues; NameOfPostProcessing postpro_modal ;
    Operation {
      Print [EigenValues, OnElementsOf PrintPoint, Format TimeTable, File "EigenValues.txt"];
      }
}
{ Name postop_eigenvectors_pos; NameOfPostProcessing postpro_modal ;
    Operation {
        Print [ EigenVectors   , OnElementsOf Omega, File "EigenVectors.pos" , EigenvalueLegend];
    }
}

{ Name postop_eigenvectors_txt; NameOfPostProcessing postpro_modal ;
    Operation {
    Print [ EigenVectors , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
    { Niy-1, Nix-1} , Format TimeTable, File "EigenVectors.txt" ];
    }
}
{ Name postop_norm_eigenvectors; NameOfPostProcessing postpro_modal ;
    Operation {
      Print [Normalization[Omega], OnElementsOf PrintPoint, Format TimeTable, File "NormsEigenVectors.txt"];
      }
}

{ Name postop_coefs_mu; NameOfPostProcessing postpro_modal ;
    Operation {
      Print [Coefs_mu[Omega], OnElementsOf PrintPoint, Format TimeTable, File "Coefs_mu.txt"];
      }
}

{ Name postop_V_incl; NameOfPostProcessing postpro_modal ;
    Operation {
      Print [V_incl[Omega], OnElementsOf PrintPoint, Format SimpleTable, LastTimeStepOnly, File "V_incl.txt"];
      }
}


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
    /* { Name postop_fields_txt; NameOfPostProcessing postpro ;
        Operation {
              Print [ v , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
              { Niy-1, Nix-1} ,Format SimpleTable, File "v.txt" ];
              Print [ vx , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
              { Niy-1, Nix-1} ,Format SimpleTable, File "vx.txt" ];
              Print [ vy , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
              { Niy-1, Nix-1} ,Format SimpleTable, File "vy.txt" ];


        	}
    } */

}
