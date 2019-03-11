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

  	// Points
  	PrintPoint      = Region[10000];
}
// #############################################################################

Function{
  If (inclusion_flag)
    /* epsilonr[host]         = Complex[eps_host_re,eps_host_im] * TensorDiag[1,1,1]; */
    epsilonr[incl]           = Complex[eps_incl_re,eps_incl_im] * TensorDiag[1,1,1];
    If (aniso)
        epsilonr_xx[]  = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}];
        epsilonr_yy[]  = Complex[ScalarField[XYZ[], 0, 1 ]{2}, ScalarField[XYZ[], 0, 1 ]{3}];
        epsilonr_zz[]  = Complex[ScalarField[XYZ[], 0, 1 ]{4}, ScalarField[XYZ[], 0, 1 ]{5}];

        epsilonr[host] =  TensorDiag[epsilonr_xx[],epsilonr_yy[],epsilonr_zz[]];
      Else
        epsilonr[host]         = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}] * TensorDiag[1,1,1];
      EndIf
  Else
  If (aniso)
      epsilonr_xx[]  = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}];
      epsilonr_yy[]  = Complex[ScalarField[XYZ[], 0, 1 ]{2}, ScalarField[XYZ[], 0, 1 ]{3}];
      epsilonr_zz[]  = Complex[ScalarField[XYZ[], 0, 1 ]{4}, ScalarField[XYZ[], 0, 1 ]{5}];

      epsilonr[Omega] =  TensorDiag[epsilonr_xx[],epsilonr_yy[],epsilonr_zz[]];
    Else
      epsilonr[Omega]         = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}] * TensorDiag[1,1,1];
    EndIf
  EndIf

    xsi[Omega] = TensorDiag[1/CompYY[epsilonr[]], 1/CompXX[epsilonr[]],1];


		If (y_flag)
			E[] = Vector[0,1,0] ;
		Else
			E[] = Vector[1,0,0] ;
		EndIf
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
    		Quantity {{ Name u; Type Local; NameOfSpace Hgrad;}}
		Equation {
		Galerkin { [xsi[]*Dof{d u} , {d u}];
                 		In Omega; Jacobian JVol; Integration Int_1; }
                 Galerkin { [ xsi[] * E[] , {d u} ];
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
              { Name epsilonr_xx; Value { Local { [CompXX[ epsilonr[]]] ; In Omega; Jacobian JVol; } } }
              { Name epsilonr_yy; Value { Local { [CompYY[ epsilonr[]]] ; In Omega; Jacobian JVol; } } }
              { Name solution; Value { Local { [ {u}] ; In Omega; Jacobian JVol; } } }
							{ Name vx; Value { Local {[CompX[{d u}]] ; In Omega; Jacobian JVol; } } }
							{ Name vy; Value { Local {[CompY[{d u}]] ; In Omega; Jacobian JVol; } } }
              { Name Ix; Value { Integral { [CompX[xsi[]*{d u}]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name Iy; Value { Integral { [CompY[xsi[]*{d u}]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
							{ Name I_inveps_xx; Value { Integral { [CompXX[xsi[]]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
          		{ Name I_inveps_yy; Value { Integral { [CompYY[xsi[]]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name V; Value { Integral { [1]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }

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
		Print [ V[Omega], OnElementsOf PrintPoint, Format SimpleTable, File "Vol.txt"];
	EndIf
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

}
