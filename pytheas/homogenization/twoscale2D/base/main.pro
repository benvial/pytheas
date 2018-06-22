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
		Galerkin { [1/epsilonr[]*Dof{d u} , {d u}];
                 		In Omega; Jacobian JVol; Integration Int_1; }
                 Galerkin { [ 1/epsilonr[] * E[] , {d u} ];
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
              { Name solution; Value { Local { [ {u}] ; In Omega; Jacobian JVol; } } }
							{ Name vx; Value { Local {[CompX[{d u}]] ; In Omega; Jacobian JVol; } } }
							{ Name vy; Value { Local {[CompY[{d u}]] ; In Omega; Jacobian JVol; } } }
              { Name Ix; Value { Integral { [CompX[1/epsilonr[]*{d u}]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name Iy; Value { Integral { [CompY[1/epsilonr[]*{d u}]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
							{ Name I_inveps; Value { Integral { [CompXX[1/epsilonr[]]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
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
		Print [ I_inveps[Omega], OnElementsOf PrintPoint, Format SimpleTable, File "I_inveps.txt"];
		Print [ V[Omega], OnElementsOf PrintPoint, Format SimpleTable, File "Vol.txt"];
	EndIf
	}
}


}
