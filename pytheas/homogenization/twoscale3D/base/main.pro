Include "parameters.dat";
Group {
  // Domains
	host         = Region[1002];
	If (inclusion_flag)
		incl         = Region[1001];
		Omega          = Region[{host,incl}];
	Else
		Omega          = Region[{host}];
	EndIf

	// Boundaries
	SurfBlochXm   = Region[505];
	SurfBlochXp  = Region[506];
	SurfBlochYm   = Region[501];
	SurfBlochYp  = Region[502];
	SurfBlochZm   = Region[504];
	SurfBlochZp  = Region[503];

	PrintPoint      = Region[10000];
}

Function{

	If (inclusion_flag)
	epsilonr[incl]            = Complex[eps_incl_re,eps_incl_im] * TensorDiag[1,1,1];
		If (coupling_flag)
	 	eps_xx[host] = Complex[ScalarField[XYZ[], 0, 1]{0}  ,ScalarField[XYZ[], 0, 1 ]{1} ];
		eps_yy[host] = Complex[ScalarField[XYZ[], 0, 1]{2}  ,ScalarField[XYZ[], 0, 1 ]{3} ];
		eps_zz[host] = Complex[ScalarField[XYZ[], 0, 1]{4}  ,ScalarField[XYZ[], 0, 1 ]{5} ];
		epsilonr[host]    = TensorDiag[eps_xx[],eps_yy[],eps_zz[]];
		Else
			epsilonr[host]           = Complex[eps_host_re,eps_host_im] * TensorDiag[1,1,1];
		EndIf
	Else
		epsilonr[Omega]    = Complex[ScalarField[XYZ[], 0, 1]{0}  ,ScalarField[XYZ[], 0, 1 ]{1} ]  * TensorDiag[1,1,1];
	EndIf

	EX[] = Vector[1,0,0] ;
	EY[] = Vector[0,1,0] ;
	EZ[] = Vector[0,0,1] ;
        }



Constraint {

		{Name Bloch;
		    Case {
              { Region SurfBlochXp; Type LinkCplx ; RegionRef SurfBlochXm;
                Coefficient Complex[1.0,0.0]; Function Vector[$X-dx,$Y,$Z] ;
              }
              { Region SurfBlochYp; Type LinkCplx ; RegionRef SurfBlochYm;
                Coefficient Complex[1.0,0.0]; Function Vector[$X,$Y-dy,$Z] ;
              }
              { Region SurfBlochZp; Type LinkCplx ; RegionRef SurfBlochZm;
                Coefficient Complex[1.0,0.0]; Function Vector[$X,$Y,$Z-dz] ;
              }
			 }
		}
}


Jacobian {
			{ Name JVol ; Case { { Region All ; Jacobian Vol ; } } }
			{ Name JSur ; Case { { Region All ; Jacobian Sur ; } }   }
			{ Name JLin ; Case { { Region All ; Jacobian Lin ; } } }
}

Integration {
  { Name Int_1 ;
    Case {
      { Type Gauss ;
        Case {
          { GeoElement Point       ; NumberOfPoints   4 ; }
          { GeoElement Line        ; NumberOfPoints  32 ; }
          { GeoElement Triangle    ; NumberOfPoints  16 ; } //1, 3, 4, 6, 7, 12, 13, 16
          { GeoElement Tetrahedron ; NumberOfPoints  29 ; }
          { GeoElement Prism       ; NumberOfPoints  51 ; }
	}
      }
    }
  }
}

FunctionSpace {
  { Name Hgrad; Type Form0;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_Node;    Support Region[Omega]; Entity NodesOf[Omega]; }
//    { Name sn2; NameOfCoef un2; Function BF_Node_2E; Support Region[Omega]; Entity EdgesOf[Omega]; }
     }
    Constraint {
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Bloch; }
//    { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
//    { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Bloch; }
   }
  }
}

Formulation {
	    {Name annex_x; Type FemEquation;
    		Quantity {{ Name ux; Type Local; NameOfSpace Hgrad;}}
					Equation {
						Galerkin { [epsilonr[]*Dof{d ux} , {d ux}];
						In Omega; Jacobian JVol; Integration Int_1; }
						Galerkin { [ -epsilonr[] * EX[] , {d ux} ];
						In Omega; Jacobian JVol; Integration Int_1; }
            }
      }
      {Name annex_y; Type FemEquation;
			Quantity {{ Name uy; Type Local; NameOfSpace Hgrad;}}
					Equation {
							Galerkin { [epsilonr[]*Dof{d uy} , {d uy}];
							In Omega; Jacobian JVol; Integration Int_1; }
							Galerkin { [ -epsilonr[] * EY[] , {d uy} ];
							In Omega; Jacobian JVol; Integration Int_1; }
					}
			}
				{Name annex_z; Type FemEquation;
					Quantity {{ Name uz; Type Local; NameOfSpace Hgrad;}}
						Equation {
							Galerkin { [epsilonr[]*Dof{d uz} , {d uz}];
							In Omega; Jacobian JVol; Integration Int_1; }
							Galerkin { [ -epsilonr[] * EZ[] , {d uz} ];
							In Omega; Jacobian JVol; Integration Int_1; }
						}
					}

        }



Resolution {
  { Name annex_x;
    System {
      { Name Sx; NameOfFormulation annex_x; Type ComplexValue;}
    }
    Operation {
      Generate[Sx] ;
      // Print[S];
      Solve[Sx] ;
      SaveSolution[Sx] ;
    }
  }
    { Name annex_y;
    System {
      { Name Sy; NameOfFormulation annex_y; Type ComplexValue;}
    }
    Operation {
      Generate[Sy] ;
      // Print[S];
      Solve[Sy] ;
      SaveSolution[Sy] ;
    }
  }
      { Name annex_z;
    System {
      { Name Sz; NameOfFormulation annex_z; Type ComplexValue;}
    }
    Operation {
      Generate[Sz] ;
      // Print[S];
      Solve[Sz] ;
      SaveSolution[Sz] ;
    }
  }
}

PostProcessing {
    { Name postpro_x; NameOfFormulation annex_x; NameOfSystem Sx;
            Quantity {
               /* { Name solutionx; Value { Local { [Re[ {ux}] ]  ; In Omega; Jacobian JVol; } } } */
              { Name Phixx; Value { Integral { [CompX[epsilonr[]*{d ux}]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name Phixy; Value { Integral { [CompY[epsilonr[]*{d ux}]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name Phixz; Value { Integral { [CompZ[epsilonr[]*{d ux}]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name Int_eps_xx; Value { Integral { [CompXX[epsilonr[]]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              /* { Name V; Value { Integral { [1]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } } */
		     }
    }
    { Name postpro_y; NameOfFormulation annex_y; NameOfSystem Sy;
            Quantity {
               /* { Name solutiony; Value { Local { [Re[ {uy}] ] ; In Omega; Jacobian JVol; } } } */
              { Name Phiyx; Value { Integral { [CompX[epsilonr[]*{d uy}]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name Phiyy; Value { Integral { [CompY[epsilonr[]*{d uy}]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name Phiyz; Value { Integral { [CompZ[epsilonr[]*{d uy}]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
							{ Name Int_eps_yy; Value { Integral { [CompYY[epsilonr[]]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }

		     }
    }
    { Name postpro_z; NameOfFormulation annex_z; NameOfSystem Sz;
            Quantity {
               /* { Name solutionz; Value { Local { [Re[ {uz}] ] ; In Omega; Jacobian JVol; } } } */
              { Name Phizx; Value { Integral { [CompX[epsilonr[]*{d uz}]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name Phizy; Value { Integral { [CompY[epsilonr[]*{d uz}]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name Phizz; Value { Integral { [CompZ[epsilonr[]*{d uz}]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
							{ Name Int_eps_zz; Value { Integral { [CompZZ[epsilonr[]]]  ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }

		     }
    }

}

PostOperation {
{ Name postop_x; NameOfPostProcessing postpro_x ;
Operation {
	/* Print [ solutionx , OnElementsOf Omega, File "solutionx.pos" ]; */
	Print [ Phixx[Omega], OnElementsOf PrintPoint, Format TimeTable, File "Phixx.txt"];
	Print [ Phixy[Omega], OnElementsOf PrintPoint, Format TimeTable, File "Phixy.txt"];
	Print [ Phixz[Omega], OnElementsOf PrintPoint, Format TimeTable, File "Phixz.txt"];
	Print [ Int_eps_xx[Omega], OnElementsOf PrintPoint, Format TimeTable, File "Int_eps_xx.txt"];
	/*Print [ V[Omega], OnElementsOf PrintPoint, Format TimeTable, File "V.txt"];*/
	}
}

{ Name postop_y; NameOfPostProcessing postpro_y ;
Operation {
 	/* Print [ solutiony , OnElementsOf Omega, File "solutiony.pos" ]; */
	Print [ Phiyx[Omega], OnElementsOf PrintPoint, Format TimeTable, File "Phiyx.txt"];
	Print [ Phiyy[Omega], OnElementsOf PrintPoint, Format TimeTable, File "Phiyy.txt"];
	Print [ Phiyz[Omega], OnElementsOf PrintPoint, Format TimeTable, File "Phiyz.txt"];
	Print [ Int_eps_yy[Omega], OnElementsOf PrintPoint, Format TimeTable, File "Int_eps_yy.txt"];

}
}
{ Name postop_z; NameOfPostProcessing postpro_z ;
Operation {
 	/* Print [ solutionz , OnElementsOf Omega, File "solutionz.pos" ]; */
	Print [ Phizz[Omega], OnElementsOf PrintPoint, Format TimeTable, File "Phizz.txt"];
	Print [ Phizy[Omega], OnElementsOf PrintPoint, Format TimeTable, File "Phizy.txt"];
	Print [ Phizx[Omega], OnElementsOf PrintPoint, Format TimeTable, File "Phizx.txt"];
	Print [ Int_eps_zz[Omega], OnElementsOf PrintPoint, Format TimeTable, File "Int_eps_zz.txt"];

}
}

}
