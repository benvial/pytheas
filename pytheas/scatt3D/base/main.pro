Include "parameters.dat";
// #############################################################################
Group {
	// SubDomains
	host         = Region[1];
	des         = Region[2];
	sph         = Region[3];
	PMLcorners   = Region[4];
	PMLxy   = Region[5];
	PMLyz   = Region[6];
	PMLxz   = Region[7];
	PMLx   = Region[8];
	PMLy   = Region[9];
	PMLz   = Region[10];
	// Boundaries
	// Lines
	// Domains
	PML = Region[{PMLcorners,PMLxy,PMLxz, PMLyz, PMLx, PMLy, PMLz}];
	Omega_nosource = Region[{PML,host,sph}];
	Omega_source   = Region[des];
	Omega_no_pml   = Region[{host,des,sph}];
	Omega_design   = Region[des];
	Omega_target   = Region[sph];
	Omega          = Region[{PML, Omega_no_pml}];
	// Points
  PrintPoint	=  Region[11];

	// SurfNeumann    = Region[{SurfBlochXm,SurfBlochXp,SurfBlochYm,SurfBlochYp}];
}


// #############################################################################
Function{
	Freq             = cel/lambda0;
	omega0           = 2*Pi*Freq;
	k0               = 2.*Pi/lambda0;
	Ae               = 1;
	Ah               = Ae*Sqrt[epsilon0/mu0];
	alpha0           = k0*Sin[theta_0]*Cos[phi_0];
	beta0            = k0*Sin[theta_0]*Sin[phi_0];
	gamma0           = k0*Cos[theta_0];
	Pinc             = 0.5*Ae*Ae*Sqrt[epsilon0/mu0] * Cos[theta_0];


	epsilon[host]  = Complex[eps_host_re , eps_host_im] * TensorDiag[1,1,1];
	epsilon[sph]   = Complex[eps_sph_re , eps_sph_im] * TensorDiag[1,1,1];
	epsilon[Omega_design]    = Complex[ScalarField[XYZ[], 0, 1]{0}  ,ScalarField[XYZ[], 0, 1 ]{1} ]  * TensorDiag[1,1,1];
	/* epsilon[Omega_design]  = Complex[eps_des_re , eps_des_im] * TensorDiag[1,1,1]; */

	mu[Omega_no_pml]       = TensorDiag[1,1,1];
	

// PML parameters
	sx[]            = Complex[a_pml,-b_pml];
	sy[]            = Complex[a_pml,-b_pml];
	sz[]            = Complex[a_pml,-b_pml];
	pml_tens[]      = TensorDiag[$3*$2/$1,$1*$3/$2,$1*$2/$3];

	epsilon[PMLcorners] = eps_host_re*pml_tens[sx[],sy[],sz[]];
	mu[PMLcorners] = pml_tens[sx[],sy[],sz[]];
	epsilon[PMLx] = eps_host_re*pml_tens[sx[],1,1];
	mu[PMLx] = pml_tens[sx[],1,1];
	epsilon[PMLy] = eps_host_re*pml_tens[1,sy[],1];
	mu[PMLy] = pml_tens[1,sy[],1];
	epsilon[PMLz] = eps_host_re*pml_tens[1,1,sz[]];
	mu[PMLz] = pml_tens[1,1,sz[]];
	epsilon[PMLxy] = eps_host_re*pml_tens[sx[],sy[],1];
	mu[PMLxy] = pml_tens[sx[],sy[],1];
	epsilon[PMLxz] = eps_host_re*pml_tens[sx[],1,sz[]];
	mu[PMLxz] = pml_tens[sx[],1,sz[]];
	epsilon[PMLyz] = eps_host_re*pml_tens[1,sy[],sz[]];
	mu[PMLyz] = pml_tens[1,sy[],sz[]];
	
	epsilon_annex[Omega_nosource]=epsilon[];
	epsilon_annex[Omega_source]=Complex[eps_host_re , eps_host_im] * TensorDiag[1,1,1];
	
	Propp[]         = Complex[ Cos[alpha0*X[]+beta0*Y[]+gamma0*(Z[]-0.0)] , Sin[alpha0*X[]+beta0*Y[]+gamma0*(Z[]-0.0)] ];
	
	Ex0[] = Propp[]*Ae *(Cos[psi_0] * Cos[theta_0] * Cos[phi_0]	- Sin[psi_0] * Sin[phi_0] );
	Ey0[] = Propp[]*Ae * (	Cos[(psi_0)] * Cos[(theta_0)] * Sin[(phi_0)]	+ Sin[(psi_0)] * Cos[(phi_0)]);
	Ez0[] = Propp[]*Ae * (-Cos[(psi_0)] * Sin[(theta_0)]);
	E0[]            = Vector[Ex0[],Ey0[],Ez0[]];
	
	source[]         = k0^2*(epsilon[]-epsilon_annex[])*E0[];//(nm/1.e-9)^2*

	// Topology optimization
	  coef_obj[] =  1/(4/3*Pi*R_sph^3);
		objective[] = coef_obj[] * SquNorm[$1 + E0[] ];
		adj_source_int[] = -2 * coef_obj[] * Conj[$1 + E0[]];
		db_deps[] = -k0^2*E0[];
		dA_deps[] = k0^2;
		dEq_deps[] = db_deps[] - dA_deps[] * ($1);
		/* source_adj[] = Vector[0,0,0]; */

}
// #############################################################################
Constraint {
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
          { GeoElement Point       ; NumberOfPoints   4 ; }
          { GeoElement Line        ; NumberOfPoints  32 ; }
          { GeoElement Triangle    ; NumberOfPoints  16 ; } //1, 3, 4, 6, 7, 12, 13, 16
          { GeoElement Tetrahedron ; NumberOfPoints  29 ; }
					{ GeoElement Hexahedron  ; NumberOfPoints  34 ; }
          { GeoElement Prism       ; NumberOfPoints  51 ; }
        }
      }
    }
  }
}
// #############################################################################
FunctionSpace {
  { Name Hcurl; Type Form1;
    BasisFunction {
      { Name sn; NameOfCoef un; Function BF_Edge;
        Support Region[{Omega}]; Entity EdgesOf[All]; }
     /* { Name sn2; NameOfCoef un2; Function BF_Edge_2E;
        Support Region[{Omega}]; Entity EdgesOf[All]; } */
				If (el_order== 2)
							{ Name sn3; NameOfCoef un3; Function BF_Edge_3F_b;
				         Support Region[Omega]; Entity FacetsOf[Omega]; }
				      { Name sn4; NameOfCoef un4; Function BF_Edge_3F_c;
				         Support Region[Omega]; Entity FacetsOf[Omega]; }
				      { Name sn5; NameOfCoef un5; Function BF_Edge_4E;
				         Support Region[Omega]; Entity EdgesOf[Omega]; }
				 EndIf
    }
    Constraint {

      }
  }
}

// #############################################################################
Formulation {{Name helmholtz_vector; Type FemEquation;
    		Quantity {{ Name u; Type Local; NameOfSpace Hcurl;}}
		Equation { Galerkin {[-1/mu[]*Dof{Curl u} , {Curl u}];
                 		In Omega; Jacobian JVol; Integration Int_1;  }
                   	   Galerkin { [k0^2*epsilon[]*Dof{u} , {u}];
                 		In Omega; Jacobian JVol; Integration Int_1;  }
                   	   Galerkin {[ ($Source ? source[] : 0) , {u}];
                 		In Omega; Jacobian JVol; Integration Int_1;  }
											Galerkin {[ ($SourceAdj ? ComplexVectorField [XYZ[]]{2} : 0) , {u}];
									 In Omega_target; Jacobian JVol; Integration Int_1;  }
								}
            }
        }
// #############################################################################
Resolution {
  { Name helmholtz_vector;
    System {
      { Name Maxwell; NameOfFormulation helmholtz_vector; Type ComplexValue; Frequency Freq; }
    }
    Operation {
	Evaluate[$Source = 1, $SourceAdj = 0];
	Evaluate[$t0=GetWallClockTime[]];
	Generate[Maxwell]; Solve[Maxwell]; SaveSolution[Maxwell];
	Evaluate[$t1=GetWallClockTime[]-$t0];
	Print[{$t1}, Format "direct: %g s"];

	If (adjoint_flag)

	PostOperation[postop_int_objective];
	PostOperation[postop_dEq_deps];
	PostOperation[postop_source_adj];
	Evaluate[$Source = 0, $SourceAdj = 1];
	/* Evaluate[$source_adj = VectorField[XYZ[]]{2}];//Complex[VectorField[XYZ[]]{2}, VectorField[XYZ[]]{3}]]; */
	/* Evaluate[$source_adj = Vector[1,0,0]]; */

	/* PostOperation[postop_source_adj_test]; */

	Evaluate[$t0=GetWallClockTime[]];
	GenerateRHSGroup[Maxwell, Omega_target]; SolveAgain[Maxwell] ; SaveSolution[Maxwell] ;
	Evaluate[$t1=GetWallClockTime[]-$t0];
 	Print[{$t1}, Format "adjoint: %g s"];
	  PostOperation[postop_adjoint];
	EndIf

    }
  }
}


        // Hinc[] : Complex[0,1] * 1/omega0 * 1/mu0 * Curl Einc[];
        // H_d    : Complex[0,1] * 1/omega0 * 1/mu0 * {Curl u};
// #############################################################################
PostProcessing {
    { Name postpro_diff; NameOfFormulation helmholtz_vector;NameOfSystem Maxwell;
            Quantity {

              { Name Etot; Value { Local { [ {u}+E0[] ]; In Omega; Jacobian JVol; } } }
							{ Name Etotx; Value { Local { [ CompX[{u}+E0[]]  ]; In Omega; Jacobian JVol; } } }
							{ Name Etoty; Value { Local { [ CompY[{u}+E0[]]  ]; In Omega; Jacobian JVol; } } }
							{ Name Etotz; Value { Local { [ CompZ[{u}+E0[]]  ]; In Omega; Jacobian JVol; } } }
              { Name Edif; Value { Local { [ {u}  ]; In Omega; Jacobian JVol; } } }
							{ Name Edifx; Value { Local { [ CompX[{u}]  ]; In Omega; Jacobian JVol; } } }
							{ Name Edify; Value { Local { [ CompY[{u}]  ]; In Omega; Jacobian JVol; } } }
							{ Name Edifz; Value { Local { [ CompZ[{u}]  ]; In Omega; Jacobian JVol; } } }
              { Name E0 ; Value { Local { [    E0[]  ]; In Omega; Jacobian JVol; } } }
							{ Name E0x; Value { Local { [ CompX[E0[]]  ]; In Omega; Jacobian JVol; } } }
							{ Name E0y; Value { Local { [ CompY[E0[]]  ]; In Omega; Jacobian JVol; } } }
							{ Name E0z; Value { Local { [ CompZ[E0[]]  ]; In Omega; Jacobian JVol; } } }
							{ Name objective_local; Value { Local { [objective[{u}]]; In Omega; Jacobian JVol; } } }
							
              { Name source; Value { Local { [ source[]]; In Omega; Jacobian JVol; } } }

							{ Name sadj_int_x_re  ; Value { Integral { [Re[ CompX[adj_source_int[{u}]] ]] ; In Omega_target    ; Integration Int_1 ; Jacobian JVol ; } } }
							{ Name sadj_int  ; Value { Integral { [ adj_source_int[{u}] ] ; In Omega_target    ; Integration Int_1 ; Jacobian JVol ; } } }
							{ Name Adj   ; Value { Local { [ {u}   ] ; In Omega; Jacobian JVol; } } }
							{ Name sadj_int_re  ; Value { Integral { [Re[ adj_source_int[{u}] ]] ; In Omega_target    ; Integration Int_1 ; Jacobian JVol ; } } }
							{ Name sadj_int_im  ; Value { Integral { [Im[ adj_source_int[{u}] ]] ; In Omega_target    ; Integration Int_1 ; Jacobian JVol ; } } }
							{ Name int_objective  ; Value { Integral { [objective[{u}]] ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
							/* { Name source_adj_test  ; Value { Integral { [$source_adj] ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } } */
							/* { Name source_adj_test; Value { Local { [ $source_adj     ]; In Omega; Jacobian JVol; } } } */
							{ Name source_adj_test; Value { Local { [ComplexVectorField [XYZ[]]{2}    ]; In Omega; Jacobian JVol; } } }

              { Name po_r ; Value { Local { [ 0.5*Im[ 1/omega0 * 1/mu0*{u}/\Conj[{Curl u}] ] ]; In Omega; Jacobian JVol; } } }
						    //// Joule Losses Distribution = 1/2*sigma*|E_tot|^2
						  { Name normalized_losses1 ; Value { Integral { [  epsilon0*omega0 * 0.5*Fabs[Im[CompXX[epsilon[]]]]*(SquNorm[{u}+E0[]]) / (Pinc) ] ; In Omega_source ; Integration Int_1 ; Jacobian JVol ; } } }
							{ Name epsilon_annex; Value { Local { [ CompXX[epsilon_annex[]] ]; In Omega; Jacobian JVol; } } }

              { Name epsilon; Value { Local { [ CompXX[epsilon[]] ]; In Omega; Jacobian JVol; } } }
							{ Name dEq_deps   ; Value { Local { [dEq_deps[{u}]] ; In Omega; Jacobian JVol; } } }

}}}

// #############################################################################
PostOperation {
	{ Name postop_fields_pos; NameOfPostProcessing postpro_diff ;
			Operation {
				
		Print[ objective_local , OnElementsOf Omega, File "objective_local.pos"];
		Print[ Etot , OnElementsOf Omega, File "Etot.pos"];
		/* Print[ E0x , OnElementsOf Omega, File "E0x.pos"]; */
		/* Print[ E0y , OnElementsOf Omega, File "E0y.pos"]; */
		/* Print[ E0z , OnElementsOf Omega, File "E0z.pos"]; */
		/* Print[ E0 , OnElementsOf Omega, File "E0.pos"]; */
		/* Print[ Etotx , OnElementsOf Omega, File "Etotx.pos"]; */
		/* Print[ Etoty , OnElementsOf Omega, File "Etoty.pos"]; */
		/* Print[ Etotz , OnElementsOf Omega, File "Etotz.pos"]; */
		/* Print[ Edif , OnElementsOf Omega, File "E.pos"]; */
		/* Print[ Edifx , OnElementsOf Omega, File "Edifx.pos"]; */
		/* Print[ Edify , OnElementsOf Omega, File "Edify.pos"]; */
		/* Print[ Edifz , OnElementsOf Omega, File "Edifz.pos"]; */
	/* Print[ epsilon , OnElementsOf Omega, File "epsilon.pos"]; */
/* Print[ epsilon_annex , OnElementsOf Omega, File "epsilon_annex.pos"]; */
			}
	}

    { Name postopQ; NameOfPostProcessing postpro_diff ;
        Operation {
	    Print[ normalized_losses1[Omega_source] , OnGlobal, File "Q.txt", Format Table ];
        }
    }
		{ Name postop_epsilon; NameOfPostProcessing postpro_diff ;
				Operation {
			Print [ epsilon , OnElementsOf des, File "epsilon.pos"];
				}
		}

		{ Name postop_source_adj; NameOfPostProcessing postpro_diff ;
				Operation {
					Print[sadj_int, OnElementsOf Omega_target, StoreInField  2];
					/* Print[sadj_int_re, OnElementsOf Omega_target, StoreInField  2];
					Print[sadj_int_im, OnElementsOf Omega_target, StoreInField  3]; */
						/* Print[sadj_int_x_re, OnElementsOf Omega_target, File "sadj_int_x_re.pos"];
							Print[sadj_int_x_re, OnElementsOf Omega_target, StoreInField  4];
					 */
						/* Print[sadj_int_re, OnElementsOf Omega_target, File "sadj_int_re.pos"]; */


			}
		}

		/* { Name postop_source_adj_test; NameOfPostProcessing postpro_diff ;
				Operation {
							Print[source_adj_test, OnElementsOf Omega_target, File "source_adj_test.pos"];
			}
		} */

		{ Name postop_int_objective; NameOfPostProcessing postpro_diff ;
       Operation {
         /*Print [ u  , OnElementsOf Omega, File "u.pos" ];*/
       Print[ int_objective[Omega_target],  OnElementsOf PrintPoint, File "objective.txt" , Format SimpleTable ];
       }
    }

		{ Name postop_adjoint; NameOfPostProcessing postpro_diff ;
        Operation {
        If (nodes_flag)
          Print[Adj, OnElementsOf Omega_design ,  Format , LastTimeStepOnly, File "adjoint.txt" ];
        Else
          Print[Adj, OnElementsOf Omega_design , Depth 0, Format SimpleTable, LastTimeStepOnly, File "adjoint.txt" ];
        EndIf
      }
    }

		    { Name postop_dEq_deps; NameOfPostProcessing postpro_diff ;
		        Operation {
		        If (nodes_flag)
		          Print[dEq_deps, OnElementsOf Omega_design ,  Format NodeTable, File "dEq_deps.txt" ];
		        Else
		          Print[dEq_deps, OnElementsOf Omega_design ,  Depth 0, Format SimpleTable, File "dEq_deps.txt" ];
		        EndIf
		      }
		    }

}
// #############################################################################
