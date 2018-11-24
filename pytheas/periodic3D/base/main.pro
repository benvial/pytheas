Include "parameters.dat";
// #############################################################################
Group {
	// SubDomains
	PMLbot         = Region[1000];
	L_6            = Region[2000];
	L_5            = Region[3000];
	L_4            = Region[4000];
	L_3            = Region[5000];  // interp
	L_2            = Region[6000];
	L_1            = Region[7000];
	PMLtop         = Region[8000];
	// Boundaries
	SurfBlochXm    = Region[750];
	SurfBlochXp    = Region[751];
	SurfBlochYm    = Region[760];
	SurfBlochYp    = Region[761];
	SurfDirichlet  = Region[770];
	SurfBloch      = Region[{SurfBlochXm,SurfBlochXp,SurfBlochYm,SurfBlochYp}];
	// Lines
	LineNeumann    = Region[10001];
	// Domains
	Omega          = Region[{PMLbot,L_6,L_5,L_4,L_3,L_2,L_1,PMLtop}];
	Omega_nosource = Region[{PMLbot,L_6,L_1,PMLtop}];
	Omega_source   = Region[{L_2,L_3,L_4,L_5}];
	Omega_no_pml   = Region[{L_1,L_2,L_3,L_4,L_5, L_6}];
	Omega_super    = Region[{L_1,L_2,L_3,L_4,L_5,PMLtop}];
	Omega_subs     = Region[{L_6,PMLbot}];
	Omega_design   = Region[{L_3}];
	Omega_target   = Region[{L_2}];
	// Points
  PrintPoint	=  Region[10000];

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
	Pinc             =  0.5*Ae*Ae*Sqrt[epsilon0/mu0] * Cos[theta_0];

	epsilon[L_1]  = Complex[eps_L1_re , eps_L1_im] * TensorDiag[1,1,1];
	epsilon[L_2]  = Complex[eps_L2_re , eps_L2_im] * TensorDiag[1,1,1];
	epsilon[L_4]  = Complex[eps_L4_re , eps_L4_im] * TensorDiag[1,1,1];
	epsilon[L_5]  = Complex[eps_L5_re , eps_L5_im] * TensorDiag[1,1,1];
	epsilon[L_6]  = Complex[eps_L6_re , eps_L6_im] * TensorDiag[1,1,1];
	epsilon[Omega_design]    = Complex[ScalarField[XYZ[], 0, 1]{0}  ,ScalarField[XYZ[], 0, 1 ]{1} ]  * TensorDiag[1,1,1];
	/* epsilon[Omega_design]          = Complex[eps_L3_re , eps_L3_im] * TensorDiag[1,1,1]; */


	mu[Omega_no_pml]       = TensorDiag[1,1,1];
	nu[Omega_no_pml]       = TensorDiag[1,1,1];
/*
	mu[L_3]       = TensorDiag[1,1,1];
	nu[L_3]       = TensorDiag[1,1,1]; */

	gamma[Omega_super] = Complex[gamma_super_re , gamma_super_im];
	gamma[Omega_subs]  = Complex[gamma_subs_re  , gamma_subs_im ];
	EXpj[Omega_subs]   = Complex[Expj_subs_re , Expj_subs_im];
	EXmj[Omega_subs]   = Complex[Exmj_subs_re , Exmj_subs_im];
	EYpj[Omega_subs]   = Complex[Eypj_subs_re , Eypj_subs_im];
	EYmj[Omega_subs]   = Complex[Eymj_subs_re , Eymj_subs_im];
	EZpj[Omega_subs]   = Complex[Ezpj_subs_re , Ezpj_subs_im];
	EZmj[Omega_subs]   = Complex[Ezmj_subs_re , Ezmj_subs_im];
	EXpj[Omega_super]  = Complex[Expj_super_re , Expj_super_im];
	EXmj[Omega_super]  = Complex[Exmj_super_re , Exmj_super_im];
	EYpj[Omega_super]  = Complex[Eypj_super_re , Eypj_super_im];
	EYmj[Omega_super]  = Complex[Eymj_super_re , Eymj_super_im];
	EZpj[Omega_super]  = Complex[Ezpj_super_re , Ezpj_super_im];
	EZmj[Omega_super]  = Complex[Ezmj_super_re , Ezmj_super_im];

	Propp[]         = Ae * Complex[ Cos[alpha0*X[]+beta0*Y[]+gamma[]*(Z[]-0.0)] , Sin[alpha0*X[]+beta0*Y[]+gamma[]*(Z[]-0.0)] ];
	Propm[]         = Ae * Complex[ Cos[alpha0*X[]+beta0*Y[]-gamma[]*(Z[]-0.0)] , Sin[alpha0*X[]+beta0*Y[]-gamma[]*(Z[]-0.0)] ];
	EXp[]           = EXpj[]*Propp[];
	EXm[]           = EXmj[]*Propm[];
	EYp[]           = EYpj[]*Propp[];
	EYm[]           = EYmj[]*Propm[];
	EZp[]           = EZpj[]*Propp[];
	EZm[]           = EZmj[]*Propm[];
	HXm[]           = -1/(omega0*mu0)*( beta0  * EZm[]  + gamma[]* EYm[]);
	HYm[]           = -1/(omega0*mu0)*(-gamma[]* EXm[]  - alpha0 * EZm[]);
	HZm[]           = -1/(omega0*mu0)*( alpha0 * EYm[]  - beta0  * EXm[]);
	HXp[]           = -1/(omega0*mu0)*( beta0  * EZp[]  - gamma[]* EYp[]);
	HYp[]           = -1/(omega0*mu0)*( gamma[]* EXp[]  - alpha0 * EZp[]);
	HZp[]           = -1/(omega0*mu0)*( alpha0 * EYp[]  - beta0  * EXp[]);
	EXcm[]          = EXp[]+EXm[];
	EYcm[]          = EYp[]+EYm[];
	EZcm[]          = EZp[]+EZm[];
	HXcm[]          = HXm[]+HXp[];
	HYcm[]          = HYm[]+HYp[];
	HZcm[]          = HZm[]+HZp[];

	Ep[]            = Vector[EXp[],EYp[],EZp[]];
	Em[]            = Vector[EXm[],EYm[],EZm[]];
	Ecm[]           = Vector[EXcm[],EYcm[],EZcm[]];
	dephX[]         = Complex[ Cos[alpha0*period_x] , Sin[alpha0*period_x] ];
	dephY[]         = Complex[ Cos[beta0 *period_y] , Sin[beta0 *period_y] ];

// PML parameters
	sx              = 1.;
	sz[]            = Complex[a_pml,-b_pml];
	sy              = 1.;

	epsilon[PMLbot] = eps_L6_re*TensorDiag[sz[]*sy/sx,sx*sz[]/sy,sx*sy/sz[]];
	epsilon[PMLtop] = eps_L1_re*TensorDiag[sz[]*sy/sx,sx*sz[]/sy,sx*sy/sz[]];
	mu[PMLbot]      =            TensorDiag[sz[]*sy/sx,sx*sz[]/sy,sx*sy/sz[]];
	mu[PMLtop]      =            TensorDiag[sz[]*sy/sx,sx*sz[]/sy,sx*sy/sz[]];
	nu[PMLbot]      =            TensorDiag[sx/(sz[]*sy),sy/(sx*sz[]),sz[]/(sx*sy)];
	nu[PMLtop]      =            TensorDiag[sx/(sz[]*sy),sy/(sx*sz[]),sz[]/(sx*sy)];


	epsilon_annex[PMLbot]       = eps_L6_re*TensorDiag[sz[]*sy/sx,sx*sz[]/sy,sx*sy/sz[]];
	epsilon_annex[PMLtop]       = eps_L1_re*TensorDiag[sz[]*sy/sx,sx*sz[]/sy,sx*sy/sz[]];
	epsilon_annex[Omega_source] = Complex[eps_L1_re , eps_L1_im] * TensorDiag[1,1,1];
	epsilon_annex[L_1]          = Complex[eps_L1_re , eps_L1_im] * TensorDiag[1,1,1];
	epsilon_annex[L_6]          = Complex[eps_L6_re , eps_L6_im] * TensorDiag[1,1,1];

	source[]         = k0^2*(epsilon[]-epsilon_annex[])*Ecm[];//(nm/1.e-9)^2*

	// postpro diffraction efficiencies
  For i In {0:nb_slice-1}
		zcut_sub~{i} = hh_L6+thick_L6-scan_dist - i*(thick_L6-2*scan_dist)/(nb_slice-1);
		zcut_sup~{i} = hh_L1+scan_dist           + i*(thick_L1-2*scan_dist)/(nb_slice-1);
	EndFor


	// Topology optimization

	  coef_obj[] =  1/(period_x*period_y*thick_L2);
		objective[] = coef_obj[] * SquNorm[$1 + Ecm[] ];
		adj_source_int[] = -2 * coef_obj[] * Conj[$1 + Ecm[]];
		db_deps[] = -k0^2*Ecm[];
		dA_deps[] = k0^2;
		dEq_deps[] = db_deps[] - dA_deps[] * ($1);

		source_adj[] = Vector[0,0,0];



}
// #############################################################################
Constraint {
        {Name Dirichlet; Type Assign;
                Case {
                        { Region SurfDirichlet; Value 0.; }
                }
        }
        {Name BlochX;
                Case {
                        { Region SurfBlochXp; Type LinkCplx ; RegionRef SurfBlochXm;
                        Coefficient dephX[]; Function Vector[$X-period_x,$Y,$Z] ;
                        }
                }
}
        {Name BlochY;
                Case {
                        { Region SurfBlochYp; Type LinkCplx ; RegionRef SurfBlochYm;
                        Coefficient dephY[]; Function Vector[$X,$Y-period_y,$Z] ;
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
// #############################################################################
FunctionSpace {
  { Name Hcurl; Type Form1;
    BasisFunction {
      { Name sn; NameOfCoef un; Function BF_Edge;
        Support Region[{Omega}]; Entity EdgesOf[All]; }
     { Name sn2; NameOfCoef un2; Function BF_Edge_2E;
        Support Region[{Omega}]; Entity EdgesOf[All]; }
				If (el_order== 2)
							{ Name sn3; NameOfCoef un3; Function BF_Edge_3F_b;
				         Support Region[Omega]; Entity FacetsOf[Omega, Not SurfBloch]; }
				      { Name sn4; NameOfCoef un4; Function BF_Edge_3F_c;
				         Support Region[Omega]; Entity FacetsOf[Omega, Not SurfBloch]; }
				      { Name sn5; NameOfCoef un5; Function BF_Edge_4E;
				         Support Region[Omega]; Entity EdgesOf[Omega, Not SurfBloch]; }
				 EndIf
    }
    Constraint {
      { NameOfCoef un;  EntityType EdgesOf ; NameOfConstraint BlochX; }
      { NameOfCoef un;  EntityType EdgesOf ; NameOfConstraint BlochY; }
      { NameOfCoef un;  EntityType EdgesOf ; NameOfConstraint Dirichlet; }
       { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint BlochX; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint BlochY; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }

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

              { Name Ecal; Value { Local { [ {u}       ]; In Omega; Jacobian JVol; } } }
              { Name Etot; Value { Local { [ {u}+Ecm[] ]; In Omega; Jacobian JVol; } } }
              { Name Edif; Value { Local { [ {u}+Em[]  ]; In Omega; Jacobian JVol; } } }
              { Name Ecm ; Value { Local { [    Ecm[]  ]; In Omega; Jacobian JVol; } } }
              { Name source; Value { Local { [ source[]]; In Omega; Jacobian JVol; } } }
              { Name Ecal_x; Value { Local { [ CompX[{u}] ]; In Omega; Jacobian JVol; } } }
              { Name Ecal_y; Value { Local { [ CompY[{u}] ]; In Omega; Jacobian JVol; } } }
              { Name Ecal_z; Value { Local { [ CompZ[{u}] ]; In Omega; Jacobian JVol; } } }
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
						  { Name normalized_losses1 ; Value { Integral { [  epsilon0*omega0 * 0.5*Fabs[Im[CompXX[epsilon[]]]]*(SquNorm[{u}+Ecm[]]) / (Pinc*period_x*period_y) ] ; In Omega_source ; Integration Int_1 ; Jacobian JVol ; } } }

              { Name epsilon; Value { Local { [ CompXX[epsilon[]] ]; In Omega; Jacobian JVol; } } }
							{ Name dEq_deps   ; Value { Local { [dEq_deps[{u}]] ; In Omega; Jacobian JVol; } } }

}}}

// #############################################################################
PostOperation {
	{ Name postop_fields_pos; NameOfPostProcessing postpro_diff ;
			Operation {
		Print[ Etot , OnElementsOf Omega, File "Etot.pos"];
		Print[ Ecal , OnElementsOf Omega, File "E.pos"];
			}
	}

    { Name postopQ; NameOfPostProcessing postpro_diff ;
        Operation {
	    Print[ normalized_losses1[Omega_source] , OnGlobal, File "Q.txt", Format Table ];
        }
    }
		{ Name postop_epsilon; NameOfPostProcessing postpro_diff ;
				Operation {
			Print [ epsilon , OnElementsOf L_3, File "epsilon.pos"];
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


{ Name Ed; NameOfPostProcessing postpro_diff ;
    Operation {
		For i In {0:nb_slice-1}
		Print [ Etot , OnPlane { {-period_x/2,-period_y/2,zcut_sub~{i}} {-period_x/2,period_y/2,zcut_sub~{i}} {period_x/2,-period_y/2,zcut_sub~{i}} } {ninterv_integ,ninterv_integ} , File > "Etot_XYcut.out",Smoothing ,Format SimpleTable];
  		Print [ Edif , OnPlane { {-period_x/2,-period_y/2,zcut_sup~{i}} {-period_x/2,period_y/2,zcut_sup~{i}} {period_x/2,-period_y/2,zcut_sup~{i}} } {ninterv_integ,ninterv_integ} , File > "Edif_XYcut.out",Smoothing, Format SimpleTable];
		EndFor
        }
}
}
// #############################################################################
