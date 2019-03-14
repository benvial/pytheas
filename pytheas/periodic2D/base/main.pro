Include "parameters.dat";


// #############################################################################
// #############################################################################

Group {
    // Domains
    pmlbot          = Region[1000];
    sub             = Region[2000];
    layer1          = Region[3000];
    design          = Region[4000];
    layer2          = Region[5000];
    sup             = Region[6000];
    pmltop          = Region[7000];
    If (inclusion_flag)
      incl          = Region[8000];
      Omega_source    = Region[{layer1,design,layer2, incl}];
      Omega_i = Region[{sup , layer1, design, layer2, incl}];
      Omega_i_bar = Region[{pmltop , pmlbot, sub}];
      sub_bar = Region[{pmltop, pmlbot, sup , layer1, design, layer2, incl}];
    Else
      Omega_source    = Region[{layer1,design,layer2}];
      Omega_i = Region[{sup , layer1, design, layer2}];
      Omega_i_bar = Region[{pmltop , pmlbot, sub}];
      sub_bar = Region[{pmltop, pmlbot, sup , layer1, design, layer2}];
    EndIf
    Omega_nosource  = Region[{pmltop,pmlbot,sup,sub}];
    Omega           = Region[{Omega_source,Omega_nosource}];




    /* Omega_target    = Region[{layer1,layer2, design}]; */
    Omega_target    = Region[{sup}];
    Omega_design    = Region[{design}];
    // Boundaries
    SurfBlochLeft   = Region[101];
    SurfBlochRight  = Region[102];
    SurfDirichlet   = Region[110];
    SurfContinuity  = Region[120];
    // Points
    PrintPoint	=  Region[10000];
}
// #############################################################################

Function{
    j[] = Complex[0.0, 1.0];
    Freq   = cel/lambda0;
    omega0 = 2.*Pi*cel/lambda0;
    k0     = 2.*Pi/lambda0;
    ksearch = 2.*Pi/lambda0search;
    k_sup  = 2.*Pi*Sqrt[eps_sup_re]/lambda0;
    k_sub  = 2.*Pi*Sqrt[eps_sub_re]/lambda0;
    alpha_sup = -k_sup*Sin[theta];
    beta_sup  =   k_sup*Cos[theta];
    beta_sub  = Sqrt[k_sub*k_sub-alpha_sup*alpha_sup];
    PW[] = $1* Exp[ j[]* ( $2 * X[] + $3 * Y[] ) ] ;  // plane wave PW[amplitude, kx, ky]
    If (TE_flag)
        R[]   =  (beta_sup-beta_sub)/(beta_sup+beta_sub);
        T[]   =  (2.*beta_sup)/(beta_sup+beta_sub);
        Pinc  =  0.5*A*A*Sqrt[epsilon0/mu0] * Cos[theta];
    Else
        beta_S_sup[] = beta_sup/Complex[eps_sup_re,eps_sup_im];
        beta_S_sub[] = beta_sub/Complex[eps_sub_re,eps_sub_im];
        R[]          = (beta_S_sup[]-beta_S_sub[])/(beta_S_sup[]+beta_S_sub[]);
        T[]          = (2.*beta_S_sup[])/(beta_S_sup[]+beta_S_sub[]);
        Pinc         =  0.5*A*A*Sqrt[mu0/epsilon0] * Cos[theta];
    EndIf

    deph[] = Complex[ Cos[alpha_sup*d], Sin[alpha_sup*d] ];

    // PML parameters
    sx               =   1.;
    sy[]             =   Complex[a_pml, -b_pml];
    sz               =   1.;

    // Permittivities
    epsilonr[sup]            = Complex[eps_sup_re, eps_sup_im] * TensorDiag[1,1,1];
    epsilonr[sub]            = Complex[eps_sub_re,eps_sub_im] * TensorDiag[1,1,1];
    epsilonr[layer1]         = Complex[eps_layer1_re,eps_layer1_im] * TensorDiag[1,1,1];

    epsilonr[layer2]         = Complex[eps_layer2_re,eps_layer2_im] * TensorDiag[1,1,1];
    epsilonr[pmltop]         = eps_sup_re*TensorDiag[sz*sy[]/sx,sx*sz/sy[],sx*sy[]/sz];
    epsilonr[pmlbot]         = eps_sub_re*TensorDiag[sz*sy[]/sx,sx*sz/sy[],sx*sy[]/sz];

    If (aniso)
        epsilonr_xx[]  = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}];
        epsilonr_yy[]  = Complex[ScalarField[XYZ[], 0, 1 ]{2}, ScalarField[XYZ[], 0, 1 ]{3}];
        epsilonr_zz[]  = Complex[ScalarField[XYZ[], 0, 1 ]{4}, ScalarField[XYZ[], 0, 1 ]{5}];
    EndIf

    If (inclusion_flag)

      epsilonr[incl]           = Complex[eps_incl_re,eps_incl_im] * TensorDiag[1,1,1];
      epsilonr_annex[incl]   = Complex[eps_sup_re,eps_sup_im] * TensorDiag[1,1,1];
      If (aniso)
          epsilonr[design] =  TensorDiag[epsilonr_xx[],epsilonr_yy[],epsilonr_zz[]];
        Else
          epsilonr[design]         = Complex[eps_des_re,eps_des_im] * TensorDiag[1,1,1];
        EndIf
    Else

    If (aniso)
        epsilonr[design] =  TensorDiag[epsilonr_xx[],epsilonr_yy[],epsilonr_zz[]];
      Else
      epsilonr[design]         = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}] * TensorDiag[1,1,1];

      EndIf

    EndIf

    epsilonr_annex[sup]      = Complex[eps_sup_re,eps_sup_im] * TensorDiag[1,1,1];
    epsilonr_annex[sub]      = Complex[eps_sub_re,eps_sub_im] * TensorDiag[1,1,1];
    epsilonr_annex[layer1]   = Complex[eps_sup_re,eps_sup_im] * TensorDiag[1,1,1];
    epsilonr_annex[design]   = Complex[eps_sup_re,eps_sup_im] * TensorDiag[1,1,1];
    epsilonr_annex[layer2]   = Complex[eps_sup_re,eps_sup_im] * TensorDiag[1,1,1];
    epsilonr_annex[pmltop]   = eps_sup_re*TensorDiag[sz*sy[]/sx,sx*sz/sy[],sx*sy[]/sz];
    epsilonr_annex[pmlbot]   = eps_sub_re*TensorDiag[sz*sy[]/sx,sx*sz/sy[],sx*sy[]/sz];

    // Permeabilities
    mur[pmltop]              = TensorDiag[sz*sy[]/sx,sx*sz/sy[],sx*sy[]/sz];
    mur[sup]                 = TensorDiag[1,1,1];
    mur[layer1]              = TensorDiag[1,1,1];
    mur[layer2]              = TensorDiag[1,1,1];
    mur[design]              =  Complex[mu_des_re,mu_des_im] *TensorDiag[1,1,1];
    mur[sub]                 = TensorDiag[1,1,1];
    mur[pmlbot]              = TensorDiag[sz*sy[]/sx,sx*sz/sy[],sx*sy[]/sz];
    mur_annex[Omega]         = TensorDiag[1,1,1];
    If (inclusion_flag)
      mur[incl]                 =  Complex[mu_incl_re,mu_incl_im] * TensorDiag[1,1,1];
    EndIf



    // Fields
    u_i[Omega_i] = PW[A, alpha_sup, beta_sup];
    u_i[Omega_i_bar]   = 0.;

    u_r[Omega_i] = PW[R[] , alpha_sup, - beta_sup];
    u_r[Omega_i_bar]   = 0.;

    u_t[sub]   =  PW[ T[] , alpha_sup,  beta_sub];
    u_t[sub_bar] = 0.;

    u_1[] 	      = u_i[]+u_r[]+u_t[];
    u_1_d[]       = u_r[]+u_t[];

    r_ref[] = rtar;//Complex[0.9, 0];
    /* r_ref[] = 0.5; */

    u_ref[Omega] = PW[r_ref[], alpha_sup, - beta_sup];


    If (TE_flag)
        dual_i[]      = Vector[-beta_sup, alpha_sup, 0] *u_i[] / (omega0*mu0*CompXX[mur_annex[]])  ;
        dual_t[]      = Vector[-beta_sub, alpha_sup, 0] *u_t[] / (omega0*mu0*CompXX[mur_annex[]])  ;
        dual_r[]      = Vector[ beta_sup, alpha_sup, 0] *u_r[] / (omega0*mu0*CompXX[mur_annex[]])  ;
        source_eps[]      =  k0^2*(epsilonr[]-epsilonr_annex[])*u_1[];
        xi[] = TensorDiag[CompYY[mur[]],CompXX[mur[]],CompZZ[mur[]]];
        source_r[]    =  j[]*(1/mur_annex[]-1/xi[])* u_r[]*TensorDiag[alpha_sup, -beta_sup,0.];
        source_i[]    =  j[]*(1/mur_annex[]-1/xi[])* u_i[]*TensorDiag[alpha_sup, beta_sup,0.];
        source_mu[]      =  source_r[]+ source_i[];
        weight[] = epsilonr[];
    Else
        dual_i[]      = Vector[ beta_sup, -alpha_sup, 0] *u_i[] / (omega0*epsilon0*CompXX[epsilonr_annex[]]);
        dual_t[]      = Vector[ beta_sub,  -alpha_sup, 0] *u_t[] / (omega0*epsilon0*CompXX[epsilonr_annex[]]);
        dual_r[]      = Vector[-beta_sup,  -alpha_sup, 0] *u_r[] / (omega0*epsilon0*CompXX[epsilonr_annex[]]);
        source_mu[]      =  k0^2*(mur[]-mur_annex[])*u_1[];
        xi[] = TensorDiag[CompYY[epsilonr[]],CompXX[epsilonr[]],CompZZ[epsilonr[]]];
        source_r[]    =  j[]*(1/epsilonr_annex[]-1/xi[])* u_r[]*TensorDiag[alpha_sup, -beta_sup,0.];
        source_i[]    =  j[]*(1/epsilonr_annex[]-1/xi[])* u_i[]*TensorDiag[alpha_sup, beta_sup,0.];
        source_eps[]      =  source_r[]+ source_i[];
        weight[] = mur[];
    EndIf

    dual_1_d[] 	       =  dual_r[] + dual_t[];
    dual_1[] 	         =  dual_i[] + dual_1_d[];

    For i In {0:nb_slice-1}
      ycut_sub~{i} = ycut_sub_min + i*(ycut_sub_max-ycut_sub_min)/(nb_slice-1);
      ycut_sup~{i} = ycut_sup_min + i*(ycut_sup_max-ycut_sup_min)/(nb_slice-1);
    EndFor




    If (TE_flag)
          coef_Q[] = 0.5 * epsilon0*omega0*Fabs[Im[CompZZ[epsilonr[]]]   ]  / (Pinc*d);
          dual[] = Vector[ CompY[ $1 ], - CompX[ $1 ], 0 ]/(j[]*omega0*mu0*CompXX[mur[]]); // TE case H = dual[{d u}]
          dual_tot[] = dual[$1] + dual_1[] * CompXX[mur_annex[]] /CompXX[mur[]] ;
          absorption[]  =  coef_Q[] * SquNorm[$1 + u_1[]] ;
    Else
          dual[] = Vector[ CompY[ $1 ], - CompX[ $1 ], 0 ]/(-j[]*omega0*epsilon0*CompXX[epsilonr[]]); // TM case E = dual[{d u}]
          dual_tot[] = -dual[$1] + dual_1[] * CompXX[epsilonr_annex[]] /CompXX[epsilonr[]] ;
          coef_Q[] = 0.5 * epsilon0*omega0  / (Pinc*d);
          absorption[]  = - coef_Q[] * ( Im[CompXX[ epsilonr[]]]  * SquNorm[CompX[dual_tot[$1]] ] + Im[CompYY[ epsilonr[]]] *SquNorm[CompY[dual_tot[$1]] ] );

    EndIf

    // Topology optimization

    /* DefineFunction[ adj_source ]; */

    If (TE_flag)
      /* coef_obj[] = coef_Q[];  // 1/(h_sub*d); */
      /* objective[] = absorption[$1];
      adj_source_int[] = - 2 * coef_obj[] *  Conj[($1 + u_1[]) ];  */

      coef_obj[] =  1/(h_sup*d);


      /* objective[] = coef_obj[] *SquNorm[$1 + u_1_d[] - u_ref[]] ;
      adj_source_int[] = - 2 * coef_obj[] *  Conj[($1 + u_1_d[] - u_ref[]) ]; */


      objective[] = coef_obj[] * SquNorm[SquNorm[$1 + u_1_d[]] - SquNorm[u_ref[]]] ;
      adj_source_int[] = - 4 * coef_obj[] * (SquNorm[$1 + u_1_d[]] - SquNorm[u_ref[]]) * Conj[$1 + u_1_d[]] ;


      db_deps[] = -k0^2*u_1[];
      dA_deps[] = k0^2;
      dEq_deps[] = db_deps[] - dA_deps[] * ($1);
    Else
    coef_obj[] =  1/(h_sub*d);
      objective[] = absorption[$1] ;
      adj_source_int[] = - 2 * coef_obj[] *  Conj[($1 + u_1[]) ];
      beta_Q[] = coef_obj[] /(omega0*epsilon0* (Re[CompXX[epsilonr[]]]^2 + Im[CompXX[epsilonr[]]]^2 ) );
      adj_source[] =  - 2 * beta_Q[] * ( Conj[$1 + dual_1[] * CompXX[epsilonr_annex[]] /CompXX[epsilonr[]] ]  ) ; //*h_des*d / ElementVol[]; //d_objective_du *ElementVol[]
      db_deps[] = -j[]/CompXX[epsilonr[]]^2* u_i[]*TensorDiag[alpha_sup, beta_sup,0.];
      dA_deps[] = -1/CompXX[epsilonr[]]^2;
      dEq_deps[] = db_deps[] - dA_deps[] * ($1);




    EndIf

// Electrostatics
/*source_electrostat[Region[{layer2}]] = 1e3 ;*/
h2 = h_layer2/2 + h_des + h_layer1;
h1 = h_layer1/2;
/*source_electrostat[Region[{layer1}]] = 1e3 * Exp[-$X^2/0.1^2 - ($Y-h1)^2/0.1^2];
source_electrostat[Region[{layer2}]] = -1e3 * Exp[-$X^2/0.1^2 - ($Y-h2)^2/0.1^2];
source_electrostat[Region[{pmltop,sup, design, sub, pmlbot}]] =0;*/
source_electrostat[] = 0;

}

// #############################################################################

Constraint {
    {Name Dirichlet; Type Assign;
        Case {
            { Region SurfDirichlet; Value 0.; }
        }
    }
    {Name Bloch;
        Case {
            { Region SurfBlochRight; Type LinkCplx ; RegionRef SurfBlochLeft; Coefficient deph[]; Function Vector[$X-d,$Y,$Z] ;
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
  { Name Hgradu; Type Form0;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_Node;    Support Region[Omega]; Entity NodesOf[Omega]; }
      { Name sn2; NameOfCoef un2; Function BF_Node_2E; Support Region[Omega]; Entity EdgesOf[Omega]; }
     }
    Constraint {
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Bloch; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Bloch; }
   }
  }
 { Name Hgraduadj; Type Form0;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_Node;    Support Region[Omega]; Entity NodesOf[Omega]; }
      { Name sn2; NameOfCoef un2; Function BF_Node_2E; Support Region[Omega]; Entity EdgesOf[Omega]; }
     }
    Constraint {
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Bloch; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Bloch; }
   }
  }
  { Name Hgradfilter; Type Form0;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_Node;    Support Region[Omega_design]; Entity NodesOf[Omega_design]; }
      { Name sn2; NameOfCoef un2; Function BF_Node_2E; Support Region[Omega_design]; Entity EdgesOf[Omega_design]; }
     }
    Constraint {
   }
  }
}

// #############################################################################

Formulation{
  /*------------ Diffraction problem -----------------*/
    {Name helmoltz_scalar; Type FemEquation;
        Quantity {
        { Name u; Type Local; NameOfSpace Hgradu;}
        }
        Equation {
        If (TE_flag)
            Galerkin { [k0^2*CompZZ[epsilonr[]]*Dof{u} , {u}];
            In Omega; Jacobian JVol; Integration Int_1;  }
            Galerkin { [-1/xi[]*Dof{d u} , {d u}];
            In Omega; Jacobian JVol; Integration Int_1; }
            Galerkin { [ ($Source ? source_eps[] : 0) , {u}];
            In Omega_source; Jacobian JVol; Integration Int_1;  }
            Galerkin { [ ($Source ? source_mu[] : 0) , {d u}];
            In Omega_source; Jacobian JVol; Integration Int_1;  }
            Galerkin { [ ($SourceAdj ? Complex[ScalarField[XYZ[], 0, 1 ]{2}, ScalarField[XYZ[], 0, 1 ]{3}] : 0) , {u}];
            In Omega_target; Jacobian JVol; Integration Int_1;  }
            /* Galerkin { [ ($SourceAdj ? ElementVol[] * adj_source_int[{u}] : 0) , {u}];
            In Omega_target; Jacobian JVol; Integration Int_1;  } */
        Else
            Galerkin { [k0^2*CompZZ[mur[]]*Dof{u} , {u}];
            In Omega; Jacobian JVol; Integration Int_1;  }
            Galerkin { [-1/xi[]*Dof{d u} , {d u}];
            In Omega; Jacobian JVol; Integration Int_1; }
            Galerkin { [ ($Source ? source_mu[] : 0) , {u}];
            In Omega_source; Jacobian JVol; Integration Int_1;  }
            Galerkin { [ ($Source ? source_eps[] : 0) , {d u}];
            In Omega_source; Jacobian JVol; Integration Int_1;  }
            Galerkin { [ ($SourceAdj ? Complex[ScalarField[XYZ[], 0, 1 ]{2}, ScalarField[XYZ[], 0, 1 ]{3}] : 0) , {d u}];
            In Omega; Jacobian JVol; Integration Int_1;  }
        EndIf
        }
    }
  /* ----------- Modal analysis --------------------*/
  {Name helmoltz_scalar_modal; Type FemEquation;
      Quantity {
      { Name u; Type Local; NameOfSpace Hgradu;}
      }
      Equation {
      If (TE_flag)
          Galerkin {  DtDtDof[ -CompZZ[epsilonr[]]*Dof{u} , {u}];
          In Omega; Jacobian JVol; Integration Int_1;  }
          Galerkin { [-1/TensorDiag[CompYY[mur[]],CompXX[mur[]],CompXX[mur[]]]*Dof{d u} , {d u}];
          In Omega; Jacobian JVol; Integration Int_1; }
      Else
          Galerkin { DtDtDof[ -CompZZ[mur[]]*Dof{u} , {u}];
          In Omega; Jacobian JVol; Integration Int_1;  }
          Galerkin { [-1/TensorDiag[CompYY[epsilonr[]],CompXX[epsilonr[]],CompXX[epsilonr[]]]*Dof{d u} , {d u}];
          In Omega; Jacobian JVol; Integration Int_1; }
      EndIf
      }
  }
  /*------------ Electrostatics problem -----------------*/
    {Name electrostat; Type FemEquation;
        Quantity {
        { Name u; Type Local; NameOfSpace Hgradu;}
        }
        Equation {
            Galerkin { [-CompZZ[epsilonr[]]*Dof{d u} , {d u}];
            In Omega; Jacobian JVol; Integration Int_1; }
            Galerkin { [ source_electrostat[] , {u}];
            In Omega_source; Jacobian JVol; Integration Int_1;  }

        }
    }

}

// #############################################################################

Resolution {
  /*------------ Diffraction problem -----------------*/
    { Name helmoltz_scalar;
        System {
        { Name S; NameOfFormulation helmoltz_scalar; Type ComplexValue; Frequency Freq;}
        }
        Operation {
        Evaluate[$Source = 1, $SourceAdj = 0];
        Printf [ "Helmoltz equation" ];
        Generate[S] ;Solve[S] ;SaveSolution[S] ;

        If (adjoint_flag)
          PostOperation[postop_int_objective];
          PostOperation[postop_dEq_deps];
          PostOperation[postop_source_adj];

          Evaluate[$Source = 0, $SourceAdj = 1];
          /* Generate[S] ;Solve[S] ;SaveSolution[S] ; */
          GenerateRHSGroup[S, Omega_target]; SolveAgain[S] ; SaveSolution[S] ;

          PostOperation[postop_adjoint];
        EndIf

        }

    }
  /* ----------- Modal analysis --------------------*/
  { Name helmoltz_scalar_modal;
      System {
      { Name Smodal; NameOfFormulation helmoltz_scalar_modal; Type ComplexValue;}
      }
      Operation {
      GenerateSeparate[Smodal]; EigenSolve[Smodal,neig,ksearch^2,0]; SaveSolutions[Smodal];
      }

  }
  /*------------ electrostat problem -----------------*/
    { Name electrostat;
        System {
        { Name S; NameOfFormulation electrostat; Type ComplexValue;}
        }
        Operation {
        Printf [ "Electrostatic equation" ];
        Generate[S] ;Solve[S] ;SaveSolution[S] ;

        }

    }
}

// #############################################################################

PostProcessing {
    /*------------ Diffraction problem -----------------*/
    { Name postpro; NameOfFormulation helmoltz_scalar;
        Quantity {
            { Name dEq_deps   ; Value { Local { [dEq_deps[{u}]] ; In Omega; Jacobian JVol; } } }
            { Name check_adj_source   ; Value { Local { [adj_source_int[{u}]] ; In Omega; Jacobian JVol; } } }
            { Name objective   ; Value { Local { [objective[{u}]] ; In Omega; Jacobian JVol; } } }
            { Name int_objective  ; Value { Integral { [objective[{u}]] ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
            { Name u   ; Value { Local { [ {u}    ] ; In Omega; Jacobian JVol; } } }
            { Name u1   ; Value { Local { [ u_1[]    ] ; In Omega; Jacobian JVol; } } }
            { Name epsilonr   ; Value { Local { [ CompZZ[epsilonr[] ]   ] ; In Omega; Jacobian JVol; } } }
            { Name u_diff   ; Value { Local { [ {u}+u_1_d[]     ] ; In Omega; Jacobian JVol; } } }
            { Name u_tot    ; Value { Local { [ {u}+u_1[]       ] ; In Omega; Jacobian JVol; } } }
            { Name u_tot_sqnorm    ; Value { Local { [ ({u}+u_1[])*Conj[({u}+u_1[]) ]       ] ; In Omega; Jacobian JVol; } } }
            { Name vx_diff   ; Value { Local { [ CompX[dual_1_d[] + dual[{d u}] ]      ] ; In Omega; Jacobian JVol; } } }
            { Name vy_diff    ; Value { Local { [  CompY[dual_1_d[] + dual[{d u}]  ]    ] ; In Omega; Jacobian JVol; } } }
            { Name vx_tot    ; Value { Local { [  CompX[dual_1[] + dual[{d u}]  ]         ] ; In Omega; Jacobian JVol; } } }
            { Name vy_tot   ; Value { Local { [  CompY[dual_1[]   + dual[{d u}]  ]  ] ; In Omega; Jacobian JVol; } } }
            { Name v_tot   ; Value { Local { [  dual_1[]   + dual[{d u}]   ] ; In Omega; Jacobian JVol; } } }
            { Name u_adj   ; Value { Local { [ {u}   ] ; In Omega; Jacobian JVol; } } }
            { Name sadj_int_re  ; Value { Integral { [Re[ adj_source_int[{u}] ]] ; In Omega_target    ; Integration Int_1 ; Jacobian JVol ; } } }
            { Name sadj_int_im  ; Value { Integral { [Im[ adj_source_int[{u}] ]] ; In Omega_target    ; Integration Int_1 ; Jacobian JVol ; } } }
            If (TE_flag)
                { Name Q  ; Value { Integral { [ absorption[{u}] ] ; In Omega_source    ; Integration Int_1 ; Jacobian JVol ; } } }
                { Name abso_density   ; Value { Local { [ absorption[{u}] ]; In Omega_source; Jacobian JVol; } } }
            Else
                { Name Q ; Value { Integral { [ absorption[{d u}] ] ; In Omega_source ; Integration Int_1 ; Jacobian JVol ; } } }
                { Name abso_density   ; Value { Local {[ absorption[{d u}] ]; In Omega_source; Jacobian JVol; } } }
            EndIf
            /*If (adjoint_flag)
                { Name obj  ; Value { Integral { [  obj_dens[] ] ; In Omega_target    ; Integration Int_1 ; Jacobian JVol ; } } }
            EndIf*/


            }

    }

    /* ----------- Modal analysis --------------------*/
    { Name postpro_modal; NameOfFormulation helmoltz_scalar_modal;
      Quantity {
        { Name EigenValues;  Value { Local{ [$EigenvalueImag]; In PrintPoint; Jacobian JVol; } } }
        { Name EigenVectors ;   Value { Local { [     {u}  ] ; In Omega; Jacobian JVol; } } }
        { Name Normalization ; Value { Integral { [  CompZZ[weight[]]*{u}*{u} ] ; In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }

            }
    }

    /* ----------- electrostatic analysis --------------------*/
    { Name postpro_electrostat; NameOfFormulation electrostat;
      Quantity {
      { Name source_electrostat   ; Value { Local { [ source_electrostat[]    ] ; In Omega; Jacobian JVol; } } }
      { Name u   ; Value { Local { [ {u}    ] ; In Omega; Jacobian JVol; } } }
      { Name e; Value { Local { [ -{d u} ]; In Omega; Jacobian JVol; } } }
      { Name norm_e; Value { Local { [  Norm[{d u}] ]; In Omega; Jacobian JVol; } } }
      { Name epsilonr   ; Value { Local { [ CompZZ[epsilonr[] ]   ] ; In Omega; Jacobian JVol; } } }
      }
    }
}


// #############################################################################

PostOperation {
    /*------------ Diffraction problem -----------------*/
    { Name postop_u_target; NameOfPostProcessing postpro ;
        Operation {
        Print[u, OnElementsOf Omega_target ,  Format NodeTable, File "u_target.txt" ];
        Print [ epsilonr  , OnElementsOf Omega, File "epsilonr.pos" ];
        /*Print [ u  , OnElementsOf Omega, File "u.pos" ];*/
      }
      }
      { Name postop_check_adj_source; NameOfPostProcessing postpro ;
          Operation {
              Print [ check_adj_source  , OnElementsOf Omega, File "check_adj_source.pos" ];
        }
    }
    { Name postop_obj_target; NameOfPostProcessing postpro ;
        Operation {
        Print[objective, OnElementsOf Omega_target ,  Format NodeTable, File "obj_target.txt" ];
      }
    }
    { Name postop_adjoint; NameOfPostProcessing postpro ;
        Operation {
        If (nodes_flag)
          Print[u_adj, OnElementsOf Omega_design ,  Format NodeTable, File "adjoint.txt" ];
        Else
          Print[u_adj, OnElementsOf Omega_design , Depth 0, Format SimpleTable, File "adjoint.txt" ];
        EndIf
      }
    }
    { Name postop_source_adj; NameOfPostProcessing postpro ;
        Operation {
          Print[sadj_int_re, OnElementsOf Omega_target, StoreInField  2];
          Print[sadj_int_im, OnElementsOf Omega_target, StoreInField  3];
      }
    }




    { Name postop_dEq_deps; NameOfPostProcessing postpro ;
        Operation {
        If (nodes_flag)
          Print[dEq_deps, OnElementsOf Omega_design ,  Format NodeTable, File "dEq_deps.txt" ];
        Else
          Print[dEq_deps, OnElementsOf Omega_design ,  Depth 0, Format SimpleTable, File "dEq_deps.txt" ];
        EndIf
      }
    }

    { Name postop_int_objective; NameOfPostProcessing postpro ;
       Operation {
         /*Print [ u  , OnElementsOf Omega, File "u.pos" ];*/
       Print[ int_objective[Omega_target],  OnElementsOf PrintPoint, File "objective.txt" , Format SimpleTable ];
       }
    }

    { Name postop_absorption; NameOfPostProcessing postpro ;
       Operation {
       Print[ Q[Omega_source],  OnElementsOf PrintPoint, File "Q.txt" , Format SimpleTable ];
       }
    }


    { Name postop_fields_cuts; NameOfPostProcessing postpro ;
      Operation {
      For i In {0:nb_slice-1}
  	    Print[u_diff , OnLine {{-d/2,ycut_sup~{i},0}{d/2, ycut_sup~{i}, 0}} {npt_integ-1}, File > "sup_field_cuts.out",	Format SimpleTable];
  	    Print[u_tot  , OnLine {{-d/2,ycut_sub~{i},0}{d/2, ycut_sub~{i}, 0}} {npt_integ-1}, File > "sub_field_cuts.out" ,	Format SimpleTable];
      EndFor
      }
    }

    { Name postop_fields_pos; NameOfPostProcessing postpro ;
        Operation {
            Print [ u   , OnElementsOf Omega, File "u.pos" ];
            Print [ epsilonr   , OnElementsOf Omega, File "epsilonr.pos" ];
            Print [ u_tot   , OnElementsOf Omega, File "u_tot.pos" ];
            /* Print [ u_diff   , OnElementsOf Omega, File "u_diff.pos" ];
            Print [ vx_diff   , OnElementsOf Omega, File "vx_diff.pos" ];
            Print [ vy_diff   , OnElementsOf Omega, File "vy_diff.pos" ];
            Print [ vx_tot   , OnElementsOf Omega, File "vx_tot.pos" ];
            Print [ vy_tot   , OnElementsOf Omega, File "vy_tot.pos" ];
            Print [ v_tot   , OnElementsOf Omega, File "v_tot.pos" ];
            Print [ abso_density   , OnElementsOf Omega, File "abso_density.pos" ]; */

        }
    }
    { Name postop_fields_txt; NameOfPostProcessing postpro ;
        Operation {
        Print [ u , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
        { Niy-1, Nix-1} ,Format SimpleTable, File "u.txt" ];
        Print [ u_diff , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
        { Niy-1, Nix-1} ,Format SimpleTable, File "u_diff.txt" ];
        Print [ vx_diff , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
        { Niy-1, Nix-1} ,Format SimpleTable, File "vx_diff.txt" ];
        Print [ vy_diff , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
        { Niy-1, Nix-1} ,Format SimpleTable, File "vy_diff.txt" ];
        Print [ u_tot , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
        { Niy-1, Nix-1} ,Format SimpleTable, File "u_tot.txt" ];
        Print [ vx_tot , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
        { Niy-1, Nix-1} ,Format SimpleTable, File "vx_tot.txt" ];
        Print [ vy_tot , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
        { Niy-1, Nix-1} ,Format SimpleTable, File "vy_tot.txt" ];
        Print [ abso_density , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
        { Niy-1, Nix-1} ,Format SimpleTable, File "abso_density.txt" ];
      }
    }


    /* ----------- Modal analysis --------------------*/
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



    /* ----------- electrostat analysis --------------------*/
    /* { Name postop_fields_electrostat_pos; NameOfPostProcessing postpro_electrostat ;
        Operation {
            Print [ source_electrostat   , OnElementsOf Omega, File "source_electrostat.pos" ];
            Print [ u   , OnElementsOf Omega_source, File "potential.pos" ];
            Print [ epsilonr   , OnElementsOf Omega_source, File "epsilonr.pos" ];
            Print [ e   , OnElementsOf Omega_source, File "static_electric_field.pos" ];
            Print [ norm_e   , OnElementsOf Omega_source, File "norm_static_electric_field.pos" ];

        }
    }
    { Name postop_fields_electrostat_txt; NameOfPostProcessing postpro_electrostat ;
        Operation {
          Print [ u , OnPlane    { { domX_L, domY_B, 0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
          { Niy-1, Nix-1} ,Format SimpleTable, File "potential.txt" ];
          Print [ e , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
          { Niy-1, Nix-1} ,Format SimpleTable, File "static_electric_field.txt" ];
          Print [ norm_e , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
          { Niy-1, Nix-1} ,Format SimpleTable, File "norm_static_electric_field.txt" ];
        }
    } */




}
