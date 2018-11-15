
Include "parameters.dat";

// #############################################################################
// #############################################################################

Group {
    // Domains
    pmlC           = Region[1];
  	pmlTB          = Region[2];
  	pmlLR          = Region[3];
  	host            = Region[{4,7}];
  	design            = Region[5];

    If (target_flag)
      target_circle = Region[7];
    EndIf
    If (inclusion_flag)
      incl          = Region[6];
      Omega_source    = Region[{incl,design}];
      Omega_i = Region[{host,design,incl}];
    Else
      Omega_source    = Region[{design}];
      Omega_i = Region[{host,design}];
    EndIf
    Omega_pml = Region[{pmlC,pmlTB,pmlLR}];
    Omega_nosource  = Region[{host,pmlC,pmlTB,pmlLR}];
    Omega           = Region[{Omega_source,Omega_nosource}];

    Box_B = Region[100];
    Box_R = Region[200];
    Box_T = Region[300];
    Box_L = Region[400];

    // Points
    PrintPoint	=  Region[10000];
    TargetPoint	=  Region[10000];

    If (target_flag)
      Omega_target    = Region[{target_circle}];
    Else
      Omega_target    = Region[{host}];
    EndIf
    /*Omega_target    = Region[{layer1}];*/
    Omega_design    = Region[{design}];



}
// #############################################################################

Function{

    j[] = Complex[0.0, 1.0];
    Freq   = cel/lambda0;
    omega0 = 2.*Pi*cel/lambda0;
    k0     = 2.*Pi/lambda0;
    ksearch = 2.*Pi/lambda0search;
    alpha0 = -k0*Sin[theta];
    beta0  =   k0*Cos[theta];
    Pinc   =  0.5*A*A*Sqrt[epsilon0/mu0] * Cos[theta];
    PW[] = $1* Exp[ j[]* ( $2 * X[] + $3 * Y[] ) ] ;  // plane wave PW[amplitude, kx, ky]


    // PML parameters

    sx[pmlC]         =   Complex[a_pml,-b_pml];
    sx[pmlTB]        =   1.;
    sx[pmlLR]        =   Complex[a_pml,-b_pml];

    sy[pmlC]  	  =   Complex[a_pml,-b_pml];
    sy[pmlTB] 	  =   Complex[a_pml,-b_pml];
    sy[pmlLR] 	  =   1.;
    sz[Omega_pml]      =   1.;



    // Permittivities
    eps[] = -3 -0.01*j[];
    mu[] = -3-0.01*j[];


    If (inclusion_flag)
      epsilonr[design]         = Complex[eps_des_re,eps_des_im] * TensorDiag[1,1,1];
      epsilonr[incl]           = Complex[eps_incl_re,eps_incl_im] * TensorDiag[1,1,1];
      epsilonr_annex[incl]   = Complex[eps_host_re,eps_host_im] * TensorDiag[1,1,1];
    Else
      epsilonr[design]         = Complex[ScalarField[XYZ[], 0, 1 ]{0}, ScalarField[XYZ[], 0, 1 ]{1}] * TensorDiag[1,1,1];
      /* epsilonr[design]         = eps [] *TensorDiag[1,1,1]; */
      /* epsilonr[design]         = TensorDiag[-1.1,-1.1, 1]; */
    EndIf
    epsilonr[host]            = Complex[eps_host_re, eps_host_im] * TensorDiag[1,1,1];
    epsilonr[Omega_pml]      = eps_host_re*TensorDiag[sz[]*sy[]/sx[],sx[]*sz[]/sy[],sx[]*sy[]/sz[]];

    epsilonr_annex[design]    =  Complex[eps_host_re, eps_host_im] * TensorDiag[1,1,1];
    epsilonr_annex[host]      = Complex[eps_host_re, eps_host_im] * TensorDiag[1,1,1];
    epsilonr_annex[Omega_pml]  = eps_host_re*TensorDiag[sz[]*sy[]/sx[],sx[]*sz[]/sy[],sx[]*sy[]/sz[]];


    // Permeabilities
    mur[Omega_pml]              = TensorDiag[sz[]*sy[]/sx[],sx[]*sz[]/sy[],sx[]*sy[]/sz[]];
    mur[host]                 = TensorDiag[1,1,1];
    mur[design]              = TensorDiag[1,1,1];
    /* mur[design]              = mu[]*TensorDiag[1,1,1]; */
    mur_annex[Omega]         = TensorDiag[1,1,1];
    If (inclusion_flag)
      mur[incl]                 = TensorDiag[1,1,1];
    EndIf

    // Fields

    If (ls_flag)
      hankel2[] = Jn[$1,$2] -j[]*Yn[$1,$2];
      Rho[] = Sqrt[(X[]-xs)^2 + (Y[]-ys)^2];
      GF[] = -j[]/4 * hankel2[0, k0*Rho[]];
      Rho_tar[] = Sqrt[(x_target-xs)^2 + (y_target-ys)^2];
      GF_tar[] = -j[]/4 * hankel2[0, k0*Rho_tar[]];
      u_i[Omega_i]=GF[];
      A_beam[] = 1;
      grad_u_i[] = -k0 * j[]/8 * (hankel2[-1, k0*Rho[]] - hankel2[1, k0*Rho[]]) / Rho[] * Vector[X[]-xs, Y[]-ys, 0];
    Else
      If (beam_flag)
        /* Xrot[] = $X * Sin[theta] + $Y * Cos[theta]; */
        Yrot[] = $X * Cos[theta] + $Y * Sin[theta];
        A_beam[] = Exp[- (Yrot[]^2)/(2 * waist^2) ];
        pw[] =  PW[A, alpha0, beta0];
        u_i[Omega_i] = A_beam[] *pw[];
        gradpw[] = j[]*pw[] * TensorDiag[alpha0, beta0,0.];
        gradAbeam[] = -1/waist^2 * Yrot[] * A_beam[]* TensorDiag[Cos[theta], Sin[theta],0.];
        grad_u_i[]  = A_beam[] * gradpw[] + gradAbeam[] * pw[];
      Else
        u_i[Omega_i] = PW[A, alpha0, beta0];
        grad_u_i[] = j[]*u_i[] * Vector[alpha0, beta0,0.];
        A_beam[] = 1;
      EndIf

    EndIf

    u_i[Omega_pml] = 0.;

    If (ls_flag)
      grad_A_beam[Omega_i] = 0 * TensorDiag[1, 1, 0.] ;
    Else
      If (beam_flag)
        grad_A_beam[Omega_i] = 0* TensorDiag[1, 1, 0.] ;
      Else
        grad_A_beam[Omega_i] =  0* TensorDiag[1, 1, 0.];
      EndIf

    EndIf
    grad_A_beam[Omega_pml] =  0* TensorDiag[1, 1, 0.];


    /*FIXME: for the beam case there is an extra term in source
    coming from the gradient of spatially varying amplitude  */

    /* TE: */
    /* (Hx, Hy) = (dxEz, -dxEz )/(j*omega*mu0*mur) */


    If (TE_flag)

        If (ls_flag)
          dual_i[]       = - Vector[Y[] - ys, -X[] + xs, 0] * k0* ( hankel2[-1, k0*Rho[]] - hankel2[1, k0*Rho[]])/(  8 * Rho[] ) * 1 / (omega0*mu0*CompXX[mur_annex[]]) ;
        Else
          dual_i[]       = Vector[beta0, -alpha0, 0] * u_i[] / (omega0*mu0*CompXX[mur_annex[]])  ;
        EndIf

        source_eps[]   =  k0^2*(epsilonr[]-epsilonr_annex[])*u_i[];
        source_mu[]    =  (1/mur_annex[]-1/mur[])  * grad_u_i[] ;
        weight[] = epsilonr[];
    Else
        dual_i[]        = Vector[ beta0, -alpha0, 0] *u_i[] / (omega0*epsilon0*CompXX[epsilonr_annex[]]);
        source_eps[]    =  (1/epsilonr_annex[]-1/epsilonr[])* grad_u_i[];
        source_mu[]     =  k0^2*(mur[]-mur_annex[])*u_i[];
        weight[] = mur[];

    EndIf

    source[]  =  source_eps[] + source_mu[] ;




    // angular loop
    it=0;
    For th In {0:2*Pi:2*Pi/(Ni_theta-1)}
        alpha0_t = -k0*Sin[th];
        beta0_t  =   k0*Cos[th];
        u_i~{it}[]= PW[A, alpha0_t, beta0_t];

          source~{it}[]  =  (epsilonr[]-epsilonr_annex[])*u_i~{it}[];


        it = it+1;
    EndFor

    modR[] = Sqrt[X[]^2 + Y[]^2];
    phiR[] = Atan2[-X[], Y[]];

    im=0;
    For m In {-M_fs:M_fs:1}
        source_fs~{im}[]  = j[]^m*(epsilonr[]-epsilonr_annex[])*Jn[m, k0*modR[]] * Exp[-j[]*m*phiR[]];
        im = im+1;
    EndFor


/*
    For i In {0:nb_slice-1}
      ycut_sub~{i} = ycut_sub_min + i*(ycut_sub_max-ycut_sub_min)/(nb_slice-1);
      ycut_sup~{i} = ycut_sup_min + i*(ycut_sup_max-ycut_sup_min)/(nb_slice-1);
    EndFor*/

    coef_Q[] = 0.5 * epsilon0*omega0*Fabs[Im[CompZZ[epsilonr[]]]   ]  ;
    /* coef_obj[] = coef_Q[];  // 1/(h_sub*d); */


    If (TE_flag)
          dual[] = Vector[ CompY[ $1 ], - CompX[ $1 ], 0 ]/(j[]*omega0*mu0*CompXX[mur[]]); // TE case H = dual[{d u}]
          dual_tot[] = dual[$1] + dual_i[] * CompXX[mur_annex[]] /CompXX[mur[]] ;
          absorption[]  =  coef_Q[] * SquNorm[$1 + u_i[]] ;
    Else
          dual[] = Vector[ CompY[ $1 ], - CompX[ $1 ], 0 ]/(-j[]*omega0*epsilon0*CompXX[epsilonr[]]); // TM case E = dual[{d u}]
          dual_tot[] = -dual[$1] + dual_i[] * CompXX[epsilonr_annex[]] /CompXX[epsilonr[]] ;
          absorption[]  =  coef_Q[] * ( SquNorm[  CompX[dual_tot[$1]] ] + SquNorm[CompY[dual_tot[$1]] ] );

    EndIf


    // near to far field
    If (ls_flag)
      n2f_field[]= $1 + u_i[] ;
      n2f_field_dual[] = dual_tot[$1];
      ui_tar[] = GF_tar[];
    Else
      n2f_field[]= $1;
      n2f_field_dual[] = dual[$1];
      ui_tar[] = 1;
      /* n2f_field[]= $1 + u_i[] ;
      n2f_field_dual[] = dual_tot[$1]; */
    EndIf

    // Topology optimization
    coef_obj[] = 1/(SquNorm[ui_tar[]]*Pi*r_target^2);
    objective[] = coef_obj[] *SquNorm[$1 + u_i[]] ;
    adj_source_int[] =  -2 * coef_obj[] * (Conj[ $1 + u_i[]]); //d_objective_du *ElementVol[]

    If (TE_flag)
      db_deps[] = -k0^2*u_i[];
      dA_deps[] = k0^2;
    Else
      db_deps[] = 1/CompXX[epsilonr[]]^2* grad_u_i[];
      dA_deps[] = -1/CompXX[epsilonr[]]^2;

    EndIf
    dEq_deps[] = db_deps[] - dA_deps[] * ($1);


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
   }
  }
 { Name Hgraduadj; Type Form0;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_Node;    Support Region[Omega]; Entity NodesOf[Omega]; }
      { Name sn2; NameOfCoef un2; Function BF_Node_2E; Support Region[Omega]; Entity EdgesOf[Omega]; }
     }
    Constraint {
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
            Galerkin { [-1/TensorDiag[CompYY[mur[]],CompXX[mur[]],CompZZ[mur[]]]*Dof{d u} , {d u}];
            In Omega; Jacobian JVol; Integration Int_1; }
            Galerkin { [ ($Source ? source_eps[] : 0) , {u}];
            In Omega_source; Jacobian JVol; Integration Int_1;  }
            Galerkin { [ ($Source ? source_mu[] : 0) , {d u}];
            In Omega_source; Jacobian JVol; Integration Int_1;  }
            /* Galerkin { [ ($SourceAdj ? adj_source[{u}] : 0) , {u}];
            In Omega_target; Jacobian JVol; Integration Int_1;  } */
            Galerkin { [ ($SourceAdj ? Complex[ScalarField[XYZ[], 0, 1 ]{2}, ScalarField[XYZ[], 0, 1 ]{3}] : 0) , {u}];
            In Omega_target; Jacobian JVol; Integration Int_1;  }
        Else
            Galerkin { [k0^2*CompZZ[mur[]]*Dof{u} , {u}];
            In Omega; Jacobian JVol; Integration Int_1;  }
            Galerkin { [-1/TensorDiag[CompYY[epsilonr[]],CompXX[epsilonr[]],CompZZ[epsilonr[]]]*Dof{d u} , {d u}];
            In Omega; Jacobian JVol; Integration Int_1; }
            Galerkin { [ ($Source ? source_eps[] : 0) , {d u}];
            In Omega_source; Jacobian JVol; Integration Int_1;  }
            Galerkin { [ ($Source ? source_mu[] : 0) , {u}];
            In Omega_source; Jacobian JVol; Integration Int_1;  }
            Galerkin { [ ($SourceAdj ? Complex[ScalarField[XYZ[], 0, 1 ]{2}, ScalarField[XYZ[], 0, 1 ]{3}] : 0) , {u}];
            In Omega_target; Jacobian JVol; Integration Int_1;  }
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

        Generate[S] ;Solve[S] ;SaveSolution[S] ;

        If (adjoint_flag)
          PostOperation[postop_int_objective];
          PostOperation[postop_dEq_deps];
          PostOperation[postop_source_adj];
          Evaluate[$Source = 0, $SourceAdj = 1];
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

}



// #############################################################################

PostProcessing {
    /*------------ Diffraction problem -----------------*/
    { Name postpro; NameOfFormulation helmoltz_scalar;
        Quantity {
          { Name check_adj_source   ; Value { Local { [adj_source_int[{u}]] ; In Omega; Jacobian JVol; } } }
            { Name objective   ; Value { Local { [objective[{u}]] ; In Omega; Jacobian JVol; } } }
            { Name int_objective  ; Value { Integral { [objective[{u}]] ; In Omega_target    ; Integration Int_1 ; Jacobian JVol ; } } }
            { Name u   ; Value { Local { [{u}  ] ; In Omega; Jacobian JVol; } } }
            { Name u_i   ; Value { Local { [ u_i[]    ] ; In Omega; Jacobian JVol; } } }
            { Name epsilonr   ; Value { Local { [ CompZZ[epsilonr[] ]   ] ; In Omega; Jacobian JVol; } } }
            { Name epsilonr_annex   ; Value { Local { [ CompZZ[epsilonr_annex[] ]   ] ; In Omega; Jacobian JVol; } } }
            { Name u_diff   ; Value { Local { [ {u}     ] ; In Omega; Jacobian JVol; } } }
            { Name u_tot    ; Value { Local { [ {u}+u_i[]       ] ; In Omega; Jacobian JVol; } } }
            { Name u_sqnorm    ; Value { Local { [ ({u})*Conj[({u}) ]       ] ; In Omega; Jacobian JVol; } } }
            { Name u_tot_norm    ; Value { Local { [Sqrt[ ({u}+u_i[])*Conj[({u}+u_i[]) ] ]       ] ; In Omega; Jacobian JVol; } } }
            { Name source    ; Value { Local { [CompZZ[source[]]]     ; In Omega; Jacobian JVol; } } }
            { Name u_tot_sqnorm    ; Value { Local { [ ({u}+u_i[])*Conj[({u}+u_i[]) ]       ] ; In Omega; Jacobian JVol; } } }
            { Name vx_diff   ; Value { Local { [ CompX[dual[{d u}] ]      ] ; In Omega; Jacobian JVol; } } }
            { Name vy_diff    ; Value { Local { [  CompY[ dual[{d u}]  ]    ] ; In Omega; Jacobian JVol; } } }
            { Name vx_tot    ; Value { Local { [  CompX[dual_i[] + dual[{d u}]  ]         ] ; In Omega; Jacobian JVol; } } }
            { Name vy_tot   ; Value { Local { [  CompY[dual_i[]   + dual[{d u}]  ]  ] ; In Omega; Jacobian JVol; } } }
            { Name v_tot   ; Value { Local { [  dual_i[]   + dual[{d u}]   ] ; In Omega; Jacobian JVol; } } }
            { Name u_int    ; Value { Integral { [ {u}      ] ; In Omega; Integration Int_1; Jacobian JVol; } } }
            { Name u_elvol   ; Value { Local { [ {u} * ElementVol[]   ] ; In Omega; Jacobian JVol; } } }
            { Name sadj_int_re  ; Value { Integral { [Re[ adj_source_int[{u}] ]] ; In Omega_target    ; Integration Int_1 ; Jacobian JVol ; } } }
            { Name sadj_int_im  ; Value { Integral { [Im[ adj_source_int[{u}] ]] ; In Omega_target    ; Integration Int_1 ; Jacobian JVol ; } } }
            { Name u_adj   ; Value { Local { [ {u} * ElementVol[]   ] ; In Omega; Jacobian JVol; } } }
            If (TE_flag)
              { Name dEq_deps   ; Value { Local { [dEq_deps[{u}]] ; In Omega; Jacobian JVol; } } }
              { Name Q  ; Value { Integral { [ absorption[{u}] ] ; In Omega_source    ; Integration Int_1 ; Jacobian JVol ; } } }
                { Name abso_density   ; Value { Local { [ absorption[{u}] ]; In Omega_source; Jacobian JVol; } } }
            Else
              /* { Name dEq_deps   ; Value { Local { [dEq_deps[{d u}]] ; In Omega; Jacobian JVol; } } } */
              { Name dEq_deps_x   ; Value { Local { [CompX[dEq_deps[{d u}]]] ; In Omega; Jacobian JVol; } } }
              { Name dEq_deps_y   ; Value { Local { [CompY[dEq_deps[{d u}]]] ; In Omega; Jacobian JVol; } } }
              { Name Q ; Value { Integral { [ absorption[{d u}] ] ; In Omega_source ; Integration Int_1 ; Jacobian JVol ; } } }
              { Name abso_density   ; Value { Local {[ absorption[{d u}] ]; In Omega_source; Jacobian JVol; } } }
            EndIf
            { Name n2f_field   ; Value { Local { [n2f_field[{u}]  ] ; In Omega; Jacobian JVol; } } }
            { Name n2f_field_dual_x   ; Value { Local { [CompX[n2f_field_dual[{d u}]] ] ; In Omega; Jacobian JVol; } } }
            { Name n2f_field_dual_y   ; Value { Local { [CompY[n2f_field_dual[{d u}]] ] ; In Omega; Jacobian JVol; } } }

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

          For it In {0:Ni_theta-1}
              { Name CouplingCoeffs~{it};  Value { Integral { [  CompZZ[source~{it}[]]*{u}     ] ;
               In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
          EndFor
          For im In {0:2*M_fs}
            { Name CouplingCoeffsFS~{im};  Value { Integral { [  CompZZ[source_fs~{im}[]]*{u}     ] ;
             In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }
          EndFor
            { Name mode_coupling;  Value { Integral { [  CompZZ[source[]]*{u}     ] ;
             In Omega    ; Integration Int_1 ; Jacobian JVol ; } } }

             { Name mode_coupling_int;  Value { Local { [  CompZZ[source[]]*{u}     ] ;
              In Omega    ; Jacobian JVol ; } } }

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
        If (TE_flag)
          If (nodes_flag)
            Print[dEq_deps, OnElementsOf Omega_design ,  Format NodeTable, File "dEq_deps.txt" ];
          Else
            Print[dEq_deps, OnElementsOf Omega_design ,   Depth 0, Format SimpleTable, File "dEq_deps.txt" ];
          EndIf
        Else
          If (nodes_flag)
            Print[dEq_deps_x, OnElementsOf Omega_design ,  Format NodeTable, File "dEq_deps_x.txt" ];
            Print[dEq_deps_y, OnElementsOf Omega_design ,  Format NodeTable, File "dEq_deps_y.txt" ];
          Else
            Print[dEq_deps_x, OnElementsOf Omega_design ,   Depth 0, Format SimpleTable, File "dEq_deps_x.txt" ];
            Print[dEq_deps_y, OnElementsOf Omega_design ,   Depth 0, Format SimpleTable, File "dEq_deps_y.txt" ];
          EndIf
          /* Print[dEq_deps_x, OnElementsOf Omega_design ,   File "dEq_deps_x.pos" ];
          Print[dEq_deps_y, OnElementsOf Omega_design ,  File "dEq_deps_y.pos" ]; */

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


    { Name postop_fields_box; NameOfPostProcessing postpro ;
      Operation {
  	    Print[u , OnLine {{-hx_des/2,-hy_des/2,0}{hx_des/2, -hy_des/2, 0}} {Nibox_x-1}, File  "field_box_B.out",	Format TimeTable];
  	    Print[u , OnLine {{-hx_des/2,hy_des/2,0}{hx_des/2, hy_des/2, 0}} {Nibox_x-1}, File  "field_box_T.out",	Format TimeTable];
        Print[u , OnLine {{-hx_des/2,-hy_des/2,0}{-hx_des/2, hy_des/2, 0}} {Nibox_y-1}, File  "field_box_L.out",	Format TimeTable];
        Print[u , OnLine {{hx_des/2,-hy_des/2,0}{hx_des/2, hy_des/2, 0}} {Nibox_y-1}, File  "field_box_R.out",	Format TimeTable];
      }
    }

    { Name postop_fields_n2f; NameOfPostProcessing postpro ;
      Operation {
        Print[n2f_field , OnLine {{Xn2f_L,Yn2f_B,0}{Xn2f_R, Yn2f_B, 0}} {Nin2f_x-1}, File  "field_n2f_B.out",	Format TimeTable];
        Print[n2f_field , OnLine {{Xn2f_L,Yn2f_T,0}{Xn2f_R, Yn2f_T, 0}} {Nin2f_x-1}, File  "field_n2f_T.out",	Format TimeTable];
        Print[n2f_field , OnLine {{Xn2f_L,Yn2f_B,0}{Xn2f_L, Yn2f_T, 0}} {Nin2f_y-1}, File  "field_n2f_L.out",	Format TimeTable];
        Print[n2f_field , OnLine {{Xn2f_R,Yn2f_B,0}{Xn2f_R, Yn2f_T, 0}} {Nin2f_y-1}, File  "field_n2f_R.out",	Format TimeTable];
        Print[n2f_field_dual_x , OnLine {{Xn2f_L,Yn2f_B,0}{Xn2f_R, Yn2f_B, 0}} {Nin2f_x-1}, File  "field_dual_x_n2f_B.out",	Format TimeTable];
        Print[n2f_field_dual_x , OnLine {{Xn2f_L,Yn2f_T,0}{Xn2f_R, Yn2f_T, 0}} {Nin2f_x-1}, File  "field_dual_x_n2f_T.out",	Format TimeTable];
        Print[n2f_field_dual_x , OnLine {{Xn2f_L,Yn2f_B,0}{Xn2f_L, Yn2f_T, 0}} {Nin2f_y-1}, File  "field_dual_x_n2f_L.out",	Format TimeTable];
        Print[n2f_field_dual_x , OnLine {{Xn2f_R,Yn2f_B,0}{Xn2f_R, Yn2f_T, 0}} {Nin2f_y-1}, File  "field_dual_x_n2f_R.out",	Format TimeTable];
        Print[n2f_field_dual_y , OnLine {{Xn2f_L,Yn2f_B,0}{Xn2f_R, Yn2f_B, 0}} {Nin2f_x-1}, File  "field_dual_y_n2f_B.out",	Format TimeTable];
        Print[n2f_field_dual_y , OnLine {{Xn2f_L,Yn2f_T,0}{Xn2f_R, Yn2f_T, 0}} {Nin2f_x-1}, File  "field_dual_y_n2f_T.out",	Format TimeTable];
        Print[n2f_field_dual_y , OnLine {{Xn2f_L,Yn2f_B,0}{Xn2f_L, Yn2f_T, 0}} {Nin2f_y-1}, File  "field_dual_y_n2f_L.out",	Format TimeTable];
        Print[n2f_field_dual_y , OnLine {{Xn2f_R,Yn2f_B,0}{Xn2f_R, Yn2f_T, 0}} {Nin2f_y-1}, File  "field_dual_y_n2f_R.out",	Format TimeTable];
      }
    }


    /*{ Name postop_fields_cuts; NameOfPostProcessing postpro ;
      Operation {
      For i In {0:nb_slice-1}
  	    Print[u_diff , OnLine {{-d/2,ycut_sup~{i},0}{d/2, ycut_sup~{i}, 0}} {npt_integ-1}, File > "sup_field_cuts.out",	Format SimpleTable];
  	    Print[u_tot  , OnLine {{-d/2,ycut_sub~{i},0}{d/2, ycut_sub~{i}, 0}} {npt_integ-1}, File > "sub_field_cuts.out" ,	Format SimpleTable];
      EndFor
      }
    }*/

    { Name postop_fields_pos; NameOfPostProcessing postpro ;
        Operation {
            Print [ u   , OnElementsOf Omega, File "u.pos" ];
            /* Print [ u_int   , OnElementsOf Omega, File "u_int.pos" ];
            Print [ u_elvol   , OnElementsOf Omega, File "u_elvol.pos" ]; */

            Print [ epsilonr   , OnElementsOf design, File "epsilonr.pos" ];
            /* Print [ epsilonr_annex   , OnElementsOf Omega, File "epsilonr_annex.pos" ]; */
            /* Print [ u_diff   , OnElementsOf Omega, File "u_diff.pos" ]; */
             /*Print [ u_sqnorm   , OnElementsOf Omega_pml, File "u_sqnorm.pos" ]; */
            /* Print [ source   , OnElementsOf Omega, File "source.pos" ]; */

            /*Print [ vx_diff   , OnElementsOf Omega, File "vx_diff.pos" ];*/
            /*Print [ vy_diff   , OnElementsOf Omega, File "vy_diff.pos" ];*/
            Print [ u_tot_norm   , OnElementsOf Omega, File "u_tot_norm.pos" ];

            Print [ u_tot   , OnElementsOf Omega, File "u_tot.pos" ];
            Print [ u_i   , OnElementsOf Omega, File "u_i.pos" ];
            /*Print [ vx_tot   , OnElementsOf Omega, File "vx_tot.pos" ];*/
            /*Print [ vy_tot   , OnElementsOf Omega, File "vy_tot.pos" ];*/
            /*Print [ v_tot   , OnElementsOf Omega, File "v_tot.pos" ];*/
            /*Print [ abso_density   , OnElementsOf Omega, File "abso_density.pos" ];*/

        }
    }
    { Name postop_fields_txt; NameOfPostProcessing postpro ;
        Operation {
        /* Print [ u , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
        { Niy-1, Nix-1} ,Format SimpleTable, File "u.txt" ];
        Print [ u_diff , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
        { Niy-1, Nix-1} ,Format SimpleTable, File "u_diff.txt" ];
        Print [ vx_diff , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
        { Niy-1, Nix-1} ,Format SimpleTable, File "vx_diff.txt" ];
        Print [ vy_diff , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
        { Niy-1, Nix-1} ,Format SimpleTable, File "vy_diff.txt" ]; */
        Print [ u_tot , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
        { Niy-1, Nix-1} ,Format SimpleTable, File "u_tot.txt" ];
        /* Print [ vx_tot , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
        { Niy-1, Nix-1} ,Format SimpleTable, File "vx_tot.txt" ];
        Print [ vy_tot , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
        { Niy-1, Nix-1} ,Format SimpleTable, File "vy_tot.txt" ];
        Print [ abso_density , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
        { Niy-1, Nix-1} ,Format SimpleTable, File "abso_density.txt" ]; */
      }
    }


    { Name postop_field_on_point; NameOfPostProcessing postpro ;
        Operation {
            Print [u_tot, OnPoint {xpp, ypp, 0}, Format SimpleTable, File "u_tot_point.txt"];
            Print [u_i, OnPoint {xpp, ypp, 0}, Format SimpleTable, File "u_i_point.txt"];
            Print [u, OnPoint {xpp, ypp, 0}, Format SimpleTable, File "u_point.txt"];

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

    { Name postop_coupling_coeffs_angle; NameOfPostProcessing postpro_modal ;
        Operation {
          For it In {0:Ni_theta-1}
            Print [CouplingCoeffs~{it}[Omega], OnElementsOf PrintPoint, Format TimeTable, File >"coupling_coeffs.txt"];
          EndFor
          }
    }
    { Name postop_coupling_coeffs_fourrier_series; NameOfPostProcessing postpro_modal ;
        Operation {
          For im In {0:2*M_fs}
            Print [CouplingCoeffsFS~{im}[Omega], OnElementsOf PrintPoint, Format TimeTable, File >"coupling_coeffs_fs.txt"];
          EndFor
          }
    }


  /* ----------- QNM expansion --------------------*/
  { Name postop_mode_coupling; NameOfPostProcessing postpro_modal ;
      Operation {
          Print [mode_coupling[Omega], OnElementsOf PrintPoint, Format TimeTable, File "mode_coupling.txt"];
        }
  }

  { Name postop_mode_coupling_int; NameOfPostProcessing postpro_modal ;
      Operation {
      Print [ mode_coupling_int , OnPlane    { { domX_L,domY_B,0 } { domX_L,domY_T,0 } { domX_R,domY_B,0 } }
      { Niy-1, Nix-1} , Format TimeTable, File "mode_coupling_int.txt" ];
      }
  }



}




// #############################################################################
