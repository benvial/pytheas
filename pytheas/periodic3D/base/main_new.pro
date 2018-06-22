Include "parameters.dat";

Group {

	// SubDomains
	If (incl_flag)
		nr = 4;
		rod_5 = Region[8];
		rod_4= Region[7];
		rod_3 = Region[6];
		rod_2 = Region[5];
		rods = Region[{rod_2, rod_3, rod_4, rod_5}];
	Else
		nr = 0;
	EndIf
	L_5            = Region[1];
	L_4            = Region[2];
	L_3            = Region[3];  // interp
	L_2            = Region[4];
	/* index =  {(5+nr):(13+nr)}; */
	/* For j In {0:8}
	  Printf("Volume number = %g", index(j));
	EndFor */
	/* PMLbot         = Region[{index(3)}];
	L_6            = Region[{index(1)}];
	L_1            = Region[{index(0)}];
	PMLtop         = Region[{index(2)}];
	SurfBlochXm    = Region[{index(4)}];
	SurfBlochXp    = Region[{index(5)}];
	SurfBlochYm    = Region[{index(6)}];
	SurfBlochYp    = Region[{index(7)}];
	SurfDirichlet  = Region[{index(8)}];
	 */
	PMLbot         = Region[8];
	L_6            = Region[6];
	L_1            = Region[5];
	PMLtop         = Region[7];
	SurfBlochXm    = Region[9];
	SurfBlochXp    = Region[10];
	SurfBlochYm    = Region[11];
	SurfBlochYp    = Region[12];
	SurfDirichlet  = Region[13];

	// Boundaries

	SurfBloch      = Region[{SurfBlochXm,SurfBlochXp,SurfBlochYm,SurfBlochYp}];
	// Lines
	/* LineNeumann    = Region[10001]; */
	// Domains
	Omega_interp   = Region[{L_3}];
	Omega_nosource = Region[{PMLbot,L_6,L_1,PMLtop}];
	Omega_subs     = Region[{L_6,PMLbot}];
	Omega_source   = Region[{L_2,L_3,L_4,L_5}];
	Omega_super    = Region[{L_1,L_2,L_3,L_4,L_5,PMLtop}];
	If (incl_flag)
		Omega_super    = Region[{Omega_super, rods}];
		Omega_source   = Region[{Omega_source, rods}];
	EndIf
	Omega_nopml = Region[{Omega_source, L_1, L_6}];
	
	Omega          = Region[{Omega_source, Omega_nosource}];

	// SurfNeumann    = Region[{SurfBlochXm,SurfBlochXp,SurfBlochYm,SurfBlochYp}];
}



Function{
	Freq             = cel/lambda0;
	omega0           = 2*Pi*Freq;
	k0               = 2.*Pi/lambda0;
	Ae               = 1;
	Ah               = Ae*Sqrt[epsilon0/mu0];
	alpha0           = k0*Sin[theta0]*Cos[phi0];
	beta0            = k0*Sin[theta0]*Sin[phi0];
	gamma0           = k0*Cos[theta0];
	Pinc             =  0.5*Ae*Ae*Sqrt[epsilon0/mu0] * Cos[theta0];
	For i In {1:6}
		epsilon[L~{i}]  = Complex[eps_re_L~{i} , eps_im_L~{i}] * TensorDiag[1,1,1];
	EndFor
	/* For i In {4:6}
		epsilon[L~{i}]  = Complex[eps_re_L~{i} , eps_im_L~{i}] * TensorDiag[1,1,1];
	EndFor
	 */
	If (incl_flag)
		epsilon[L_3]  = Complex[eps_re_L_3 , eps_im_L_3] * TensorDiag[1,1,1];
		
		For i In {2:5}
			eps_re_rod~{i} = 4.;
			eps_im_rod~{i} = -0.0;
			epsilon[rod~{i}]    = Complex[eps_re_rod~{i} , eps_im_rod~{i}] * TensorDiag[1,1,1];
			
		EndFor
		/* i=5;
		eps_re_rod~{i} = 1.;
		eps_im_rod~{i} = -0.;
		epsilon[rod~{i}]    = Complex[eps_re_rod~{i} , eps_im_rod~{i}] * TensorDiag[1,1,1];
		 */
 	Else
		/* epsilon[Omega_interp]    = Complex[ScalarField[XYZ[], 0, 1]{0}  ,ScalarField[XYZ[], 0, 1 ]{1} ]  * TensorDiag[1,1,1]; */
		/* epsilon[Omega_interp]          = Complex[eps_re_L_3 , eps_im_L_3] * TensorDiag[1,1,1]; */

	EndIf
	mu[Omega_nopml]       = TensorDiag[1,1,1];

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

	a_pml           = 1.;
	b_pml           = 1.;
	sx              = 1.;
	sz[]            = Complex[a_pml,-b_pml];
	sy              = 1.;

// Why this does not work?
	epsilon[PMLbot] = eps_re_L_6*TensorDiag[sz[]*sy/sx,sx*sz[]/sy,sx*sy/sz[]];
	epsilon[PMLtop] = eps_re_L_1*TensorDiag[sz[]*sy/sx,sx*sz[]/sy,sx*sy/sz[]];
	mu[PMLbot]      =            TensorDiag[sz[]*sy/sx,sx*sz[]/sy,sx*sy/sz[]];
	mu[PMLtop]      =            TensorDiag[sz[]*sy/sx,sx*sz[]/sy,sx*sy/sz[]];

	epsilon_annex[PMLbot]       = eps_re_L_6*TensorDiag[sz[]*sy/sx,sx*sz[]/sy,sx*sy/sz[]];
	epsilon_annex[PMLtop]       = eps_re_L_1*TensorDiag[sz[]*sy/sx,sx*sz[]/sy,sx*sy/sz[]];
	epsilon_annex[Omega_source] = Complex[eps_re_L_1 , eps_im_L_1] * TensorDiag[1,1,1];
	epsilon_annex[L_1]          = Complex[eps_re_L_1 , eps_im_L_1] * TensorDiag[1,1,1];
	epsilon_annex[L_6]          = Complex[eps_re_L_6 , eps_im_L_6] * TensorDiag[1,1,1];

	source[]         = (omega0/cel)^2*(epsilon[]-epsilon_annex[])*Ecm[];//(nm/1.e-9)^2*

	//////////////////////////////////////
  For i In {0:nb_slice-1}
		zcut_sub~{i} = hh_L_6+thick_L_6-scan_dist - i*(thick_L_6-2*scan_dist)/(nb_slice-1);
		zcut_sup~{i} = hh_L_1+scan_dist           + i*(thick_L_1-2*scan_dist)/(nb_slice-1);

	EndFor


}

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

  /*   { NameOfCoef un3; EntityType FacetsOf ; NameOfConstraint BlochX; }
     { NameOfCoef un3; EntityType FacetsOf ; NameOfConstraint BlochY; }
     { NameOfCoef un3; EntityType FacetsOf ; NameOfConstraint Dirichlet; }
     { NameOfCoef un4; EntityType FacetsOf ; NameOfConstraint BlochX; }
     { NameOfCoef un4; EntityType FacetsOf ; NameOfConstraint BlochY; }
     { NameOfCoef un4; EntityType FacetsOf ; NameOfConstraint Dirichlet; }
     { NameOfCoef un5; EntityType EdgesOf ; NameOfConstraint BlochX; }
     { NameOfCoef un5; EntityType EdgesOf ; NameOfConstraint BlochY; }
     { NameOfCoef un5; EntityType EdgesOf ; NameOfConstraint Dirichlet; } */
      }
  }
}


Formulation {{Name helmholtz_vector; Type FemEquation;
    		Quantity {{ Name u; Type Local; NameOfSpace Hcurl;}}
		Equation { Galerkin {[-1/mu[]*Dof{Curl u} , {Curl u}];
                 		In Omega; Jacobian JVol; Integration Int_1;  }
                   	   Galerkin { [(omega0/cel)^2*epsilon[]*Dof{u} , {u}];
                 		In Omega; Jacobian JVol; Integration Int_1;  }
                   	   Galerkin { [source[] , {u}];
                 		In Omega; Jacobian JVol; Integration Int_1;  }
                }
            }
        }

Resolution {
  { Name helmholtz_vector;
    System {
      { Name M; NameOfFormulation helmholtz_vector; Type ComplexValue; Frequency Freq; }
    }
    Operation {
	/*GmshOpen["tmp.pos"];*/
	Generate[M]; Solve[M]; SaveSolution[M];

    }
  }
}



        // Hinc[] : Complex[0,1] * 1/omega0 * 1/mu0 * Curl Einc[];
        // H_d    : Complex[0,1] * 1/omega0 * 1/mu0 * {Curl u};

PostProcessing {
    { Name get_Ed; NameOfFormulation helmholtz_vector;NameOfSystem M;
            Quantity {

		//// E diffracted 3 components, Im and Re parts


              { Name ex_re_t; Value { Local { [Re[  CompX[{u}+Ecm[]] ]]; In Omega; Jacobian JVol; } } }
              { Name ey_re_t; Value { Local { [Re[  CompY[{u}+Ecm[]] ]]; In Omega; Jacobian JVol; } } }
              { Name ez_re_t; Value { Local { [Re[  CompZ[{u}+Ecm[]] ]]; In Omega; Jacobian JVol; } } }
              { Name ex_im_t; Value { Local { [Im[  CompX[{u}+Ecm[]] ]]; In Omega; Jacobian JVol; } } }
              { Name ey_im_t; Value { Local { [Im[  CompY[{u}+Ecm[]] ]]; In Omega; Jacobian JVol; } } }
              { Name ez_im_t; Value { Local { [Im[  CompZ[{u}+Ecm[]] ]]; In Omega; Jacobian JVol; } } }

              { Name ex_re_r; Value { Local { [Re[  CompX[{u}+Em[]] ]]; In Omega; Jacobian JVol; } } }
              { Name ey_re_r; Value { Local { [Re[  CompY[{u}+Em[]] ]]; In Omega; Jacobian JVol; } } }
              { Name ez_re_r; Value { Local { [Re[  CompZ[{u}+Em[]] ]]; In Omega; Jacobian JVol; } } }
              { Name ex_im_r; Value { Local { [Im[  CompX[{u}+Em[]] ]]; In Omega; Jacobian JVol; } } }
              { Name ey_im_r; Value { Local { [Im[  CompY[{u}+Em[]] ]]; In Omega; Jacobian JVol; } } }
              { Name ez_im_r; Value { Local { [Im[  CompZ[{u}+Em[]] ]]; In Omega; Jacobian JVol; } } }


              { Name ex_t; Value { Local { [ CompX[{u}+Ecm[]] ]; In Omega; Jacobian JVol; } } }
              { Name ey_t; Value { Local { [ CompY[{u}+Ecm[]] ]; In Omega; Jacobian JVol; } } }
              { Name ez_t; Value { Local { [ CompZ[{u}+Ecm[]] ]; In Omega; Jacobian JVol; } } }
              { Name ex_r; Value { Local { [ CompX[{u}+Em[]] ]; In Omega; Jacobian JVol; } } }
              { Name ey_r; Value { Local { [ CompY[{u}+Em[]] ]; In Omega; Jacobian JVol; } } }
              { Name ez_r; Value { Local { [ CompZ[{u}+Em[]] ]; In Omega; Jacobian JVol; } } }

              { Name Ecal; Value { Local { [ {u}       ]; In Omega; Jacobian JVol; } } }
              { Name Etot; Value { Local { [ {u}+Ecm[] ]; In Omega; Jacobian JVol; } } }
              { Name Edif; Value { Local { [ {u}+Em[]  ]; In Omega; Jacobian JVol; } } }
              { Name Ecm ; Value { Local { [    Ecm[]  ]; In Omega; Jacobian JVol; } } }
              { Name source; Value { Local { [ source[]]; In Omega; Jacobian JVol; } } }
              { Name Ecal_x; Value { Local { [ CompX[{u}] ]; In Omega; Jacobian JVol; } } }
              { Name Ecal_y; Value { Local { [ CompY[{u}] ]; In Omega; Jacobian JVol; } } }
              { Name Ecal_z; Value { Local { [ CompZ[{u}] ]; In Omega; Jacobian JVol; } } }

              { Name po_r ; Value { Local { [ 0.5*Im[ 1/omega0 * 1/mu0*{u}/\Conj[{Curl u}] ] ]; In Omega; Jacobian JVol; } } }
		    //// Joule Losses Distribution = 1/2*sigma*|E_tot|^2
		  { Name normalized_losses1 ; Value { Integral { [  epsilon0*omega0 * 0.5*Fabs[Im[CompXX[epsilon[]]]]*(SquNorm[{u}+Ecm[]]) / (Pinc*period_x*period_y) ] ; In Omega_source ; Integration Int_1 ; Jacobian JVol ; } } }



		  // VERIFS
		  { Name cooX; Value { Local { [X[]]; In Omega; Jacobian JVol; } } }
		  { Name cooY; Value { Local { [Y[]]; In Omega; Jacobian JVol; } } }
		  { Name cooXpYs2; Value { Local { [X[]+Y[]/2]; In Omega; Jacobian JVol; } } }
							{ Name ex_d_re; Value { Local { [Re[  CompX[{u}] ]]; In Omega; Jacobian JVol; } } }
  						{ Name ey_d_re; Value { Local { [Re[  CompY[{u}] ]]; In Omega; Jacobian JVol; } } }
              { Name ez_d_re; Value { Local { [Re[  CompZ[{u}] ]]; In Omega; Jacobian JVol; } } }
              { Name ex_d_im; Value { Local { [Im[  CompX[{u}] ]]; In Omega; Jacobian JVol; } } }
              { Name ey_d_im; Value { Local { [Im[  CompY[{u}] ]]; In Omega; Jacobian JVol; } } }
              { Name ez_d_im; Value { Local { [Im[  CompZ[{u}] ]]; In Omega; Jacobian JVol; } } }

              { Name ex_cm_re; Value { Local { [Re[  CompX[Ecm[]] ]]; In Omega; Jacobian JVol; } } }
              { Name ex_cm_im; Value { Local { [Im[  CompX[Ecm[]] ]]; In Omega; Jacobian JVol; } } }
              { Name ey_cm_re; Value { Local { [Re[  CompY[Ecm[]] ]]; In Omega; Jacobian JVol; } } }
              { Name ey_cm_im; Value { Local { [Im[  CompY[Ecm[]] ]]; In Omega; Jacobian JVol; } } }
              { Name ez_cm_re; Value { Local { [Re[  CompZ[Ecm[]] ]]; In Omega; Jacobian JVol; } } }
              { Name ez_cm_im; Value { Local { [Im[  CompZ[Ecm[]] ]]; In Omega; Jacobian JVol; } } }

              { Name e_cm_norm; Value { Local { [Norm[Ecm[]]]; In Omega; Jacobian JVol; } } }
              { Name etot_norm; Value { Local { [Norm[{u}]]; In Omega; Jacobian JVol; } } }

              { Name srcx_re; Value { Local { [Re[  CompX[source[]] ]]; In Omega; Jacobian JVol; } } }
              { Name srcx_im; Value { Local { [Im[  CompX[source[]] ]]; In Omega; Jacobian JVol; } } }
              { Name srcy_re; Value { Local { [Re[  CompY[source[]] ]]; In Omega; Jacobian JVol; } } }
              { Name srcy_im; Value { Local { [Im[  CompY[source[]] ]]; In Omega; Jacobian JVol; } } }
              { Name srcz_re; Value { Local { [Re[  CompZ[source[]] ]]; In Omega; Jacobian JVol; } } }
              { Name srcz_im; Value { Local { [Im[  CompZ[source[]] ]]; In Omega; Jacobian JVol; } } }

              { Name epsilon; Value { Local { [ CompXX[epsilon[]] ]; In Omega; Jacobian JVol; } } }
/*                { Name dummy2; Value { Local { [Re[  CompXX[epsilon1[]] ]]; In Omega; Jacobian JVol; } } }*/
/*//             H diffracted 3 components, Im and Re parts                 */
/*                { Name hx_d_re; Value { Local { [ -Im[CompX[ 1/omega0*1/mu0 * {Curl u} ] ] ]; In Omega; Jacobian JVol; } } }*/
/*                { Name hx_d_im; Value { Local { [  Re[CompX[ 1/omega0*1/mu0 * {Curl u} ] ] ]; In Omega; Jacobian JVol; } } }*/
/*                { Name hy_d_re; Value { Local { [ -Im[CompY[ 1/omega0*1/mu0 * {Curl u} ] ] ]; In Omega; Jacobian JVol; } } }*/
/*                { Name hy_d_im; Value { Local { [  Re[CompY[ 1/omega0*1/mu0 * {Curl u} ] ] ]; In Omega; Jacobian JVol; } } }*/
/*                { Name hz_d_re; Value { Local { [ -Im[CompZ[ 1/omega0*1/mu0 * {Curl u} ] ] ]; In Omega; Jacobian JVol; } } }*/
/*                { Name hz_d_im; Value { Local { [  Re[CompZ[ 1/omega0*1/mu0 * {Curl u} ] ] ]; In Omega; Jacobian JVol; } } }*/
/*//             Poynting diffracted field => used in reflexion coefficient (real quantity)*/
                // { Name pox; Value { Local { [ 0.5*CompX[ Im[ 1/omega0 * 1/mu0*{u}/\Conj[{Curl u}] ] ] ]; In Omega; Jacobian JVol; } } }*/
                // { Name poy; Value { Local { [ 0.5*CompY[ Im[ 1/omega0 * 1/mu0*{u}/\Conj[{Curl u}] ] ] ]; In Omega; Jacobian JVol; } } }
                // { Name poz; Value { Local { [ 0.5*CompZ[ Im[ 1/omega0 * 1/mu0*{u}/\Conj[{Curl u}] ] ] ]; In Omega; Jacobian JVol; } } }
/*//             To be integrated via GMSH Integrate Plugin (use the cutplane view) */
/*                { Name poz_r_int           ; Value {    Local { [ 1/period_x * 1/period_y * CompZ[ 0.5*Re[ Complex[0,1]* 1/omega0 * 1/mu0 * {u}/\Conj[{Curl u}] ] ] / Pinc ]; In Omega; Jacobian JVol; } } }*/
/*                { Name poz_t_int           ; Value {    Local { [ 1/period_x * 1/period_y * CompZ[ 0.5*Re[ (Einc[]+{u})/\Conj[ Complex[0,1]/(omega0*mu0)*{Curl u} + Hinc[]] ] ] / Pinc ]; In Omega; Jacobian JVol; } } }*/
/*                { Name joule_losses        ; Value {    Local { [  epsilon0*omega0 * 0.5*Fabs[Im[epsilon_rode[]]]*(SquNorm[{u}+Einc[]]) / (Pinc*period_x*period_y)  ] ; In Omega; Jacobian JVol; } } }*/
/*//             R (using E_d) then T (using E_d+E_tot). Makes use of the fact that we know a priori that the efficiencies are constant along the substrate */
/*                {  Name int_vol_div_P_Super   ; Value { Integral { [ 1/Pinc*1/period_x*1/period_y*1/(h_supb-h_supa)*CompZ[ 0.5*Re[ Complex[0,1]* 1/omega0 * 1/mu0 * {u}/\Conj[{Curl u}] ] ]   ]             ; In Super ; Integration Int_1 ; Jacobian JVol ; } } }*/
/*                {  Name int_vol_div_P_Subs    ; Value { Integral { [ 1/Pinc*1/period_x*1/period_y*1/(h_subb-h_suba)*CompZ[ 0.5*Re[ (Einc[]+{u})/\Conj[ Complex[0,1]/(omega0*mu0)*{Curl u} + Hinc[]] ] ]   ] ; In Subs  ; Integration Int_1 ; Jacobian JVol ; } } }*/
/*//             R (using E_d) then T (using E_d+E_tot) with 2 plane cuts in the superstrate and substrate*/
/*//                 {  Name int_sur_div_P_Sup1 ; Value { Integral { [ 1/Pinc*1/period_x*1/period_y*CompZ[ 0.5*Re[ Complex[0,1]* 1/omega0 * 1/mu0 * {u}/\Conj[{Curl u}] ] ]   ]             ; In SurfIntegSup1 ; Integration Int_1 ; Jacobian JSur ; } } }*/
/*//                 {  Name int_sur_div_P_Sup2 ; Value { Integral { [ 1/Pinc*1/period_x*1/period_y*CompZ[ 0.5*Re[ Complex[0,1]* 1/omega0 * 1/mu0 * {u}/\Conj[{Curl u}] ] ]   ]             ; In SurfIntegSup2 ; Integration Int_1 ; Jacobian JSur ; } } }*/
/*//                 {  Name int_sur_div_P_Sub1 ; Value { Integral { [ 1/Pinc*1/period_x*1/period_y*CompZ[ 0.5*Re[ (Einc[]+{u})/\Conj[ Complex[0,1]/(omega0*mu0)*{Curl u} + Hinc[]] ] ]   ] ; In SurfIntegSub1 ; Integration Int_1 ; Jacobian JSur ; } } }*/
/*//                 {  Name int_sur_div_P_Sub2 ; Value { Integral { [ 1/Pinc*1/period_x*1/period_y*CompZ[ 0.5*Re[ (Einc[]+{u})/\Conj[ Complex[0,1]/(omega0*mu0)*{Curl u} + Hinc[]] ] ]   ] ; In SurfIntegSub2 ; Integration Int_1 ; Jacobian JSur ; } } }*/
/*//             Computes the volume of a region */
/*                {  Name blabla              ; Value { Integral { [1] ; In Rode ; Integration Int_1 ; Jacobian JVol ; } } }*/

}}}




PostOperation {


    { Name postopQ; NameOfPostProcessing get_Ed ;
        Operation {
	    Print[ normalized_losses1[Omega_source] , OnGlobal, File "Q.txt", Format Table ];
        }
    }
		{ Name postop_epsilon; NameOfPostProcessing get_Ed ;
				Operation {
			/* Print [ epsilon , OnElementsOf Omega, File "epsilon.pos"]; */
			Print [ ex_d_re  , OnElementsOf Omega, File "ex_d_re.pos"];
		/* Print [ source , OnElementsOf Omega, File "source.pos"]; */
	/* Print [ Ecal , OnElementsOf Omega, File "Ecal.pos"]; */
		
			
				}
		}



{ Name Ed; NameOfPostProcessing get_Ed ;
    Operation {

 			      // Print [ ex_re_r  , OnElementsOf Omega, File "Ex_re_r.pos", Smoothing];
			      // Print [ ey_re_r  , OnElementsOf Omega, File "Ey_re_r.pos", Smoothing];
			      // Print [ ez_re_r  , OnElementsOf Omega, File "Ez_re_r.pos", Smoothing];
			      // Print [ ex_im_r  , OnElementsOf Omega, File "Ex_im_r.pos", Smoothing];
			      // Print [ ey_im_r  , OnElementsOf Omega, File "Ey_im_r.pos", Smoothing];
			      // Print [ ez_im_r  , OnElementsOf Omega, File "Ez_im_r.pos", Smoothing];
			      //
			      // Print [ ex_re_t  , OnElementsOf Omega, File "Ex_re_t.pos", Smoothing];
			      // Print [ ey_re_t  , OnElementsOf Omega, File "Ey_re_t.pos", Smoothing];
			      // Print [ ez_re_t  , OnElementsOf Omega, File "Ez_re_t.pos", Smoothing];
			      // Print [ ex_im_t  , OnElementsOf Omega, File "Ex_im_t.pos", Smoothing];
			      // Print [ ey_im_t  , OnElementsOf Omega, File "Ey_im_t.pos", Smoothing];
			      // Print [ ez_im_t  , OnElementsOf Omega, File "Ez_im_t.pos", Smoothing];
			      // Print [ e_cm_norm  , OnElementsOf Omega, File "e_cm_norm.pos", Smoothing];
			      // Print [ etot_norm  , OnElementsOf Omega, File "etot_norm.pos", Smoothing];

			 /*Print [ etot_norm  , OnElementsOf Omega, File "etot_norm.pos"];*/
       /*Print [ Etot , OnElementsOf Omega, File "Etot.pos"];*/
      // Print [ Edif , OnElementsOf Omega, File "Edif.pos"];
      // Print [ Ecm  , OnElementsOf Omega, File "Ecm.pos"];

		For i In {0:nb_slice-1}
			// Print [ ex_t , OnPlane { {-period_x/2,-period_y/2,zcut_sub~{i}} {-period_x/2,period_y/2,zcut_sub~{i}} {period_x/2,-period_y/2,zcut_sub~{i}} } {ninterv_integ,ninterv_integ} , File "Ex_t_XYcut.out",Format Table];
			// Print [ ey_t , OnPlane { {-period_x/2,-period_y/2,zcut_sub~{i}} {-period_x/2,period_y/2,zcut_sub~{i}} {period_x/2,-period_y/2,zcut_sub~{i}} } {ninterv_integ,ninterv_integ} , File "Ey_t_XYcut.out",Format Table];
			// Print [ ez_t , OnPlane { {-period_x/2,-period_y/2,zcut_sub~{i}} {-period_x/2,period_y/2,zcut_sub~{i}} {period_x/2,-period_y/2,zcut_sub~{i}} } {ninterv_integ,ninterv_integ} , File "Ez_t_XYcut.out",Format Table];
			// Print [ ex_r , OnPlane { {-period_x/2,-period_y/2,zcut_sup~{i}} {-period_x/2,period_y/2,zcut_sup~{i}} {period_x/2,-period_y/2,zcut_sup~{i}} } {ninterv_integ,ninterv_integ} , File "Ex_r_XYcut.out",Format Table];
			// Print [ ey_r , OnPlane { {-period_x/2,-period_y/2,zcut_sup~{i}} {-period_x/2,period_y/2,zcut_sup~{i}} {period_x/2,-period_y/2,zcut_sup~{i}} } {ninterv_integ,ninterv_integ} , File "Ey_r_XYcut.out",Format Table];
			// Print [ ez_r , OnPlane { {-period_x/2,-period_y/2,zcut_sup~{i}} {-period_x/2,period_y/2,zcut_sup~{i}} {period_x/2,-period_y/2,zcut_sup~{i}} } {ninterv_integ,ninterv_integ} , File "Ez_r_XYcut.out",Format Table];
  		Print [ Etot , OnPlane { {-period_x/2,-period_y/2,zcut_sub~{i}} {-period_x/2,period_y/2,zcut_sub~{i}} {period_x/2,-period_y/2,zcut_sub~{i}} } {ninterv_integ,ninterv_integ} , File > "Etot_XYcut.out",Smoothing ,Format Table];
  		Print [ Edif , OnPlane { {-period_x/2,-period_y/2,zcut_sup~{i}} {-period_x/2,period_y/2,zcut_sup~{i}} {period_x/2,-period_y/2,zcut_sup~{i}} } {ninterv_integ,ninterv_integ} , File > "Edif_XYcut.out",Smoothing, Format Table];
				  		// Print [ ex_re_t , OnPlane { {-period_x/2,-period_y/2,zcut_sub~{i}} {-period_x/2,period_y/2,zcut_sub~{i}} {period_x/2,-period_y/2,zcut_sub~{i}} } {ninterv_integ,ninterv_integ} , File Sprintf("Onplane_ex_re_t_XYcut_%02g.pos", i)];
				  		// Print [ ex_im_t , OnPlane { {-period_x/2,-period_y/2,zcut_sup~{i}} {-period_x/2,period_y/2,zcut_sup~{i}} {period_x/2,-period_y/2,zcut_sup~{i}} } {ninterv_integ,ninterv_integ} , File Sprintf("Onplane_ex_im_t_XYcut_%02g.pos", i)];
		EndFor
// 		Print[ normalized_losses1[L_3]           , OnGlobal, File "temp-Q1.txt", Format Table ];


			//Print [ Ecm  , OnElementsOf Omega, File "Ecm.pos", Smoothing];
			//Print [ Etot  , OnElementsOf Omega, File "Etot.pos", Smoothing];
			//Print [ Ecal  , OnElementsOf Omega, File "Ecal.pos", Smoothing];
			//Print [ Ecal_x  , OnElementsOf Omega, File "Ecal_x.pos", Smoothing];
			//Print [ Ecal_y  , OnElementsOf Omega, File "Ecal_y.pos", Smoothing];
			//Print [ Ecal_z  , OnElementsOf Omega, File "Ecal_z.pos", Smoothing];












      // Print [ cooX , OnPlane { {-period_x/2,-period_y/2,zcut_sub~{i}} {-period_x/2,period_y/2,zcut_sub~{i}} {period_x/2,-period_y/2,zcut_sub~{i}}  } {ninterv_integ,ninterv_integ} , File "cooX.out",Format Table];
      // Print [ cooY , OnPlane { {-period_x/2,-period_y/2,zcut_sub~{i}} {-period_x/2,period_y/2,zcut_sub~{i}} {period_x/2,-period_y/2,zcut_sub~{i}}  } {ninterv_integ,ninterv_integ} , File "cooY.out",Format Table];
      // Print [ cooXpYs2 , OnPlane { {-period_x/2,-period_y/2,zcut_sub~{i}} {-period_x/2,period_y/2,zcut_sub~{i}} {period_x/2,-period_y/2,zcut_sub~{i}}  } {ninterv_integ,ninterv_integ} , File "cooXpYs2.out",Format Table];
			// VERIFS
        // Print [ complex  , OnElementsOf Omega, File "complex.pos"];

        // Print [ complex , OnPlane { {0,0,0} {1,0,0} {0,1,0} } {50,50} , File "complex.out",Format Table];
        /*Print [ complex  , OnElementsOf Omega, File "complex.pos"];*/


        //Print[ int_vol_div_P_Subs[Subs]          , OnGlobal, File "temp.txt"  , Format Table ];
        //Print[ int_vol_div_P_Super[Super]        , OnGlobal, File > "temp.txt", Format Table ];
      // Print [ ex_cm_re  , OnElementsOf Omega, File "ex_cm_re.pos"];
      // Print [ ey_cm_re  , OnElementsOf Omega, File "ey_cm_re.pos"];
      // Print [ ez_cm_re  , OnElementsOf Omega, File "ez_cm_re.pos"];
      // Print [ ex_cm_im  , OnElementsOf Omega, File "ex_cm_im.pos"];
      // Print [ ey_cm_im  , OnElementsOf Omega, File "ey_cm_im.pos"];
      // Print [ ez_cm_im  , OnElementsOf Omega, File "ez_cm_im.pos"];
      //
      /* Print [ srcx_re  , OnElementsOf Omega, File "srcx_re.pos"];
      Print [ srcx_im  , OnElementsOf Omega, File "srcx_im.pos"];
      Print [ srcy_re  , OnElementsOf Omega, File "srcy_re.pos"];
      Print [ srcy_im  , OnElementsOf Omega, File "srcy_im.pos"];
      Print [ srcz_re  , OnElementsOf Omega, File "srcz_re.pos"];
      Print [ srcz_im  , OnElementsOf Omega, File "srcz_im.pos"]; */
/* 
			Print [ ex_d_re  , OnElementsOf Omega, File "ex_d_re.pos"];
      Print [ ey_d_re  , OnElementsOf Omega, File "ey_d_re.pos"];
      Print [ ez_d_re  , OnElementsOf Omega, File "ez_d_re.pos"];
      Print [ ex_d_im  , OnElementsOf Omega, File "ex_d_im.pos"];
      Print [ ey_d_im  , OnElementsOf Omega, File "ey_d_im.pos"];
      Print [ ez_d_im  , OnElementsOf Omega, File "ez_d_im.pos"]; */
/*        Print [ poz_r_int, OnElementsOf Omega, File "poz_r_int.pos"]; //, Smoothing ];*/
/*        Print [ poz_t_int, OnElementsOf Omega, File "poz_t_int.pos"]; //, Smoothing ];*/
//       Print [ dummy1  , OnElementsOf Omega, File "dummy1.pos"];
//        Print [ dummy2  , OnElementsOf Omega, File "dummy2.pos"];

/*        Print [ ex_d_im  , OnElementsOf Omega, File "map_imEX_diffacted.pos"];*/
/*        Print [ ey_d_im  , OnElementsOf Omega, File "map_imEY_diffacted.pos"];*/
/*        Print [ ez_d_im  , OnElementsOf Omega, File "map_imEZ_diffacted.pos"];*/
/*        Print [ ex_cm_re  , OnElementsOf Omega, File "map_reEX_cm.pos"];*/
/*        Print [ ey_cm_re  , OnElementsOf Omega, File "map_reEY_cm.pos"];*/
/*        Print [ ez_cm_re  , OnElementsOf Omega, File "map_reEZ_cm.pos"];*/
/*        Print [ ex_cm_im  , OnElementsOf Omega, File "map_imEX_cm.pos"];*/
/*        Print [ ey_cm_im  , OnElementsOf Omega, File "map_imEY_cm.pos"];*/
/*        Print [ ez_cm_im  , OnElementsOf Omega, File "map_imEZ_cm.pos"];*/
/*        Print [ hy_d_re  , OnElementsOf Omega, File "map_reHY_diffacted.pos"]; //, Smoothing ];*/
        }
}
}
