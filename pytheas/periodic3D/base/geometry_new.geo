
Include "parameters.dat";


lc_pmltop       = lambda0/(parmesh_pml*Sqrt[Fabs[eps_re_L_1]]);
lc_layer_1       = lambda0/(parmesh*Sqrt[Fabs[eps_re_L_1]]);
lc_layer_2       = lambda0/(parmesh_des*Sqrt[Fabs[eps_re_L_2]]);
lc_layer_3       = lambda0/(parmesh_des*Sqrt[Fabs[eps_re_L_3]]);
lc_layer_4       = lambda0/(parmesh*Sqrt[Fabs[eps_re_L_4]]);
lc_layer_5       = lambda0/(parmesh*Sqrt[Fabs[eps_re_L_5]]);
lc_layer_6       = lambda0/(parmesh*Sqrt[Fabs[eps_re_L_6]]);
lc_pmlbot       = lambda0/(parmesh_pml*Sqrt[Fabs[eps_re_L_6]]);


Point(1)  = {-period_x/2,-period_y/2, PML_bot_hh, lc_pmltop};
Point(2)  = {-period_x/2, period_y/2, PML_bot_hh, lc_pmltop};
Point(3)  = { period_x/2, period_y/2, PML_bot_hh, lc_pmltop};
Point(4)  = { period_x/2,-period_y/2, PML_bot_hh, lc_pmltop};
Point(5)  = {-period_x/2,-period_y/2, hh_L_6, lc_layer_6};
Point(6)  = {-period_x/2, period_y/2, hh_L_6, lc_layer_6};
Point(7)  = { period_x/2, period_y/2, hh_L_6, lc_layer_6};
Point(8)  = { period_x/2,-period_y/2, hh_L_6, lc_layer_6};
Point(9)  = {-period_x/2,-period_y/2, hh_L_5, lc_layer_5};
Point(10) = {-period_x/2, period_y/2, hh_L_5, lc_layer_5};
Point(11) = { period_x/2, period_y/2, hh_L_5, lc_layer_5};
Point(12) = { period_x/2,-period_y/2, hh_L_5, lc_layer_5};
Point(13) = {-period_x/2,-period_y/2, hh_L_4, lc_layer_4};
Point(14) = {-period_x/2, period_y/2, hh_L_4, lc_layer_4};
Point(15) = { period_x/2, period_y/2, hh_L_4, lc_layer_4};
Point(16) = { period_x/2,-period_y/2, hh_L_4, lc_layer_4};
Point(17) = {-period_x/2,-period_y/2, hh_L_3, lc_layer_3};
Point(18) = {-period_x/2, period_y/2, hh_L_3, lc_layer_3};
Point(19) = { period_x/2, period_y/2, hh_L_3, lc_layer_3};
Point(20) = { period_x/2,-period_y/2, hh_L_3, lc_layer_3};
Point(21) = {-period_x/2,-period_y/2, hh_L_2, lc_layer_2};
Point(22) = {-period_x/2, period_y/2, hh_L_2, lc_layer_2};
Point(23) = { period_x/2, period_y/2, hh_L_2, lc_layer_2};
Point(24) = { period_x/2,-period_y/2, hh_L_2, lc_layer_2};
Point(25) = {-period_x/2,-period_y/2, hh_L_1, lc_layer_1};
Point(26) = {-period_x/2, period_y/2, hh_L_1, lc_layer_1};
Point(27) = { period_x/2, period_y/2, hh_L_1, lc_layer_1};
Point(28) = { period_x/2,-period_y/2, hh_L_1, lc_layer_1};
Point(29) = {-period_x/2,-period_y/2, PML_top_hh, lc_layer_1};
Point(30) = {-period_x/2, period_y/2, PML_top_hh, lc_layer_1};
Point(31) = { period_x/2, period_y/2, PML_top_hh, lc_layer_1};
Point(32) = { period_x/2,-period_y/2, PML_top_hh, lc_layer_1};
Point(33) = {-period_x/2,-period_y/2, PML_top_hh+PML_top, lc_pmltop};
Point(34) = {-period_x/2, period_y/2, PML_top_hh+PML_top, lc_pmltop};
Point(35) = { period_x/2, period_y/2, PML_top_hh+PML_top, lc_pmltop};
Point(36) = { period_x/2,-period_y/2, PML_top_hh+PML_top, lc_pmltop};


Line(1)  = {2, 1};
Line(2)  = {1, 4};
Line(3)  = {4, 3};
Line(4)  = {3, 2};
Line(5)  = {2, 6};
Line(6)  = {6, 5};
Line(7)  = {5, 8};
Line(8)  = {8, 7};
Line(9)  = {7, 6};
Line(10) = {6, 10};
Line(11) = {10, 9};
Line(12) = {9, 12};
Line(13) = {12, 11};
Line(14) = {11, 10};
Line(15) = {10, 14};
Line(16) = {14, 13};
Line(17) = {13, 16};
Line(18) = {16, 15};
Line(19) = {15, 14};
Line(20) = {14, 18};
Line(21) = {18, 17};
Line(22) = {17, 20};
Line(23) = {20, 19};
Line(24) = {19, 18};
Line(25) = {18, 22};
Line(26) = {22, 21};
Line(27) = {21, 24};
Line(28) = {24, 23};
Line(29) = {23, 22};
Line(30) = {22, 26};
Line(31) = {26, 25};
Line(32) = {25, 28};
Line(33) = {28, 27};
Line(34) = {27, 26};
Line(35) = {27, 23};
Line(36) = {23, 19};
Line(37) = {19, 15};
Line(38) = {15, 11};
Line(39) = {11, 7};
Line(40) = {7, 3};
Line(41) = {1, 5};
Line(42) = {5, 9};
Line(43) = {9, 13};
Line(44) = {13, 17};
Line(45) = {17, 21};
Line(46) = {21, 25};
Line(47) = {28, 24};
Line(48) = {24, 20};
Line(49) = {20, 16};
Line(50) = {16, 12};
Line(51) = {12, 8};
Line(52) = {8, 4};
Line(53) = {26,30};
Line(54) = {30,29};
Line(55) = {29,32};
Line(56) = {32,31};
Line(57) = {31,30};
Line(58) = {30,34};
Line(59) = {34,33};
Line(60) = {33,36};
Line(61) = {36,35};
Line(62) = {35,34};
Line(63) = {25,29};
Line(64) = {29,33};
Line(65) = {28,32};
Line(66) = {32,36};
Line(67) = {27,31};
Line(68) = {31,35};


rb=0.3;
rl=0.3;
Point(70)  = { 0  , 0  , hh_L_4 , lc_layer_1};  //Centre du plot élliptique L_4
Point(71)  = { rb , 0  , hh_L_4 , lc_layer_1};  //Point du grand axe du plot élliptique  (+X) L_4
Point(72)  = { 0  , rl , hh_L_4 , lc_layer_1};  //Point du petit axe du plot élliptique  (+Y) L_4
Point(73)  = {-rb , 0  , hh_L_4 , lc_layer_1};  //Point du grand axe du plot élliptique  (-X) L_4
Point(74)  = { 0  ,-rl , hh_L_4 , lc_layer_1};  //Point du petit axe du plot élliptique  (-Y) L_4

Ellipse(300) = {71,70,71,72};
Ellipse(301) = {72,70,71,73};
Ellipse(302) = {73,70,71,74};
Ellipse(303) = {74,70,71,71};

Line Loop(176) = {300, 301, 302, 303}; // base plot L_4
Plane Surface(178) = {176}; // base plot L_4

Extrude {0, 0, thick_L_4}{ Surface{178} ; }
Extrude {0, 0, thick_L_3}{ Surface{325} ; }


//Face 1
Line Loop(105) = {1, 2, 3, 4};
Line Loop(106) = {5, 6, -41, -1};
Plane Surface(107) = {106};
Line Loop(108) = {10, 11, -42, -6};
Plane Surface(109) = {108};
Line Loop(110) = {15, 16, -43, -11};
Plane Surface(111) = {110};
Line Loop(112) = {20, 21, -44, -16};
Plane Surface(113) = {112};
Line Loop(114) = {25, 26, -45, -21};
Plane Surface(115) = {114};
Line Loop(116) = {30, 31, -46, -26};
Plane Surface(117) = {116};
Line Loop(118) = {53, 54, -63, -31};
Plane Surface(119) = {118};
Line Loop(120) = {58, 59, -64, -54};
Plane Surface(121) = {120};

//Face 2
Line Loop(122) = {40, -3, -52, 8};
Plane Surface(123) = {122};
Line Loop(124) = {39, -8, -51, 13};
Plane Surface(125) = {124};
Line Loop(126) = {38, -13, -50, 18};
Plane Surface(127) = {126};
Line Loop(128) = {37, -18, -49, 23};
Plane Surface(129) = {128};
Line Loop(130) = {36, -23, -48, 28};
Plane Surface(131) = {130};
Line Loop(132) = {35, -28, -47, 33};
Plane Surface(133) = {132};
Line Loop(134) = {-67, -33, 65, 56};
Plane Surface(135) = {134};
Line Loop(136) = {-68, -56, 66, 61};
Plane Surface(137) = {136};

//Face 3
Line Loop(138) = {41, 7, 52, -2};
Plane Surface(139) = {138};
Line Loop(140) = {42, 12, 51, -7};
Plane Surface(141) = {140};
Line Loop(142) = {43, 17, 50, -12};
Plane Surface(143) = {142};
Line Loop(144) = {44, 22, 49, -17};
Plane Surface(145) = {144};
Line Loop(146) = {45, 27, 48, -22};
Plane Surface(147) = {146};
Line Loop(148) = {46, 32, 47, -27};
Plane Surface(149) = {148};
Line Loop(150) = {63, 55, -65, -32};
Plane Surface(151) = {150};
Line Loop(152) = {64, 60, -66, -55};
Plane Surface(153) = {152};


//Face 4
Line Loop(154) = {5, -9, 40, 4};
Plane Surface(155) = {154};
Line Loop(156) = {10, -14, 39, 9};
Plane Surface(157) = {156};
Line Loop(158) = {15, -19, 38, 14};
Plane Surface(159) = {158};
Line Loop(160) = {20, -24, 37, 19};
Plane Surface(161) = {160};
Line Loop(162) = {25, -29, 36, 24};
Plane Surface(163) = {162};
Line Loop(164) = {30, -34, 35, 29};
Plane Surface(165) = {164};
Line Loop(166) = {53, -57, -67, 34};
Plane Surface(167) = {166};
Line Loop(168) = {58, -62, -68, 57};
Plane Surface(169) = {168};

Periodic Surface (139) {41, 7, 52, -2}    = (155)  {5, -9, 40, 4};
Periodic Surface (141) {42, 12, 51, -7}   = (157)  {10, -14, 39, 9};
Periodic Surface (143) {43, 17, 50, -12}  = (159)  {15, -19, 38, 14};
Periodic Surface (145) {44, 22, 49, -17}  = (161)  {20, -24, 37, 19};
Periodic Surface (147) {45, 27, 48, -22}  = (163)  {25, -29, 36, 24};
Periodic Surface (149) {46, 32, 47, -27}  = (165)  {30, -34, 35, 29};
Periodic Surface (151) {63, 55,-65, -32}  = (167)  {53, -57, -67, 34};
Periodic Surface (153) {64, 60,-66, -55}  = (169)  {58, -62, -68, 57};

Periodic Surface (107) {5, 6, -41, -1}   = (123) {-40, -8, 52, 3};
Periodic Surface (109) {10, 11, -42, -6} = (125) {-39, -13, 51, 8};
Periodic Surface (111) {15, 16, -43, -11} = (127) {-38, -18, 50, 13};
Periodic Surface (113) {20, 21, -44, -16} = (129) {-37, -23, 49, 18};
Periodic Surface (115) {25, 26, -45, -21} = (131) {-36, -28, 48, 23};
Periodic Surface (117) {30, 31, -46, -26} = (133) {-35, -33, 47, 28};
Periodic Surface (119) {53, 54, -63, -31} = (135) {67, -56, -65, 33};
Periodic Surface (121) {58, 59, -64, -54} = (137) {68, -61, -66, 56};

Physical Line(10001) = {41,42,43,44,45,46,63,64,40,39,38,37,36,35,67,68,52,51,50,49,48,47,65,66,5,10,15,20,25,30,53,58,1,6,11,16,21,26,31,54,59,2,7,12,17,22,27,32,55,60,4,9,14,19,24,29,34,57,62,3,8,13,18,23,28,33,56,61};

// Bloch X - puis +
Physical Surface(750) = {107,109,111,113,115,117,119,121};
Physical Surface(751) = {123,125,127,129,131,133,135,137};

// Bloch Y - puis +
Physical Surface(760) = {139,141,143,145,147,149,151,153};
Physical Surface(761) = {155,157,159,161,163,165,167,169};

// Diri
Physical Surface(770) = {170,192};


Plane Surface(170) = {105};

Line Loop(171) = {6, 7, 8, 9};
Plane Surface(172) = {171};
Line Loop(173) = {11, 12, 13, 14};
Plane Surface(174) = {173};
Line Loop(175) = {16,17,18,19};
Plane Surface(177) = {175, 176}; // base L_4

Line Loop(187) = {31, 32, 33, 34};
Plane Surface(188) = {187};
Line Loop(189) = {54, 55, 56, 57};
Plane Surface(190) = {189};
Line Loop(191) = {59, 60, 61, 62};
Plane Surface(192) = {191};

Line Loop(10002) = {23, 24, 21, 22};
Line Loop(10003) = {305, 306, 307, 308};
Plane Surface(10004) = {10002, 10003};
Line Loop(10005) = {28, 29, 26, 27};
Line Loop(10006) = {328, 329, 330, 327};
Plane Surface(10007) = {10005, 10006};

Surface Loop(10008) = {139, 107, 155, 123, 170, 172};
Volume(10009) = {10008}; //PML_bot
Surface Loop(10010) = {172, 125, 157, 109, 141, 174};
Volume(10011) = {10010}; //L_6 substrat
Surface Loop(10012) = {174, 143, 111, 159, 127, 177, 178};
Volume(10013) = {10012};//L_5 accroche
Surface Loop(10014) = {177, 113, 161, 129, 145, 10004, 316, 320, 324, 312};
Volume(10015) = {10014};//L_4
Surface Loop(10016) = {10004, 115, 163, 131, 147, 10007, 342, 346, 334, 338};
Volume(10017) = {10016};//L_3
Surface Loop(10018) = {10007, 347, 149, 117, 165, 133, 188};
Volume(10019) = {10018};//L_2 couverture
Surface Loop(10020) = {188, 167, 119, 151, 135, 190};
Volume(10021) = {10020};//L_1 superstrat
Surface Loop(10022) = {190, 153, 121, 169, 192, 137};
Volume(10023) = {10022};//PML_top

Physical Volume(1000) = {10009};  // PML bot
Physical Volume(2000) = {10011};  // layer L6
Physical Volume(3000) = {10013};  // layer L5
Physical Volume(4000) = {10015};  // layer L4
Physical Volume(5000) = {10017};  // layer L3
Physical Volume(6000) = {10019};  // layer L2
Physical Volume(7000) = {10021};  // layer L1
Physical Volume(8000) = {10023};  // PML top
Physical Volume(9000) = {2};  // InSin1
Physical Volume(10000) = {1};  // InSin2

Mesh.Algorithm   = 1; // // 1=MeshAdapt, 5=Delaunay, 6=Frontal
Mesh.Algorithm3D = 2; // // 1=Delaunay, 4=Frontal
Mesh.Optimize = 1;
