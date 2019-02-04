// This code was created by pygmsh v4.3.6.
SetFactory("OpenCASCADE");
Include "parameters.dat";
lc_des = lambda0/(Sqrt[eps_des_re]*parmesh_des);
lc_host = lambda0/(Sqrt[eps_host_re]*parmesh);
lc_pml = lambda0/(Sqrt[eps_host_re]*parmesh_pml);
lc_sph = lambda0/(Sqrt[eps_sph_re]*parmesh_incl);
s522 = news;
Rectangle(s522) = {-hx_des/2, -hy_des/2, -hz_des/2, hx_des, hy_des};
pts_s522[] = PointsOf{Surface{s522};};
Characteristic Length{pts_s522[]} = lc_des;
ex1[] = Extrude{0,0,hz_des}{Surface{s522};};
s524 = news;
Rectangle(s524) = {-hx_box/2, -hy_box/2, -hz_box/2, hx_box, hy_box};
pts_s524[] = PointsOf{Surface{s524};};
Characteristic Length{pts_s524[]} = lc_host;
ex2[] = Extrude{0,0,hz_box}{Surface{s524};};
s526 = news;
Rectangle(s526) = {-h_pml-hx_box/2, -h_pml-hy_box/2, -h_pml-hz_box/2, h_pml, h_pml};
pts_s526[] = PointsOf{Surface{s526};};
Characteristic Length{pts_s526[]} = lc_pml;
ex3[] = Extrude{0,0,h_pml}{Surface{s526};};
s528 = news;
Rectangle(s528) = {-h_pml-hx_box/2, -h_pml-hy_box/2, -hz_box/2, h_pml, h_pml};
pts_s528[] = PointsOf{Surface{s528};};
Characteristic Length{pts_s528[]} = lc_pml;
ex4[] = Extrude{0,0,hz_box}{Surface{s528};};
s530 = news;
Rectangle(s530) = {-h_pml-hx_box/2, -h_pml-hy_box/2, hz_box/2, h_pml, h_pml};
pts_s530[] = PointsOf{Surface{s530};};
Characteristic Length{pts_s530[]} = lc_pml;
ex5[] = Extrude{0,0,h_pml}{Surface{s530};};
s532 = news;
Rectangle(s532) = {hx_box/2, -h_pml-hy_box/2, -h_pml-hz_box/2, h_pml, h_pml};
pts_s532[] = PointsOf{Surface{s532};};
Characteristic Length{pts_s532[]} = lc_pml;
ex6[] = Extrude{0,0,h_pml}{Surface{s532};};
s534 = news;
Rectangle(s534) = {hx_box/2, -h_pml-hy_box/2, -hz_box/2, h_pml, h_pml};
pts_s534[] = PointsOf{Surface{s534};};
Characteristic Length{pts_s534[]} = lc_pml;
ex7[] = Extrude{0,0,hz_box}{Surface{s534};};
s536 = news;
Rectangle(s536) = {hx_box/2, -h_pml-hy_box/2, hz_box/2, h_pml, h_pml};
pts_s536[] = PointsOf{Surface{s536};};
Characteristic Length{pts_s536[]} = lc_pml;
ex8[] = Extrude{0,0,h_pml}{Surface{s536};};
s538 = news;
Rectangle(s538) = {-h_pml-hx_box/2, hy_box/2, -h_pml-hz_box/2, h_pml, h_pml};
pts_s538[] = PointsOf{Surface{s538};};
Characteristic Length{pts_s538[]} = lc_pml;
ex9[] = Extrude{0,0,h_pml}{Surface{s538};};
s540 = news;
Rectangle(s540) = {-h_pml-hx_box/2, hy_box/2, -hz_box/2, h_pml, h_pml};
pts_s540[] = PointsOf{Surface{s540};};
Characteristic Length{pts_s540[]} = lc_pml;
ex10[] = Extrude{0,0,hz_box}{Surface{s540};};
s542 = news;
Rectangle(s542) = {-h_pml-hx_box/2, hy_box/2, hz_box/2, h_pml, h_pml};
pts_s542[] = PointsOf{Surface{s542};};
Characteristic Length{pts_s542[]} = lc_pml;
ex11[] = Extrude{0,0,h_pml}{Surface{s542};};
s544 = news;
Rectangle(s544) = {hx_box/2, hy_box/2, -h_pml-hz_box/2, h_pml, h_pml};
pts_s544[] = PointsOf{Surface{s544};};
Characteristic Length{pts_s544[]} = lc_pml;
ex12[] = Extrude{0,0,h_pml}{Surface{s544};};
s546 = news;
Rectangle(s546) = {hx_box/2, hy_box/2, -hz_box/2, h_pml, h_pml};
pts_s546[] = PointsOf{Surface{s546};};
Characteristic Length{pts_s546[]} = lc_pml;
ex13[] = Extrude{0,0,hz_box}{Surface{s546};};
s548 = news;
Rectangle(s548) = {hx_box/2, hy_box/2, hz_box/2, h_pml, h_pml};
pts_s548[] = PointsOf{Surface{s548};};
Characteristic Length{pts_s548[]} = lc_pml;
ex14[] = Extrude{0,0,h_pml}{Surface{s548};};
s550 = news;
Rectangle(s550) = {-hx_box/2, -h_pml-hy_box/2, -h_pml-hz_box/2, hx_box, h_pml};
pts_s550[] = PointsOf{Surface{s550};};
Characteristic Length{pts_s550[]} = lc_pml;
ex15[] = Extrude{0,0,h_pml}{Surface{s550};};
s552 = news;
Rectangle(s552) = {-hx_box/2, -h_pml-hy_box/2, -hz_box/2, hx_box, h_pml};
pts_s552[] = PointsOf{Surface{s552};};
Characteristic Length{pts_s552[]} = lc_pml;
ex16[] = Extrude{0,0,hz_box}{Surface{s552};};
s554 = news;
Rectangle(s554) = {-hx_box/2, -h_pml-hy_box/2, hz_box/2, hx_box, h_pml};
pts_s554[] = PointsOf{Surface{s554};};
Characteristic Length{pts_s554[]} = lc_pml;
ex17[] = Extrude{0,0,h_pml}{Surface{s554};};
s556 = news;
Rectangle(s556) = {-hx_box/2, hy_box/2, -h_pml-hz_box/2, hx_box, h_pml};
pts_s556[] = PointsOf{Surface{s556};};
Characteristic Length{pts_s556[]} = lc_pml;
ex18[] = Extrude{0,0,h_pml}{Surface{s556};};
s558 = news;
Rectangle(s558) = {-hx_box/2, hy_box/2, -hz_box/2, hx_box, h_pml};
pts_s558[] = PointsOf{Surface{s558};};
Characteristic Length{pts_s558[]} = lc_pml;
ex19[] = Extrude{0,0,hz_box}{Surface{s558};};
s560 = news;
Rectangle(s560) = {-hx_box/2, hy_box/2, hz_box/2, hx_box, h_pml};
pts_s560[] = PointsOf{Surface{s560};};
Characteristic Length{pts_s560[]} = lc_pml;
ex20[] = Extrude{0,0,h_pml}{Surface{s560};};
s562 = news;
Rectangle(s562) = {-h_pml-hx_box/2, -hy_box/2, -h_pml-hz_box/2, h_pml, hy_box};
pts_s562[] = PointsOf{Surface{s562};};
Characteristic Length{pts_s562[]} = lc_pml;
ex21[] = Extrude{0,0,h_pml}{Surface{s562};};
s564 = news;
Rectangle(s564) = {-h_pml-hx_box/2, -hy_box/2, -hz_box/2, h_pml, hy_box};
pts_s564[] = PointsOf{Surface{s564};};
Characteristic Length{pts_s564[]} = lc_pml;
ex22[] = Extrude{0,0,hz_box}{Surface{s564};};
s566 = news;
Rectangle(s566) = {-h_pml-hx_box/2, -hy_box/2, hz_box/2, h_pml, hy_box};
pts_s566[] = PointsOf{Surface{s566};};
Characteristic Length{pts_s566[]} = lc_pml;
ex23[] = Extrude{0,0,h_pml}{Surface{s566};};
s568 = news;
Rectangle(s568) = {hx_box/2, -hy_box/2, -h_pml-hz_box/2, h_pml, hy_box};
pts_s568[] = PointsOf{Surface{s568};};
Characteristic Length{pts_s568[]} = lc_pml;
ex24[] = Extrude{0,0,h_pml}{Surface{s568};};
s570 = news;
Rectangle(s570) = {hx_box/2, -hy_box/2, -hz_box/2, h_pml, hy_box};
pts_s570[] = PointsOf{Surface{s570};};
Characteristic Length{pts_s570[]} = lc_pml;
ex25[] = Extrude{0,0,hz_box}{Surface{s570};};
s572 = news;
Rectangle(s572) = {hx_box/2, -hy_box/2, hz_box/2, h_pml, hy_box};
pts_s572[] = PointsOf{Surface{s572};};
Characteristic Length{pts_s572[]} = lc_pml;
ex26[] = Extrude{0,0,h_pml}{Surface{s572};};
s574 = news;
Rectangle(s574) = {-hx_box/2, -hy_box/2, -h_pml-hz_box/2, hx_box, hy_box};
pts_s574[] = PointsOf{Surface{s574};};
Characteristic Length{pts_s574[]} = lc_pml;
ex27[] = Extrude{0,0,h_pml}{Surface{s574};};
s576 = news;
Rectangle(s576) = {-hx_box/2, -hy_box/2, hz_box/2, hx_box, hy_box};
pts_s576[] = PointsOf{Surface{s576};};
Characteristic Length{pts_s576[]} = lc_pml;
ex28[] = Extrude{0,0,h_pml}{Surface{s576};};
vol18 = newv;
Sphere(vol18) = {x_sph, y_sph, z_sph, R_sph};
pts_vol18[] = PointsOf{Volume{vol18};};
Characteristic Length{pts_vol18[]} = lc_sph;
bo1[] = BooleanDifference{ Volume{ex2[1]}; Delete; } { Volume{vol18};Volume{ex1[1]}; Delete;};
vol19 = newv;
Sphere(vol19) = {x_sph, y_sph, z_sph, R_sph};
pts_vol19[] = PointsOf{Volume{vol19};};
Characteristic Length{pts_vol19[]} = lc_sph;
s578 = news;
Rectangle(s578) = {-hx_des/2, -hy_des/2, -hz_des/2, hx_des, hy_des};
pts_s578[] = PointsOf{Surface{s578};};
Characteristic Length{pts_s578[]} = lc_des;
ex29[] = Extrude{0,0,hz_des}{Surface{s578};};
p9 = newp;
Point(p9) = {0, 0, 0};
Physical Volume(1) = {bo1[]};
Physical Volume(2) = {ex29[1]};
Physical Volume(3) = {vol19};
Physical Volume(4) = {ex3[1], ex5[1], ex6[1], ex8[1], ex9[1], ex11[1], ex12[1], ex14[1]};
Physical Volume(5) = {ex4[1], ex7[1], ex10[1], ex13[1]};
Physical Volume(6) = {ex15[1], ex17[1], ex18[1], ex20[1]};
Physical Volume(7) = {ex21[1], ex23[1], ex24[1], ex26[1]};
Physical Volume(8) = {ex22[1], ex25[1]};
Physical Volume(9) = {ex16[1], ex19[1]};
Physical Volume(10) = {ex27[1], ex28[1]};
Physical Point(11) = {p9};
Coherence;
Coherence;
Coherence;
Mesh.Algorithm=6;
Mesh.Algorithm3D=4;