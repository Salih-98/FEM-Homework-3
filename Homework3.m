%% Linear FEM - Homework 1 - Task 2
% --- Author: Salih Sahbegovic
% --- Date: 01.02.2024.
% --- Subject: Finite Element Methods in Linear Structural Mechanics
% --- Semester: Winter Semester 2023/2024

%% Defining data given in the task formulation 
clc
clear
addpath('calculationFunctions')
% In this section we first define W, X, Y, Z based on inputed
% immatriculation number, and then based on that value we define rest of
% given parameters a,b, c,d, A, E, pmax, prescribed values of F and u23, as
% well as spring stiffness 

% -- TODO: This is really rigid coding, just in order to perform
% calculations. Add more flexibility to the code later. 

% Extract W, X, Y, Z from the immatrikulation number
% TODO: Make it so the user can input it during execution, later.

iNum = '108022241234';
% Ali
%iNum = '108022249396';
% TO compare with Abdulkareem
%iNum = '108022241234';
[W, X, Y, Z] = getWXYZ(iNum);

%% Define geometry
a = 1 + 0.05*X;
b = 4 + 0.1*Z;
c = 0.1 + 0.2*Y;
d = 3 + 0.1*W;
t1 = 0.75 + 0.05*X - 0.04*Z;
t2 = t1/1.5;

%% Material
Eb = (4 + 0.25*X)*10^6;
EA_24 = Eb * t2^2;
EI_24 = (Eb * t2^4)/12;
EA_t = (3+0.15*W)*10^5;

%% Loading
p = 80 + 2*Y;
Q = 210 + 5*Z;

%% Spring 
ks = (5+0.5*Y)*10^3;

%% Element 1 - Bernoulli beam
% Length of element 1
L1 = b;

% Change of rectangle side as a function of ksi
syms ksi
side = t1 + ((t2-t1)/L1)*L1*(1+ksi)/2;

% Area A as the function of Ksi
A = expand(side^2);
A = vpa(A,4);

% B matrix for u 
Bu = [-1/L1 1/L1];

% Jacobian 
J = L1/2;

% Element matrix for element 1 - axial part
K1_loc_axial = Bu' * Eb * A * Bu * J;
K1_loc_axial = int(K1_loc_axial, -1,1);
K1_loc_axial_reduced = eval(K1_loc_axial(2,2));

% Moment of inertia as a function of ksi
I = expand((side^4)/12);
I = vpa(I,4);

% B matrix for w
Bw = (2/L1)^2 * [-3*ksi/2 -(1+3*ksi)*L1/4];

% Numerical integration with 4 gaus points

% Gaus point 1
w1 = 0.34785; ksi1 = -0.86114;
Bw_1 = double(subs(Bw,ksi,ksi1));
I_1 = double(subs(I,ksi,ksi1));
K1_loc_bending_1 = w1 * Bw' * Eb * I * Bw * J;
K1_loc_bending_1 = double(subs(K1_loc_bending_1,ksi,ksi1));

% Gaus point 2
w2 = 0.65241; ksi2 = -0.33998;
Bw_2 = double(subs(Bw,ksi,ksi2));
I_2 = double(subs(I,ksi,ksi2));
K1_loc_bending_2 = w2 * Bw' * Eb * I * Bw * J;
K1_loc_bending_2 = double(subs(K1_loc_bending_2,ksi,ksi2));

% Gaus point 3
w3 = 0.65241; ksi3 = 0.33998;
Bw_3 = double(subs(Bw,ksi,ksi3));
I_3 = double(subs(I,ksi,ksi3));
K1_loc_bending_3 = w3 * Bw' * Eb * I * Bw * J;
K1_loc_bending_3 = double(subs(K1_loc_bending_3,ksi,ksi3));

% Gaus point 4
w4 = 0.34785; ksi4 = 0.86114;
Bw_4 = double(subs(Bw,ksi,ksi4));
I_4 = double(subs(I,ksi,ksi4));
K1_loc_bending_4 = w4 * Bw' * Eb * I * Bw * J;
K1_loc_bending_4 = double(subs(K1_loc_bending_4,ksi,ksi4));

% Sum the matrices for gaus points
K1_loc_bending = K1_loc_bending_1 + K1_loc_bending_2 + K1_loc_bending_3 + K1_loc_bending_4;

% Local element matrix 
K1_loc = [K1_loc_axial_reduced 0 0; 0 K1_loc_bending(1,:); 0 K1_loc_bending(2,:)];

% Transformation matrix
C11 = cosd(90);
C21 = cosd(270-90);
TM1 = getTransformationMatrix(C11,C21);
TM1red = applyBCtransformationMatrix(TM1,0,0,0,1,1,1);

% Global element matrix - element 1
K1red = TM1red' * K1_loc * TM1red;

%% Element 1 - Spring applied at a generic point along the beam
ksiS = -1/3;
N = [(1/4) * (1+ksiS)^2*(2-ksiS) ((1+ksiS)^2)*(1-ksiS)*L1/8];

Kspring_1 = N' * ks * N;


%% Element 2 - Bernoulli beam
L2 = sqrt(a^2 + d^2);
K2 = BernoulliElementStiffnessMatrix(EA_24,EI_24,L2);
alfa12 = atand(a/d);
alfa22 = 270-alfa12;
C12 = cosd(alfa12);
C22 = cosd(alfa22);
TM2 = getTransformationMatrix(C12,C22);

K2 = TM2' * K2 *TM2;
K2 = applyBC(K2,1,1,1,0,1,1);
%% Element 3 - truss
L3 = a + b + c;
K3_local = EA_t/(L3*2) * [1 -1; -1 1];

% Directional cosines 
x13 = 0; z13 = 0;
x14 = 0; z14 = L3;
C13 = (x14 - x13)/L3;
C23 = (z14 - z13)/L3;

% Transformation matrix
TM3 = [C13 C23 0 0; 0 0 C13 C23];
K3_global = TM3' * K3_local * TM3;
K3_global_reduced = K3_global(2,2);


%% Element 5 - Spring applied at a generic point along the beam
ksiS = 1/3;
N = [(1/4) * (1-ksiS)^2*(2+ksiS) ((1-ksiS)^2)*(1+ksiS)*L1/8];

Kspring_5 = N' * ks * N;

%% Assembly stiffness matrix
K = K2;
K(1,1) = K(1,1)+K1red(1,1) + Kspring_1(1,1);
K(1,3) = K(1,3)+K1red(1,3) + Kspring_1(1,2);
K(2,2) = K(2,2)+K1red(2,2);
K(3,1) = K(3,1)+K1red(3,1) + Kspring_1(2,1);
K(3,3) = K(3,3)+K1red(3,3) + Kspring_1(2,2);
K(4,4) = K(4,4)+K3_global_reduced;
K_global_reduced = K;

%% Distributed load vector 
syms ksi
Nu1 = (1-ksi)/2;
Nu2 = (1+ksi)/2;
Nw1 = (1-ksi)^2 * (2+ksi)/4;
Nfi1 = -(1-ksi)^2 * (1+ksi)*L2/8;
Nw2 = (1+ksi)^2*(2-ksi)/4;
Nfi2 = (1+ksi)^2*(1-ksi)*L2/8;
J = L2/2;
N = [Nu1 0 0 Nu2 0 0; 0 Nw1 Nfi1 0 Nw2 Nfi2];
p = [0;p];
rb = N' * p  * J;
rb_local = double(int(rb,ksi,-1,1));
rb_local = [rb_local(1:3); rb_local(5:6)];

Bw = [(1-ksi)^2*(2+ksi)/4 -(1-ksi)^2*(1+ksi)*L2/8 (1+ksi)^2*(2-ksi)/4 (1+ksi)^2*(1-ksi)*L2/8];
rn_local = Bw' * Q;
rn_local = double(subs(rn_local,ksi,0));
rn_local = [0; rn_local(1:2,1); rn_local(3:4,1)];
r_local = rb_local + rn_local;

r_global = applyBC(TM2,1,1,1,0,1,1)' * r_local;

%% Displacements
u = inv(K_global_reduced)*r_global;

%% Curvature Kapa
Kapa = -(2/L2)^2 * [3*ksi/2 L2*(1-3*ksi)/4 -3*ksi/2 -(1+3*ksi)*L2/4];
ksiM = -1/(2+Z);
Kapa = double(subs(Kapa,ksi,ksiM));
Kapa = Kapa*u(2:5,1);

%% Moment at ksiM
M = EI_24 * Kapa;