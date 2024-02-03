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
iNum = '108022249956';
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
L1 = b;

% Area A as the function of Ksi
syms ksi 
side = t1 + ((t2-t1)/L1)*L1*(1+ksi)/2;
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

% Element matrix for element 1 - bending part
K1_loc_bending = Bw' * Eb * I * Bw * J;
K1_loc_bending = int(K1_loc_bending,-1,1);
K1_loc_bending_reduced = eval(K1_loc_bending(1,1));

% Local element matrix 
K1_loc = [K1_loc_axial_reduced 0; 0 K1_loc_bending_reduced ];

% Transformation matrix
C11 = cosd(90);
C21 = cosd(180);
TM1 = getTransformationMatrix(C11,C21);
TM1red = applyBCtransformationMatrix(TM1,0,0,0,1,1,0);

% Global element matrix - element 1
K1red = TM1red' * K1_loc * TM1red;

%% Element 1 - Spring applied at a generic point along the beam
ksi = -1/3;
N = (1/4) * (1+ksi)^2*(2-ksi);

Kspring_1 = N * ks*N;


%% Element 2 - Bernoulli beam
L2 = sqrt(a^2 + d^2);
K2 = BernoulliElementStiffnessMatrix(EA_24,EI_24,L2);
K2red = applyBC(K2,1,1,0,1,1,1);

alfa12 = atand(a/d);
alfa22 = 270-alfa12;
C12 = cosd(alfa12);
C22 = cosd(alfa22);
TM2 = getTransformationMatrix(C12,C22);
TM2red = applyBCtransformationMatrix(TM2,1,1,0,1,1,1);

K2red = TM2red' * K2red *TM2red;

%% Element 3 - truss
L3 = a + b + c;
K3_local = EA_t/L3 * [1 -1; -1 1];

% Directional cosines 
x13 = 0; z13 = 0;
x14 = 0; z14 = L3;
C13 = (x14 - x13)/L3;
C23 = (z14 - z13)/L3;

% Transformation matrix
TM3 = [C13 C23 0 0; 0 0 C13 C23];
K3_global = TM3' * K3_local * TM3;
K3_global_reduced = K3_global(1:2,1:2);

%% Element 4 - beam
L4 = sqrt(a^2 + d^2);
K4 = BernoulliElementStiffnessMatrix(EA_24,EI_24,L2);
K4red = applyBC(K2,1,1,1,1,1,0);

alfa14 = 360-atand(a/d);
alfa24 = 270-alfa14;
C14 = cosd(alfa14);
C24 = cosd(alfa24);
TM4 = getTransformationMatrix(C14,C24);
TM4red = applyBCtransformationMatrix(TM4,1,1,1,1,1,0);

K4red = TM4red' * K4red *TM4red;

%% Element 5 - beam

L5 = b;

% Area A as the function of Ksi
syms ksi 
side = t2 + ((t1-t2)/L5)*L5*(1+ksi)/2;
A = expand(side^2);
A = vpa(A,4);

% B matrix for u 
Bu = [-1/L5 1/L5];

% Jacobian 
J = L5/2;

% Element matrix for element 1 - axial part
K5_loc_axial = Bu' * Eb * A * Bu * J;
K5_loc_axial = int(K5_loc_axial, -1,1);
K5_loc_axial_reduced = eval(K5_loc_axial(1,1));

% Moment of inertia as a function of ksi
I = expand((side^4)/12);
I = vpa(I,4);

% B matrix for w
Bw = (2/L5)^2 * [-3*ksi/2 -(1+3*ksi)*L1/4];

% Element matrix for element 5 - bending part
K5_loc_bending = Bw' * Eb * I * Bw * J;
K5_loc_bending = int(K5_loc_bending,-1,1);
K5_loc_bending_reduced = eval(K5_loc_bending(1,1));

% Local element matrix 
K5_loc = [K5_loc_axial_reduced 0; 0 K5_loc_bending_reduced ];

% Transformation matrix
C15 = cosd(90);
C25 = cosd(180);
TM5 = getTransformationMatrix(C15,C25);
TM5red = applyBCtransformationMatrix(TM5,0,0,0,1,1,0);

% Global element matrix - element 5
K5red = TM5red' * K5_loc * TM5red;

%% Element 5 - Spring applied at a generic point along the beam
ksi = 1/3;
N = (1/4) * (1-ksi)^2*(2+ksi);

Kspring_5 = N * ks*N;
