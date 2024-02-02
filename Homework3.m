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
EI_24 = Eb * t2^2/12;
EA_t = (3+0.15*W)*10^5;

%% Loading
p = 80 + 2*Y;
Q = 210 + 5*Z;

%% Spring 
ks = (5+0.5*Y)*10^3;

%% Element 1 - Bernoulli beam
L1 = b;

K1 = BernoulliElementStiffnessMatrix(EA_24, EI_24,L1);
K1red = applyBC(K1,0,0,0,1,1,0);

C11 = cosd(90);
C21 = cosd(180);
TM1 = getTransformationMatrix(C11,C21);
TM1red = applyBCtransformationMatrix(TM1,0,0,0,1,1,0);

K1red = TM1red' * K1red * TM1red;

%% Element 1 - Spring applied at a generic point along the beam
ksi = -1/3;
N = (1/4) * (1+ksi)^2*(2-ksi);

Kspring = N * ks*N;


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



