function reducedStiffMatrix = applyBC(elementStiffnessMatrix,r1,r2,r3,r4,r5,r6)
% This is a function which reduces the element stiffness matrix by applying
% the homogenous Dirichlet Boundary conditions. Homogoneous Dirichlet
% Boundary Conditions are applied by eliminating corresponding rows and
% matrices of the element stiffness matrix. 
% -- Function takes in as the arguments the element stiffness matrix, and 6
% values which specify if a one of 6 rows and columns should be
% eliminated, with 1 if the corresponding row and column should be saved, and 0 if it
% should be eliminated. 
% -- NOTE: When applying the boundary conditions, the first three columns
%and rows correspond to the DOFs of the first node, and the other three to
%the second node. This comes in to play when building the global stiffness
%matrix. 

% -- TODO: Again, this is a very rigid code made specifically for a three
% by three matrix. Should try to make it more flexible later. 

reducedStiffMatrix = elementStiffnessMatrix;
redRow = 0;
if r1 == 0
    reducedStiffMatrix(1,:) =[];
    reducedStiffMatrix(:,1) = [];

    redRow = redRow +1;
end

if r2 == 0
    reducedStiffMatrix(2-redRow, :) = [];
    reducedStiffMatrix (:, 2-redRow) = [];

    redRow = redRow +1;
end

if r3 == 0
    reducedStiffMatrix(3-redRow, :) = [];
    reducedStiffMatrix (:, 3-redRow) = [];

    redRow = redRow +1;
end
if r4 == 0
    reducedStiffMatrix(4-redRow, :) = [];
    reducedStiffMatrix(:, 4-redRow) = [];

    redRow = redRow +1;
end
if r5 == 0
    reducedStiffMatrix(5-redRow, :) = [];
    reducedStiffMatrix (:, 5-redRow) = [];


    redRow = redRow +1;
end
if r6 == 0
    reducedStiffMatrix(6-redRow, :) = [];
    reducedStiffMatrix (:, 6-redRow) = [];

end
end

