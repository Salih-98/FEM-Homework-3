function SM = BernoulliElementStiffnessMatrix(EA, EI, L)
SM1 = [EA/L 0 0 -EA/L 0 0]; 
SM2 = [0 12*EI/L^3 -6*EI/L^2 0 -12*EI/L^3 -6*EI/L^2];
SM3 = [0 -6*EI/L^2 4*EI/L 0 6*EI/L^2 2*EI/L];
SM4 = [-EA/L 0 0 EA/L 0 0];
SM5 = [0 -12*EI/L^3 6*EI/L^2 0 12*EI/L^3 6*EI/L^2];
SM6 = [0 -6*EI/L^2 2*EI/L 0 6*EI/L^2 4*EI/L];
SM = [SM1; SM2; SM3; SM4; SM5; SM6];
end

