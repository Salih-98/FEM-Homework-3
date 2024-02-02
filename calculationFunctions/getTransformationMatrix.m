function TM = getTransformationMatrix(C1,C2)
TM = [C1 C2 0 0 0 0; -C2 C1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 C1 C2 0;...
    0 0 0 -C2 C1 0; 0 0 0 0 0 1];
end

