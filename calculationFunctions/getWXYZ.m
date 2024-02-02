function [W,X,Y,Z] = getWXYZ(iNum)
% -- getWXYZ function - extracts values for task parameters from the
% immatrikulation nummer

i = length(iNum);
Z = str2double(iNum(i)); 
Y = str2double(iNum(i-1));
X = str2double(iNum(i-2));
W = str2double(iNum(i-3));
end

