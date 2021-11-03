function [V, D] = descenteig(R)
% DESCENTEIG: eigenvalue decomposition of the matrix R,
% with the eigenvaule listed in a descent order
% Usage: [V, D] = descenteig(R)
%   Input:
%     R : covariance matrix
%   Output:
%     D : Eigenvalue matrix
%     V: eigenvector matrix
%
% ------
% Author:  Jigndong Chen, Ph.D.
% Copyright 2012(c)
%
if (nargin ~= 1)
    error('Incorrect number of inputs!');
end
L = size(R, 1);
[HV, HD] = eig(R);
for i = 1 : L
    for j = i: L
        if HD(i,i) < HD(j,j)
            tempD = HD(i,i);
            HD(i,i) = HD(j,j);
            HD(j,j) = tempD;
            tempV = HV(:,i);
            HV(:,i) = HV(:,j);
            HV(:,j) = tempV;
        end
    end
end
if nargout == 1
    D = HD;
end
if nargout == 2
    D = HD;
    V = HV;
end
