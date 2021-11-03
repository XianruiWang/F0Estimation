function invR_robust = Func_inv(X, method, threshold)
%--------------------------------------------------------------------------
% calculate robust inverse of the input matrix X
% Usage:
%       invR_robust = Func_inv(X, method, threshold)
% Inputs:
%       X: input Matric
%       method: 
%           "diag": diagnol loading method
%           "svd":  singular value decomposition method
%           "eig":  eigenvalue decomposition method
%       threshold: threshold value 
% Output:
%       invR_robust: robust inverse of X
% Author: 
%       Gongping Huang, Center of Intelligent Acoustics and Immersive
%       Communications (CIAIC), 2013.
%       Modified by Xianrui Wang, 2021
% Contact:     
%       wangxianrui@mail.nwpu.edu.cn
%
% All conpyrights preserved
%--------------------------------------------------------------------------
if size(X,1)~=size(X,2)
    error("Input matrix must be square");
end
dim = size(X,1);
if nargin==2
    threshold = 1e-8;
end
%--------------------------------------------------------------------------
if method == "diag"
    X = X + threshold*trace(X)*eye(dim)/dim;
    invR_robust = inv(X);
elseif method == "svd"
    [U,S,V] = svd(X);
    s = diag(S);
    threshold = dim*eps(max(s));
    r = sum(s>threshold);
    S_inv = diag*(ones(r,1)./s(1:r));
    invR_robust = V(:, 1:r) * S_inv * U(:, 1:r)';
elseif method =="eig"
    [Q,A] = eig(X);
    maxdiagv = max(diag(A));
    threshold = min(maxdiagv*1.e-3, threshold);
    AA = diag(A);
    ipos = find(AA>threshold);
    AA(find(AA<=threshold)) = 0;
    AA(ipos) = 1./AA(ipos);
    A = diag(AA);
    invR_robust = Q * A * Q';
else
    error("unsupported method");
end
%--------------------------------------EOF---------------------------------
    