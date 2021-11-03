function R = CoMat_estimation(x, L, mode, lambda)
%--------------------------------------------------------------------------
%% estimate covariance matrix with observed signal
% Usage:
%       R = CoMat_estimation(x, L, mode, lambda)
% Inputs:
%       x: input signal of size 1*L
%       L: dimension of estimated cvariance matrix
%       mode: 
%           forward:  forward estimation
%           backward: backward estimation
%           forward-backward: forward-backward(modified) method
%           recursive: using forgetting factor to update estimation
%       lambda: forgetting factor
% Output:
%       R: estimated covariance matrix of size L*L
% Author: 
%       Xianrui Wang, Center of Intelligent Acoustics and Immersive
%       Communications (CIAIC).
% Contact:     
%       wangxianrui@mail.nwpu.edu.cn
% Date:
%       6-16,2021
% Reference:
%       Performance analysis of forward-backward
%       matched-filterbank spectral estimators
% all conpyrights preserved
%--------------------------------------------------------------------------
if nargin<3
    error('signal, matrix dimension and mode necessary, lambda optional');
end
N = length(x);
reshape(x, N, 1);
%-----------------------------use forward estimation-----------------------
if mode == "forward"
    R = zeros(L);
    x_f = x;
    for i = 1:N-L+1
        xf = x_f(i+L-1:-1:i); 
        %# forward estimation
        R = R+xf*xf';
    end
    R = R/(N-L+1);
elseif mode == "backward"
    R = zeros(L);
    x_b = conj(flip(x));               
    for i = 1:N-L+1
        xb = x_b(i+L-1:-1:i); 
        %# backward estimation
        R = R+xb*xb';
    end   
    R = R/(N-L+1);
elseif mode == "modified"   
    R = zeros(L);                      
    x_f = x;
    x_b = conj(flip(x));  
    for i = 1:N-L+1
        xf = x_f(i+L-1:-1:i);         
        xb = x_b(i+L-1:-1:i);        
        %# forward-backward estimation
        R = R + (xf*xf'+xb*xb')/2;
    end
    R = R/(N-L+1);
elseif mode == "recursive"
    if nargin==3
        lambda = 0.95;
    end
    R = 1e-2*eye(L);
    xf = zeros(L,1);
    xb = zeros(L,1);
    for i = 1:N
        xf = [x(i);xf(1:end-1)];
        xb = [xb(2:end);x(i)'];
        Rt = (xf*xf'+xb*xb')/2;
        R = lambda*R + (1-lambda)*Rt;
    end
end
%----------------------------------EOF-------------------------------------