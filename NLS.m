function f0 = NLS(x, N, pitchBounds, L)
%--------------------------------------------------------------------------
%   Fundamental frequency(Pitch) Estimation with
%   Nonlinear Least Squares(NLS) 
%
%   Usage:
%       f0 = NLS(R, N, pitchBounds, L)
%   Output:
%       f0: estimated pitch
%
%   Input:
%       x: input data
%       N: number of searching grid, (uniformly searching)
%       L: maximum model number of harmonics (i.e., order) that is expected
%          (positive integer, scalar)
%       pitchBounds: Lower and upper bounds on the fundamental frequency in
%          cycles/sample. The lower bound should not be set lower than
%          1/N and the upper bound can at most be 1 
%
%   Author:
%       Xianrui Wang, Center of Intelligent Acoustics and Immersive
%       Communications.
%
%   Contact:
%       wangxianrui@mail.nwpu.edu.cn
%   Reference:
%       Multi Pitch Estimation
%   This manuscript strictly follows the procedure in the paper 
%   All copyrights reserved, 11-2, 2021.
%--------------------------------------------------------------------------
if pitchBounds(1)<1/N||pitchBounds(2)>=1
    error('pitch should in range of  1/N~1 ');
end
fRange = pitchBounds(1):1/N:pitchBounds(2);
fLen = length(fRange);
P = zeros(size(fRange));
x = reshape(x, [], 1);
M = size(x, 1);
for iterNum = 1:fLen
    fIter = fRange(iterNum);
    %# eq. 5 construct the Vandermonde Matrix
    Z = zeros(M,L);
    m_lin = (0:M-1)';
    z = exp(-1j*2*pi*fIter*m_lin);
    for lIter = 1:L
        Z(:,lIter) = z.^lIter;      
    end
    % eq. 13 calculate spectrum
    P(iterNum) = norm(Z'*x)^2;
end
P = P./max(P);
[~, pos] = max(P);
f0 = fRange(pos);
figure;
plot(fRange,P);
legend('NLS');
%-------------------------------EOF----------------------------------------