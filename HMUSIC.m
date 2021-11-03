function f0 = HMUSIC(R, N, pitchBounds, L, K)
%--------------------------------------------------------------------------
%   Fundamental frequency(Pitch) Estimation with harmonics MUSIC 
%
%   Usage:
%       f0 = HMUSIC(R, N, pitchBounds, L, K)
%       f0 = HMUSIC(R, N, pitchBounds, L)
%   Output:
%       f0: estimated pitch
%
%   Input:
%       R: covariance matrix of size M
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
%       Multi pitch estimation
%       SUBSPACE-BASED FUNAMENTAL FREQUENCY ESTIMATION *
%   All copyrights reserved, 10-28, 2021.
%--------------------------------------------------------------------------
if pitchBounds(1)<1/N||pitchBounds(2)>=1
    error('pitch should in range of  1/N~1 ');
end
fRange = pitchBounds(1):1/N:pitchBounds(2);
fLen = length(fRange);
P = zeros(size(fRange));
M = size(R, 1);
%# eq. 10 EVD with eigenvalue in descent order
[EV,D]=eig(R);%
EVA=diag(D)';
[EVA,I]=sort(EVA);
EVA=fliplr(EVA);
EV=fliplr(EV(:,I));
for iterNum = 1:fLen
    fIter = fRange(iterNum);
    if nargin<5
        K=1;
    end
    %# eq. 5 construct the Vandermonde Matrix
    Z = zeros(M,L);
    m_lin = (0:M-1)';
    z = exp(1j*2*pi*fIter*m_lin);
    for lIter = 1:L
        Z(:,lIter) = z.^lIter;      
    end
    KL = K*L;
    %# eq. 11 construct noise subspace
    G = EV(:,KL+1:M);
    %# eq. 18 calculate pesudo-spectrum
    denom = norm(Z'*G, 'fro').^2;
    P(iterNum) = abs(1/denom);
end
P = P./max(P);
[~, pos] = max(P);
f0 = fRange(pos);
figure;
plot(fRange,P);
legend('HMUSIC');
%-------------------------------EOF----------------------------------------