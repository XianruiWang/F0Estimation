clear;clc;
%--------------------------------------------------------------------------
%% construct observed signal consists of five sinusoids,
% corrupted by complex white gaussian noise
% location of spectral lines     
fVect = [0.16, 0.32, 0.48];   
% amplitude of frequency components
coefVec = [1, 1, 1, 1, 1];
% number of frequency components
nf = length(fVect);    
% angular frequency
wVect = 2*pi*fVect;     
% number of available samples
T = 256;
nT = 0:T-1;
x = randn(T,1);
for fIndex = 1:nf
    %# sum all components and create signals
    % random phase
    phi = 2*pi*(2*rand(1,1)-1);     
    %# sum all components
    x = x + 2*coefVec(fIndex)*cos(wVect(fIndex)*nT+phi).';   
end
%# ensure the signal is real
x = real(x);
x = x-mean(x);
%--------------------------------------------------------------------------
%% estimate covariance matrix with modified method
M=128;
N=10000;
% the range of pitch, should be carefully choosen to avoid 
% doubling and halving
pitchBounds=[0.01, 0.4];
R = CoMat_estimation(x, M, 'forward');
f0_MUISC = HMUSIC(R, N, pitchBounds,3);
f0_NLS = NLS(x, N, pitchBounds,3);
f0_Capon = Capon(R, N, pitchBounds,3);