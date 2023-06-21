%   d: Number of signal sources
%   M: Number of array elements (antennas)
%   N: Signal length
%   Delta: Inter-element spacing
%   theta: Arrival angles of signals (in degrees)
%   f: Frequencies of signals (in Hz)
%   SNR: Signal-to-Noise Ratio (in dB)

% if mode does not exist, mode = "original"
if ~exist('mode','var')
    mode = "original";
end


% default values
d=2;
M=5;
N=20;
theta=[-20;30];
f=[0.1;0.3];
SNR=20;
Delta=1/2;



if mode == "original"
    ; % do nothing
elseif mode == "double sample"
    N=2*N; % double the number of samples
elseif mode == "extreme sample"
    N=100*N; % extreme the number of samples
elseif mode == "double antenna"
    M=2*M; % double the number of antennas
elseif mode == "small angle diff"
    theta=[10;20];
elseif mode == "small freq diff"
    f=[0.1;0.11];
elseif mode == "comparison"
    % d = 2 sources, M = 3, N = 20, θ = [−20, 30]T , f = [0.1, 0.12]T
    d=2;
    M=3;
    N=20;
    theta=[-20;30];
    f=[0.1;0.12];
end