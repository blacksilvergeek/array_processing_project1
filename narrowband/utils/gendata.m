function [X, A, S] = gendata(M, N, Delta, theta_in_angle, f, SNR, is_noiseless)
% Generate array signal data
% Inputs:
%   M: Number of array elements
%   N: Signal length
%   Delta: Inter-element spacing
%   theta: Arrival angles of signals (in degrees)
%   f: Frequencies of signals (in Hz)
%   SNR: Signal-to-Noise Ratio (in dB)
% Outputs:
%   X: Observation signal matrix (size: M x N)
%   A: Array response matrix (size: M x d)
%   S: Source signal matrix (size: d x N)

if nargin < 7 || isempty(is_noiseless)
    is_noiseless = false;
end

% Convert angles to radians
theta = theta_in_angle / 180 * pi;

% Number of sources
d = length(theta);

% Generate array response matrix A
A = zeros(M, d);
for m = 1:d
    angle = sin(theta(m));
    for i = 1:M
        A(i, m) = exp(1j * 2 * pi * Delta * angle * (i - 1));
    end
end

% Generate source signal matrix S
% This is given by instruction
S = zeros(d, N);
for m = 1:d

    freq = f(m);

    for k = 1:N
        S(m, k) = exp(1j * 2 * pi * freq * k);
    end
end

% Add noise to the observation signal matrix X
n = randn(M, N) + 1j * randn(M, N);
noise_power = norm(n, 'fro')^2 / (M * N);
signal_power = norm(S, 'fro')^2 / (M * N);

% SNR corresponds num 
scale_SNR = 10^(SNR / 10);

% calculate the scaling factor for noise
scale = sqrt(signal_power / (scale_SNR * noise_power));

% correct the noise
n_corrected = scale * n;
% rank(S)
% rank(A*S)
if is_noiseless
    X = A * S;
else
    X = A * S + n_corrected;
end
% X = A * S;
% X = A * S + n_corrected;

end





