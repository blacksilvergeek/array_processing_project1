clear, clc
% Get the full path of the "utils" folder
utilsFolder = fullfile(pwd, 'utils');
close all
addpath(utilsFolder);

%%
mode = "comparison";
set_parameter
is_noiseless = false;
SNR =10;

M = 3;
[X,A,S] = gendata(M,N,Delta,theta,f,SNR);


f_hat = espritfreq(X,d);
% konw theta
% Construct new model based on estimation of frequency
[X, ~,S_real]=gendata(M, N, Delta, theta, f_hat, SNR);

W=pinv(X*pinv(S_real));
W = W';

% [~, A_est,~]=gendata(M, N, Delta, theta, f_hat, SNR);
% W=pinv(A_est);
% W=W';
% S_est=W'* X;

%% scan over -90 to 90 for freq estimation
d = 1;
ths = -90:1:90;
% ths = ths * pi / 180;
ys1 = [];
ys2 = [];
for th = ths
    [~, A_real,~]=gendata(M, N, Delta, th, f(1), SNR);
    y = W(:,1)' * A_real;
    ys1 = [ys1, y];
end

for th = ths
    [~, A_real,~]=gendata(M, N, Delta, th, f(2), SNR);
    y = W(:,2)' * A_real;
    ys2 = [ys2, y];
end

% plot ys
figure
plot(ths, abs(ys1))
hold on
plot(ths, abs(ys2))
hold off
legend('source 1', 'source 2')
title('spatial responses based on estimated frequency')

% run for many times to get stable result
run_times = 1000;
ys1_avg = zeros(1, length(ths));
ys2_avg = zeros(1, length(ths));
ys1s = zeros(run_times, length(ths));
ys2s = zeros(run_times, length(ths));

for i = 1:run_times
    [X,A,S] = gendata(M,N,Delta,theta,f,SNR);
    f_hat = espritfreq(X,d);
    % konw theta
    % Construct new model based on estimation of frequency
    [X, ~,S_real]=gendata(M, N, Delta, theta, f_hat, SNR);

    W=pinv(X*pinv(S_real));
    W = W';
    for k = 1:length(ths)
        [~, A_real,~]=gendata(M, N, Delta, ths(k), f(1), SNR);
        y = W(:,1)' * A_real;
        ys1_avg(k) = ys1_avg(k) + y;
    end
    for k = 1:length(ths)
        [~, A_real,~]=gendata(M, N, Delta, ths(k), f(2), SNR);
        y = W(:,2)' * A_real;
        ys2_avg(k) = ys2_avg(k) + y;
    end
end

ys1_avg = ys1_avg / run_times;
ys2_avg = ys2_avg / run_times;
% plot
figure
plot(ths, abs(ys1_avg))
hold on
plot(ths, abs(ys2_avg))
hold off
legend('source 1', 'source 2')
title('spatial responses based on estimated frequency')



%% angle
% know frequency
mode = "comparison";
set_parameter
is_noiseless = false;
SNR =10;

M = 3;
[X,A,S] = gendata(M,N,Delta,theta,f,SNR);


%% scan over -90 to 90 for angle estimation
th_hat = esprit(X,d)
[~, A_est,~]=gendata(M, N, Delta, th_hat, f, SNR);
W=pinv(A_est);
W = W';


d = 1;
ths = -90:1:90;
% ths = ths * pi / 180;
ys1 = [];
ys2 = [];
for th = ths
    [~, A_real,~]=gendata(M, N, Delta, th, f(1), SNR);
    y = W(:,1)' * A_real;
    ys1 = [ys1, y];
end
for th = ths
    [~, A_real,~]=gendata(M, N, Delta, th, f(2), SNR);
    y = W(:,2)' * A_real;
    ys2 = [ys2, y];
end
% plot ys
figure
plot(ths, abs(ys1))
hold on
plot(ths, abs(ys2))
hold off
legend('source 1', 'source 2')
title('spatial responses based on estimated angle')
