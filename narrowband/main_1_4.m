clc, clear
close all



%% Add the folder to the MATLAB search path
% Get the full path of the "utils" folder
utilsFolder = fullfile(pwd, 'utils');
addpath(utilsFolder);


%% joint 
mode = "original";
m = 5;

is_noiseless = true;
set_parameter
[X,A,S] = gendata(M,N,Delta,theta,f,SNR,is_noiseless);
[th_hat,f_hat] = joint(X,d,m);
all_close(f_hat, f, "frequency of joint");
all_close(th_hat, theta, "directions of joint");


%% Beamformer
%% Zero-forcing beamformers based on the direction estimation
% generate data
[X,A,S] = gendata(M,N,Delta,theta,f,SNR);
angle_hat=esprit(X,d);
% Construct new model based on estimation of angles
[~, A_est, ~]=gendata(M, N, Delta, angle_hat, f, SNR);
% Zeros-forcing beamformer without noise
W=pinv(A_est);
W=W';
S_est=W'* X;

%%% test noiseless accuracy
is_noiseless = true;
[X,A,S] = gendata(M,N,Delta,theta,f,SNR,is_noiseless);
angle_hat=esprit(X,d);
% Construct new model based on estimation of angles
[~, A_est,~]=gendata(M, N, Delta, angle_hat, f, SNR);
% Zeros-forcing beamformer without noise
W=pinv(A_est);
W=W';
S_est=W'* X;
all_close(S_est, S, "beamformer based on direction estimation");




%% Zero-forcing beamformers based on the frequency estimation
% generate data
[X,A,S] = gendata(M,N,Delta,theta,f,SNR);

frequency_hat = espritfreq(X,d);
% Construct new model based on estimation of frequency
[~, A_est,~]=gendata(M, N, Delta, theta, frequency_hat, SNR);
W=pinv(A_est);
W=W';
S_est=W'* X;

%%% test noiseless accuracy
is_noiseless = true;
[X,A,S] = gendata(M,N,Delta,theta,f,SNR,is_noiseless);
frequency_hat = espritfreq(X,d);
% Construct new model based on estimation of frequency
[~, A_est,~ ]=gendata(M, N, Delta, theta, frequency_hat, SNR);
W=pinv(A_est);
W=W';
S_est=W'* X;
all_close(S_est, S, "beamformer based on frequency estimation");

%% plot for different m
ms = 1:10;
error_fs = zeros(1, length(ms));
error_ths = zeros(1, length(ms));
for k = 1:length(ms)
    m = ms(k);
    [X,A,S] = gendata(M,N,Delta,theta,f,SNR);
    [th_hat,f_hat] = joint(X,d,m);
    error_f = all_close(f_hat, f, "frequency of joint");
    error_th = all_close(th_hat, theta, "directions of joint");
    error_fs(m) = error_f;
    error_ths(m) = error_th;
end
% plot
figure
plot(ms, error_fs, 'r', 'LineWidth', 2)
hold on
plot(ms, error_ths, 'b', 'LineWidth', 2)
legend('frequency error', 'direction error')
xlabel('m')
ylabel('error (dB)')
title('error vs m')
hold off





%% function space in this file
function  [] = singular_plot(X, is_normalise)
    if nargin < 2 || isempty(is_normalise)
        is_normalise = false;
    end
    singular_values = svd(X);
    if is_normalise
        singular_values = singular_values / sum(singular_values);
    end
    % Plot the singular values
    figure
    % stem
    stem(singular_values)
    title('Singular Values of X')
    xlabel('Index')
    ylabel('Singular Value')
end

function [error_db] = all_close(data_est, data_true, data_name)
    % actually, data_est and data_true can be in any order
    % so we need to sort them first
    
    % But simple ordering doesn't work if data are complex and norm 1
    % so we'll tackle these cases (data real or not) separately

    % if data is complex, compare the arg value
    error_threshold = 1e-6;
    if ~isreal(data_true)
        % get arg value
        true_arg = angle(data_true);
        est_arg = angle(data_est);

        % reorder data by arg value
        [~, true_arg_order] = sort(true_arg);
        [~, est_arg_order] = sort(est_arg);
        data_true = data_true(true_arg_order);
        data_est = data_est(est_arg_order);
        if norm(data_true - data_est) < error_threshold
            disp("estimation for "+ data_name + " is accurate")
            % print error in dB
            disp("error in dB is " + num2str(20*log10(norm(data_true - data_est))))
        else
            disp("estimation for "+ data_name + " is not accurate")
            disp("error in dB is " + num2str(20*log10(norm(data_true - data_est))))
        end
    end

    % if data is real, compare the value
    if isreal(data_true)
        % reorder data by value
        data_true = sort(data_true);
        data_est = sort(data_est);
        if norm(data_true - data_est) < error_threshold
            disp("estimation for "+ data_name + " is accurate")
            disp("error in dB is " + num2str(20*log10(norm(data_true - data_est))))
        else
            disp("estimation for "+ data_name + " is not accurate")
            disp("error in dB is " + num2str(20*log10(norm(data_true - data_est))))
        end
    end

    error_db = 20*log10(norm(data_true - data_est));

end