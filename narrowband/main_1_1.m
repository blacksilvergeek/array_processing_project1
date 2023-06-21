clc, clear
close all
utilsFolder = fullfile(pwd, 'utils');
addpath(utilsFolder);




%% Add the folder to the MATLAB search path
% Get the full path of the "utils" folder
utilsFolder = fullfile(pwd, 'utils');
addpath(utilsFolder);

%% siganl model

mode = "original";
set_parameter % set parameter M, N, Delta, theta, f, SNR
is_noiseless = false;
[X,A,S] = gendata(M,N,Delta,theta,f,SNR,is_noiseless);
% Plot the singular values
X_ori = X;
% singular_plot(X)

% double the sample
mode = "double sample";
set_parameter
[X,A,S] = gendata(M,N,Delta,theta,f,SNR);
% Plot the singular values
figure
hold on
is_normalise = false;
singular_plot(X, is_normalise)
singular_plot(X_ori, is_normalise)
legend("double sample", "original")
grid on
hold off


% double the sample
mode = "double sample";
set_parameter
[X,A,S] = gendata(M,N,Delta,theta,f,SNR);
% Plot the singular values
figure
hold on
is_normalise = true;
singular_plot(X, is_normalise)
singular_plot(X_ori, is_normalise)
legend("double sample", "original")
grid on
hold off

% % extreme the sample
% mode = "extreme sample";
% set_parameter
% [X,A,S] = gendata(M,N,Delta,theta,f,SNR);
% % Plot the singular values
% figure
% hold on
% is_normalise = true;
% singular_plot(X, is_normalise)
% singular_plot(X_ori, is_normalise)
% legend("extreme", "original")
% hold off


mode = "double antenna";
set_parameter
[X,A,S] = gendata(M,N,Delta,theta,f,SNR);
% Plot the singular values
figure
hold on
is_normalise = false;
singular_plot(X, is_normalise)
singular_plot(X_ori, is_normalise)
legend("doubled number of antennas", "original")
grid minor
hold off

mode = "double antenna";
set_parameter
[X,A,S] = gendata(M,N,Delta,theta,f,SNR);
% Plot the singular values
figure
hold on
is_normalise = true;
singular_plot(X, is_normalise)
singular_plot(X_ori, is_normalise)
legend("doubled number of antennas", "original")
grid minor
hold off



mode = "small angle diff";
set_parameter
[X,A,S] = gendata(M,N,Delta,theta,f,SNR);
% Plot the singular values
figure
hold on
is_normalise = false;
singular_plot(X, is_normalise)
singular_plot(X_ori, is_normalise)
legend("small angle difference", "original")
grid minor
hold off

mode = "small angle diff";
set_parameter
[X,A,S] = gendata(M,N,Delta,theta,f,SNR);
% Plot the singular values
figure
hold on
is_normalise = true;
singular_plot(X, is_normalise)
singular_plot(X_ori, is_normalise)
legend("small angle difference", "original")
grid minor
hold off



mode = "small freq diff";
set_parameter
[X,A,S] = gendata(M,N,Delta,theta,f,SNR);
% Plot the singular values
figure
hold on
is_normalise = false;
singular_plot(X, is_normalise)
singular_plot(X_ori, is_normalise)
legend("small frequency difference", "original")
grid minor
hold off


mode = "small freq diff";
set_parameter
[X,A,S] = gendata(M,N,Delta,theta,f,SNR);
% Plot the singular values
figure
hold on
is_normalise = true;
singular_plot(X, is_normalise)
singular_plot(X_ori, is_normalise)
legend("small frequency difference", "original")
grid minor
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
    
    % stem
    stem(singular_values)
    title('Singular Values of X')
    xlabel('Index')
    ylabel('Singular Value')
end

function [] = all_close(data_est, data_true, data_name)
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
        else
            disp("estimation for "+ data_name + " is not accurate")
        end
    end

    % if data is real, compare the value
    if isreal(data_true)
        % reorder data by value
        data_true = sort(data_true);
        data_est = sort(data_est);
        if norm(data_true - data_est) < error_threshold
            disp("estimation for "+ data_name + " is accurate")
        else
            disp("estimation for "+ data_name + " is not accurate")
        end
    end



end