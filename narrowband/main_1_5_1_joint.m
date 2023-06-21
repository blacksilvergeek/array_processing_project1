clc, clear
close all
m = 5;


%% Add the folder to the MATLAB search path
% Get the full path of the "utils" folder
utilsFolder = fullfile(pwd, 'utils');
addpath(utilsFolder);


%% estimation of directions
mode = "comparison";
set_parameter
algo = "joint";
dbs = 0:2:20;
th_joint_averages = zeros(2, length(dbs));
th_joint_deviations = zeros(2, length(dbs));
f_joint_averages = zeros(2, length(dbs));
f_joint_deviations = zeros(2, length(dbs));



for k = 1:length(dbs)
    % set SNR
    SNR = dbs(k);
    results = mean_and_deviation(SNR, algo);
    th_joint_averages(:,k) = results(:,1);
    th_joint_deviations(:,k) = results(:,2);
    f_joint_averages(:,k) = results(:,3);
    f_joint_deviations(:,k) = results(:,4);
    
end


% plot the average and standard deviation
errorbar(dbs, f_joint_averages(1,:), f_joint_deviations(1,:), 'o-', 'LineWidth', 1.5);
hold on
errorbar(dbs, f_joint_averages(2,:), f_joint_deviations(2,:), 'o-', 'LineWidth', 1.5);
legend('frequency 1', 'frequency 2');


title('Mean and Standard Deviation');
xlabel('SNR (dB)');
ylabel('Normalized Frequency (Hz)');
hold off

figure
errorbar(dbs, th_joint_averages(1,:), th_joint_deviations(1,:), 'o-', 'LineWidth', 1.5);
hold on
errorbar(dbs, th_joint_averages(2,:), th_joint_deviations(2,:), 'o-', 'LineWidth', 1.5);
legend('direction 1', 'direction 2');
title('Mean and Standard Deviation');
xlabel('SNR (dB)');
ylabel('Estimation of Direction (degree)');




function results = mean_and_deviation(SNR,algo)
    % run for 1000 times to find the average error and deviation
    % load the parameters from workspace
    M = evalin('base', 'M');
    N = evalin('base', 'N');
    Delta = evalin('base', 'Delta');
    theta = evalin('base', 'theta');
    f = evalin('base', 'f');
    d = evalin('base', 'd');
    m = evalin('base', 'm');
    num_of_runs = 1000;
    if algo == "esprit"
        th_esprits = zeros(2, num_of_runs);
        for i = 1:num_of_runs
            [X,A,S] = gendata(M,N,Delta,theta,f,SNR);
            f_hat = esprit(X,d);
            % sort the estimation
            f_hat = sort(f_hat);
            th_esprits(:,i) = f_hat;
            % disp("run " + num2str(i) + " finished")
        end
        th_esprit_average = mean(th_esprits, 2);
        th_esprit_deviation = std(th_esprits, 0, 2);
        results = [th_esprit_average, th_esprit_deviation];
    end
    if algo == "espritfreq"
        f_esprits = zeros(2, num_of_runs);
        for i = 1:num_of_runs
            [X,A,S] = gendata(M,N,Delta,theta,f,SNR);
            f_hat = espritfreq(X,d);
            % sort the estimation
            f_hat = sort(f_hat);
            f_esprits(:,i) = f_hat;
        end
        f_esprit_average = mean(f_esprits, 2);
        f_esprit_deviation = std(f_esprits, 0, 2);
        results = [f_esprit_average, f_esprit_deviation];
    end
    if algo == "joint"
        th_joint = zeros(2, num_of_runs);
        f_joint = zeros(2, num_of_runs);
        for i = 1:num_of_runs
            [X,A,S] = gendata(M,N,Delta,theta,f,SNR);
            [th_hat,f_hat] = joint(X,d,m);
            % sort the estimation
            th_hat = sort(th_hat);
            f_hat = sort(f_hat);
            th_joint(:,i) = th_hat;
            f_joint(:,i) = f_hat;
            % disp("run " + num2str(i) + " finished")
        end
        th_joint_average = mean(th_joint, 2);
        th_joint_deviation = std(th_joint, 0, 2);
        f_joint_average = mean(f_joint, 2);
        f_joint_deviation = std(f_joint, 0, 2);
        results = [th_joint_average, th_joint_deviation, f_joint_average, f_joint_deviation];
    end

end




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
    corr = data_true*data_est';
    % diag sum
    diag_sum = real(sum(diag(corr)));
    % counter diag sum
    counter_diag_sum = real(sum(diag(fliplr(corr))));
    % if diag_sum smaller, then reverse the row of data_est
    if diag_sum < counter_diag_sum
        data_est = flipud(data_est);
    end


    if ~isreal(data_true)
        % get arg value
        true_arg = angle(data_true);
        est_arg = angle(data_est);

        % reorder data by arg value
        % [~, true_arg_order] = sort(true_arg);
        % [~, est_arg_order] = sort(est_arg);
        % data_true = data_true(true_arg_order);
        % data_est = data_est(est_arg_order);
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