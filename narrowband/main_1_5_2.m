clc, clear
close all
m = 5;
mode = "comparison";
set_parameter
is_noiseless = true;

%% Beamformer
%% Zero-forcing beamformers based on the direction estimation
% generate data


%%% test noiseless accuracy
is_noiseless = true;
[X,A,S] = gendata(M,N,Delta,theta,f,SNR, is_noiseless);

th_hat = esprit(X,d)
[~, A_est,~]=gendata(M, N, Delta, th_hat, f, SNR, is_noiseless);
W=pinv(A_est);
W = W';

S_est = W'*X;
all_close(S_est, S, "beamformer based on direction estimation");




%% Zero-forcing beamformers based on the frequency estimation
% generate data


%%% test noiseless accuracy
[X,A,S] = gendata(M,N,Delta,theta,f,SNR, is_noiseless);


f_hat = espritfreq(X,d);
% konw theta
% Construct new model based on estimation of frequency
[X, ~,S_real]=gendata(M, N, Delta, theta, f_hat, SNR, is_noiseless);

W=pinv(X*pinv(S_real));
W = W';
S_est = W'*X;
all_close(S_est, S, "beamformer based on frequency estimation");






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

    % if data is complex, check corr



    if ~isreal(data_true)
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
        % get arg value

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