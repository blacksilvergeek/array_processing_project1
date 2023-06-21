function [theta,f] = joint(X,d,m)
    % ESPRIT algorithm implementation
    % Inputs:
    %   X: Observation signal matrix (size: M x N)
    %   d: Number of sources
    % Outputs:
    %   theta: Estimated arrival angles of the sources (size: d x 1)
    %   f:     Estimated frequency of the sources (size: d x 1)


    % load Delta, d from workspace
    Delta = evalin('base', 'Delta');

    signal = X;
    [M,N] = size(signal);
    M_re = M-1;
    jthresh = 1e-6;
    % following the notation in the course notes
    X = [signal(1:end-1,1:end-1)]; % X
    Y = [signal(1:end-1,2:end)]; % PHI time
    Z = [signal(2:end,1:end-1)]; % TH
    
    % get size and rescale size


    % stack the matrices
    ALL = [X;Y;Z];
    
    [U,S,V] = svd(ALL);
    
    % unpack the matrices of singular vectors
    % since only 2 sources, we only need the first 2 columns of U
    Ux = U(1:M_re,1:d);
    Uy = U((1*M_re+1):2*M_re,1:d);
    Uz = U((2*M_re+1):end,1:d);
    
    % result of LS
    obj_y = pinv(Ux)*Uy;
    obj_z = pinv(Ux)*Uz;
    
    % joint diagonalization
    [Vj,Dj]=joint_diag([obj_y, obj_z], jthresh);

    % get the corresponding f and theta
    D_j_y = diag(Dj(1:d,1:d));
    D_j_z = diag(Dj(1:d,d+1:end));
    

    f_hat = angle(D_j_y)/(2*pi);
    th_hat = asin(angle(D_j_z)/(2*pi*Delta)) * 180 / pi;

    theta = th_hat;
    f = f_hat;
end