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

    % all stack for time 
    all_stack_t = [];
    for k=1:m
        % assign time shift of X is Xk
        Yk = signal(1:end,k:end-(m-k));
        all_stack_t = [all_stack_t; Yk];
    end
    
    
    [U,S,V] = svd(all_stack_t);
    
    % unpack the matrices of singular vectors
    % since only 2 sources, we only need the first 2 columns of U
    Ux_f = U(1:M*(m-1),1:d);
    % This is actually equal to the kron product method
    % Just a different way of presenting it
    % J_x = [eye(m-1), zeros(m-1,1)];
    % J_x = kron(J_x, eye(M));
    % Ux_f2 = J_x*U(:,1:d);
    % you can test actually Ux_f and Ux_f2 are the same

    Uy_f = U(M+1:end,1:d);
    
    % take all first M-1 rows of every M rows of U
    Ux_th = [];
    for k=1:m
        Ux_th = [Ux_th; U((k-1)*M+1:(k-1)*M+M-1,1:d)];
    end
    Uz_th = [];
    for k=1:m
        Uz_th = [Uz_th; U((k-1)*M+2:(k-1)*M+M,1:d)];
    end

    % result of LS
    obj_y = pinv(Ux_f)*Uy_f;
    obj_z = pinv(Ux_th)*Uz_th;
    
    % joint diagonalization
    [Vj,Dj]=joint_diag([obj_y, obj_z], jthresh);

    % get the corresponding f and theta
    % here D_j_y and D_j_z are diag vectors for simplicity
    D_j_y = diag(Dj(1:d,1:d)); 
    D_j_z = diag(Dj(1:d,d+1:end));
    

    f_hat = angle(D_j_y)/(2*pi);
    th_hat = asin(angle(D_j_z)/(2*pi*Delta)) * 180 / pi;

    theta = th_hat;
    f = f_hat;
end