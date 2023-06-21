function theta = esprit(X, d)
    % ESPRIT algorithm implementation
    % Inputs:
    %   X: Observation signal matrix (size: M x N)
    %   d: Number of sources
    % Outputs:
    %   theta: Estimated arrival angles of the sources (size: d x 1)
    
    % to give same notation as in the slides
    % load Delta, d from workspace
    Delta = evalin('base', 'Delta');
    

    Z = X;
    
    % expand the matrix Z
    X = Z(1:end-1, :);
    Y = Z(2:end, :);

    Z_expand = [X; Y];


    [U_hat, S_hat, V_hat] = svd(Z_expand);
    
    U_hat = U_hat(:, 1:d);

    U_x_hat = U_hat(1:end/2, :);
    U_y_hat = U_hat(end/2+1:end, :);

    
    % pesudo-inverse of U_x_hat
    U_x_hat_inv = (U_x_hat' * U_x_hat)^-1 * U_x_hat';

    obj_est = U_x_hat_inv * U_y_hat;
    [T_inv, TH] = eig(obj_est);

    % diagonalize the TH
    TH = diag(TH);

    % get the arg value of the eigenvalues
    angles = angle((TH)); % this is 2 pi * Delta * sin(theta)

    % get the estimated angles theta
    theta = asin(angles / (2 * pi * Delta));
    % rad to angle
    theta = theta * 180 / pi;
    


end