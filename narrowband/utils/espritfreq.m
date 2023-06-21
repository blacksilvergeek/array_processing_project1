function f = espritfreq(X,d)
    signal = X;

    % LOAD M from workspace
    [M, N] = size(signal);


    X = [signal(:,1:end-1); signal(:,2:end)];
    
    % svd of X
    [U, S, V] = svd(X);
    U = U(:,1:2);
    
    Ux = U(1:end-M,:);
    Uy = U((M+1):end,:);
    
    obj_opt = pinv(Ux)*Uy;
    
    % eig of A
    [V, D] = eig(obj_opt);
    f=angle(diag(D))/(2*pi);
end 




