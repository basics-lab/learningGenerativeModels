function [U_avg,V_avg] = projSGD_multivar(Y, U0, V0)
    %% Parameters
    [d,n] = size(Y);
    lambda = 100;
    batch_size = 10;
    num_batch = floor(n/batch_size);
    U = U0;
    V = V0;
    r = 3;
    U_avg = zeros(d,d);
    V_avg = zeros(d,1);
    %% Iterations
    for t=1:num_batch
        alpha = 1/(lambda*t);
        tic;
        [dU, dV] = grad_lik_function(Y(:,(batch_size*(t-1)+1):batch_size*t),U,V);
        t1 = toc;
        V_old = V;
        V = V + alpha*dV;
        U_old = U;
        U = U+alpha*dU;
        tic;
        [L,D] = eig(U);
        D = diag(min(max(1/r, diag(D)), r));
        U = L*D*L';
        t2 = toc;
        try chol(U);
        catch
        keyboard;
        end
        U_avg = U_avg + U;
        V_avg = V_avg + V;
        fprintf("Grad Compute=%.3f, Projection=%.3f\n",t1,t2);
    end
    U_avg = U_avg/num_batch;
    V_avg = V_avg/num_batch;
end

