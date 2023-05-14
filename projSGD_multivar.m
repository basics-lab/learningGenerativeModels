function [U_avg,V_avg,U,V] = projSGD_multivar(Y, U0, V0)
    %% Parameters
    [d,n] = size(Y);
    lambda = 0.5;
    batch_size = 100;
    num_batch = floor(n/batch_size);
    U = U0;
    V = V0;
    r = 4;
    U_avg = zeros(d,d);
    V_avg = zeros(d,1);
    n_epochs = 10;
    n_updates = floor(num_batch/10);
    %% Iterations
    n_fail = 0;
    for i=1:n_epochs
        fprintf("Starting epoch %i of %i\n",i, n_epochs);
        fprintf("Total of %i Batches to be run\n", num_batch)
        for t=1:num_batch
            alpha = 1/(lambda*(t + (i-1)*num_batch));
            if mod(t,n_updates) == 0
                fprintf("Total of %d/%d Batches Complete\n", t,num_batch)
            end
            tic;
            [dU, dV, fail] = grad_lik_function(Y(:,(batch_size*(t-1)+1):batch_size*t),U,V);
            if fail == true
                n_fail = n_fail+1;
                continue;
            end
            dU = dU/batch_size;
            dV = dV/batch_size;
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
        end
    end
    U_avg = U_avg/(num_batch*n_epochs - n_fail);
    V_avg = V_avg/(num_batch*n_epochs - n_fail);
end