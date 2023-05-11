function psd = isPSD(A)
    try chol(A);
        psd = true;
    catch ME
        psd = false;
    end