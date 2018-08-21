function maps = tradESPIRiT(X, r, k, c, wnsvn)
    % [maps, M, W] = tradESPIRiT(X, r, k, c, wnsvn)
    %
    % Function that calculates traditional ESPIRiT maps. 
    %
    % INPUTS:
    %	X         - k-space data [kx, ky, nc]
    %	r         - size of calibration region
    %	k         - kernel size
    %	c         - crop threshold
    %	wnsvn     - window normalized number of singular values (determines subspace size)
    % 
    % OUTPUTS
    %   maps      - espirit sensitivity maps [kx, ky, nc, nc]
    %   MSE       - estimated MSE given parameters choices [lst_k dim, lst_c dim, lst_wnsvn dim]
    %   optParam  - optimum param choices given SURE metric [opt_k, opt_c, opt_wnsnv]
    %
    % Copyright 2016. The Regents of the University of California.

    nc = size(X, 3);

    C = extractCalreg(X, r); 

    imSize = [size(X, 1), size(X, 2)];

    A = calreg2calmat(C, k); [U, S, V] = svd(A); s = diag(S);

    n = floor(wnsvn * k * k);
    vecWeight = ones([n, 1]);
    vecWeight = [vecWeight; zeros([size(V, 2) - n, 1])];

    wV = V * diag(vecWeight);

    kernel = reshape(wV, k, k, nc, size(wV, 2));
    [M, W] = kernelEig(kernel, [size(X, 1), size(X, 2)]);

    maps = cropEigvec(M, W, c);

end
