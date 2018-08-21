function [maps, MSE, optParam] = autoESPIRiT(X, r, weight, stdev, lst_k, lst_c, lst_wnsvn)
    % [maps, MSE] = autoESPIRiT(X, r, weight, stdev, lst_k, lst_c, [lst_wnsvn )
    %
    % Function that calculates ESPIRiT maps. Uses Stein's Unbiased Risk Estimate (or SURE) 
    % as a metric to determine optimal parameters from the passed in arguments.
    %
    % INPUTS:
    %	X         - k-space data [kx, ky, nc]
    %	r         - size of calibration region
    %	weight    - set to true to soft-weight subspace. set to false to sweep subspace sizes.
    %	stdev     - estimated noise standard deviation
    %	lst_k     - list of kernel sizes to choose from [k1, k2, ...]
    %	lst_c     - list of crop thresholds to choose from [c1, c2, ...]
    %	lst_wnsvn - list of window normalized number of singular values to choose from [w1, w2, ...]
    %               required if 'weight' is set to true
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

    fprintf('Parameters                                                  | MSE\n');
    if (weight)
        [MSE, minMSE, optParam] = estimateMSE(C, r, imSize, weight, stdev, lst_k, lst_c);
        fprintf('\nOptimal Parameters:\n   k: %d\n   c: %f\n', optParam(1), optParam(2));
    else
        [MSE, minMSE, optParam] = estimateMSE(C, r, imSize, weight, stdev, lst_k, lst_c, lst_wnsvn);
	      fprintf('\nOptimal Parameters:\n   k: %d\n   c: %f\n   wnsvn: %f\n', optParam(1), optParam(2), optParam(3));
    end

    k = optParam(1);
    c = optParam(2);

    A = calreg2calmat(C, k); [U, S, V] = svd(A); s = diag(S);

    if (weight == true)
        vecWeight = svWeightsSURE(s, stdev, size(A), false);
        vecWeight = [vecWeight; zeros([size(V, 2) - length(vecWeight), 1])];
    else
	      n = floor(optParam(3) * k * k);
	      vecWeight = ones([n, 1]);
	      vecWeight = [vecWeight; zeros([size(V, 2) - n, 1])];
    end

    wV = V * diag(vecWeight);

    kernel = reshape(wV, k, k, nc, size(wV, 2));
    [M, W] = kernelEig(kernel, [size(X, 1), size(X, 2)]);

    maps = cropEigvec(M, W, c);

end

function [MSE, minMSE, optParam] = estimateMSE(C, r, imSize, weight, stdev, lst_k, lst_c, lst_wnsvn)
    nc = size(C, 3);
    im = ifft2c(zpad(C, imSize(1), imSize(2), nc)); 

    I = zeros([imSize(1), imSize(2), nc, nc]);
    for idx=1:1:nc
        I(:, :, idx, idx) = ones(imSize);
    end
    
    if (weight == false)
        assert(nargin >= 8, 'If not weighting, please pass in lst_wnsnv')
        MSE = zeros([length(lst_k), length(lst_c), length(lst_wnsvn)]);
    else
        MSE = zeros([length(lst_k), length(lst_c), 1]);
    end

    kdx = 0;
    for k = lst_k
        fprintf('| Kernel size: %d\n', k);
        kdx = kdx + 1;

        A = calreg2calmat(C, k); 

        if (weight == true)
            MSE(kdx, :)	= softMSE(A, im, k, lst_c, stdev, I);
        else
            MSE(kdx, :, :) = hardMSE(A, im, k, lst_wnsvn, lst_c, stdev, I);
	      end
    end

    if (weight == true)
        [minMSE, idx] = min(MSE(:));
        [kdx, cdx] = ind2sub(size(MSE), idx);
        optParam = [lst_k(kdx), lst_c(cdx)];
    else
        [minMSE, idx] = min(MSE(:));
        [kdx, cdx, wdx] = ind2sub(size(MSE), idx);
        optParam = [lst_k(kdx), lst_c(cdx), lst_wnsvn(wdx)];
    end
end

function MSE = pointMSE(M, W, c, im, stdev, I)

    nc = size(im, 3);
    maps = cropEigvec(M, W, c);

    pimg = squeeze(sum(conj(maps) .* repmat(im, [1, 1, 1, nc]), 3)); %proj img
    proj = zeros(size(im)); %coil projections

    for idx=nc:-1:1
        proj = proj + repmat(pimg(:, :, idx), [1, 1, nc]) .* maps(:, :, :, idx);
    end

    null = im - proj; % Null projection.

    T = zeros(size(maps)); % To calculate the projection maps (M * M')
    for pdx=1:1:nc
        for qdx=1:1:nc
            T(:, :, pdx, qdx) = sum(squeeze(maps(:, :, pdx, :)) .* squeeze(conj(maps(:, :, qdx, :))), 3);
        end
    end
	
    NP = T - I; % Calculating the null projector matrix.

    div = 0; %Calculating the divergence for SURE.
    for pdx=1:1:nc
        tmp = NP(:, :, pdx, pdx);
        div = div + sum(tmp(:));
    end

    % From SURE.
    MSE = prod(size(im)) * stdev^2 + norm(null(:), 2)^2 + stdev^2 * div;

end

function MSE = softMSE(A, im, k, lst_c, stdev, I)

    nc = size(im, 3);
    MSE = zeros([length(lst_c), 1]);

    [U, S, V] = svd(A); s = diag(S);

    weights = svWeightsSURE(s, stdev, size(A), true);

    if (length(weights) < size(V, 2)) 
        weights = [weights; zeros([size(V, 2) - length(weights), 1])];
    end 

    Vt = V * diag(weights); 

    kernel = reshape(Vt, k, k, nc, size(Vt, 2));
    [M, W] = kernelEig(kernel, [size(im, 1), size(im, 2)]);

    cdx = 0;
    for c = lst_c
        cdx = cdx + 1;
        MSE(cdx) = pointMSE(M, W, c, im, stdev, I);
        fprintf('|   |   | Crop threshold: %f ------------------------ | %0.6f\n', c, MSE(cdx));
    end
end

function MSE = hardMSE(A, im, k, lst_wnsvn, lst_c, stdev, I)

    nc = size(im, 3);
    MSE = zeros([length(lst_c), length(lst_wnsvn)]);

    [U, S, V] = svd(A); s = diag(S);

    wdx = 0;
    for w = lst_wnsvn
        wdx = wdx + 1;

        fprintf('|   |wnsvn: %f\n', w);
        n = round(k * k * w);
        weights = [ones([n 1]); zeros([size(V, 2) - n, 1])];

        Vt = V * diag(weights); 

        kernel = reshape(Vt, k, k, nc, size(Vt, 2));
        [M, W] = kernelEig(kernel, [size(im, 1), size(im, 2)]);

        cdx = 0;
        for c = lst_c
            cdx = cdx + 1;
            MSE(cdx, wdx) = pointMSE(M, W, c, im, stdev, I);
            fprintf('|   |   | Crop threshold: %f ------------------------ | %f\n', c, MSE(cdx, wdx));
        end
    end
end
