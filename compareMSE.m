function [trueMSE, trueOptParam, sureFullMSE, sureFullOptParam, sureCalMSE, sureCalOptParam, gmax] = compareMSE(fileID, X, Y, r, weight, stdev, lst_k, lst_c, lst_wnsvn, gfactor)

  if (nargin < 10)
    gfactor = false;
  end

  imSize = [size(Y, 1), size(Y, 2)];
  nc     = size(Y, 3);
  C      = extractCalreg(Y, r);
  x      = ifft2c(X);
  yf     = ifft2c(Y);
  yc     = ifft2c(zpad(C, imSize(1), imSize(2), nc));
  I      = zeros([imSize(1), imSize(2), nc, nc]);
  for idx=1:1:nc
    I(:, :, idx, idx) = ones(imSize);
  end
    
  if (weight == false)
    assert(nargin >= 8, 'If not weighting, please pass in lst_wnsnv.')
    trueMSE     = zeros([length(lst_k), length(lst_c), length(lst_wnsvn)]);
    sureFullMSE = zeros([length(lst_k), length(lst_c), length(lst_wnsvn)]);
    sureCalMSE  = zeros([length(lst_k), length(lst_c), length(lst_wnsvn)]);
    if (gfactor)
      gmax      = zeros([length(lst_k), length(lst_c), length(lst_wnsvn)]);
    else
      gmax      = -1;
    end
  else
    trueMSE     = zeros([length(lst_k), length(lst_c)]);
    sureFullMSE = zeros([length(lst_k), length(lst_c)]);
    sureCalMSE  = zeros([length(lst_k), length(lst_c)]);
    if (gfactor)
      gmax      = zeros([length(lst_k), length(lst_c)]);
    else
      gmax      = -1;
    end
  end

  kdx = 0;
  for k = lst_k
    fprintf(fileID, '| k = %d\n', k);
    kdx = kdx + 1;
    A = calreg2calmat(C, k); 
    if (weight == true)
      if (gfactor)
        [trueMSE(kdx, :), sureFullMSE(kdx, :), sureCalMSE(kdx, :), gmax(kdx, :)] = ...
          softMSE(A, r, k, lst_c, stdev, I, x, yf, yc, gfactor, fileID);
      else
        [trueMSE(kdx, :), sureFullMSE(kdx, :), sureCalMSE(kdx, :), ~] = ...
          softMSE(A, r, k, lst_c, stdev, I, x, yf, yc, gfactor, fileID);
      end
    else
      if (gfactor)
        [trueMSE(kdx, :, :), sureFullMSE(kdx, :, :), sureCalMSE(kdx, :, :), gmax(kdx, :)] = ...
          hardMSE(A, r, k, lst_c, lst_wnsvn, stdev, I, x, yf, yc, gfactor, fileID);
      else
        [trueMSE(kdx, :, :), sureFullMSE(kdx, :, :), sureCalMSE(kdx, :, :), ~] = ...
          hardMSE(A, r, k, lst_c, lst_wnsvn, stdev, I, x, yf, yc, gfactor, fileID);
      end
    end
  end

  if (weight == true)
	  [trueMinMSE, idx]     = min(trueMSE(:));
    [kdx, cdx]            = ind2sub(size(trueMSE), idx);
    trueOptParam          = [lst_k(kdx), lst_c(cdx)];

    [sureFullMinMSE, idx] = min(sureFullMSE(:));
    [kdx, cdx]            = ind2sub(size(sureFullMSE), idx);
    sureFullOptParam      = [lst_k(kdx), lst_c(cdx)];
     
    [sureCalMinMSE, idx]  = min(sureCalMSE(:));
    [kdx, cdx]            = ind2sub(size(sureCalMSE), idx);
    sureCalOptParam       = [lst_k(kdx), lst_c(cdx)];

    fprintf(fileID, '\nTrue Optimal Parameters:\n   k: %d\n   c: %f\n',        trueOptParam(1),     trueOptParam(2));
    fprintf(fileID, '\nSURE Full Optimal Parameters:\n   k: %d\n   c: %f\n',   sureFullOptParam(1), sureFullOptParam(2));
    fprintf(fileID, '\nSURE CalReg Optimal Parameters:\n   k: %d\n   c: %f\n', sureCalOptParam(1),  sureCalOptParam(2));
  else
    [trueMinMSE, idx]     = min(trueMSE(:));
    [kdx, cdx, wdx]       = ind2sub(size(trueMSE), idx);
    trueOptParam          = [lst_k(kdx), lst_c(cdx), lst_wnsvn(wdx)];

    [sureFullMinMSE, idx] = min(sureFullMSE(:));
    [kdx, cdx, wdx]       = ind2sub(size(sureFullMSE), idx);
    sureFullOptParam      = [lst_k(kdx), lst_c(cdx), lst_wnsvn(wdx)];

    [sureCalMinMSE, idx]  = min(sureCalMSE(:));
    [kdx, cdx, wdx]       = ind2sub(size(sureCalMSE), idx);
    sureCalOptParam       = [lst_k(kdx), lst_c(cdx), lst_wnsvn(wdx)];

    fprintf(fileID, '\nTrue Optimal Parameters:\n   k: %d\n   c: %f\n   wnsvn: %f\n', ...
      trueOptParam(1), trueOptParam(2), trueOptParam(3));
    fprintf(fileID, '\nSURE Full Optimal Parameters:\n   k: %d\n   c: %f\n   wnsvn: %f\n', ...
      sureFullOptParam(1), sureFullOptParam(2), sureFullOptParam(3));
    fprintf(fileID, '\nSURE CalReg Optimal Parameters:\n   k: %d\n   c: %f\n   wnsvn: %f\n', ...
      sureCalOptParam(1), sureCalOptParam(2), sureCalOptParam(3));
  end
end

function [trueMSE, sureFullMSE, sureCalMSE, gmax] = pointMSE(M, W, r, c, stdev, I, x, yf, yc, gfactor);
  nc   = size(x, 3);
  maps = cropEigvec(M, W, c);
  T    = zeros(size(maps)); % To calculate the projection operator (maps * maps')
  for pdx=1:1:nc
    for qdx=1:1:nc
      T(:, :, pdx, qdx) = sum(squeeze(maps(:, :, pdx, :)) .* squeeze(conj(maps(:, :, qdx, :))), 3);
    end
  end

  NP  = T - I; % Calculating the null projector matrix.
  div = 0;     % Calculating the divergence for SURE.
  for pdx=1:1:nc
    tmp = NP(:, :, pdx, pdx);
    div = div + sum(tmp(:));
  end

  % Calculate projection.
  pimg = squeeze(sum(conj(maps) .* repmat(yf, [1, 1, 1, nc]), 3)); % Proj img.
  proj = zeros(size(yf));                                          % Coil projections.
  for idx=nc:-1:1
    proj = proj + repmat(pimg(:, :, idx), [1, 1, nc]) .* maps(:, :, :, idx);
  end

  % true MSE
  null    = x - proj;
  trueMSE = norm(null(:), 2)^2;

  % SURE Full kspace
  null = yf - proj; % Null projection.
  sureFullMSE = prod(size(yf)) * stdev^2 + norm(null(:), 2)^2 + 2 * stdev^2 * div;

  % for SURE CalReg
  T = zeros(size(maps)); % To calculate the projection operator (lowres * maps * maps' * lowres')
  for pdx=1:1:nc
    for qdx=1:1:nc
      T(:, :, pdx, qdx) = sum( ... 
             ifft2c(zpad(extractCalreg(fft2c(squeeze(maps(:, :, pdx, :))), r), size(x))) .* ... 
        conj(ifft2c(zpad(extractCalreg(fft2c(squeeze(maps(:, :, qdx, :))), r), size(x)))), 3);
    end
  end 

  div = 0; %Calculating the divergence for Cal SURE.
  for pdx=1:1:nc
    tmp = T(:, :, pdx, pdx);
    div = div + sum(tmp(:));
  end 

  % Take the last projection and extract the lower frequencies.
  proj = ifft2c(zpad(extractCalreg(fft2c(proj), r), size(yc)));
  null = yc - proj; % Null projection (of low freqs).
  sureCalMSE = (r * r * nc) * stdev^2 + norm(null(:), 2)^2 + 2 * stdev^2 * div;

  if (gfactor)
    gmax = calcGFactor(squeeze(maps(:,:,:,end)), 2, 2);
    gmax = gmax(abs(gmax(:)) > 0);
    gmax = max(gmax(:));
    if (numel(gmax) < 1)
      gmax = 0;
    end
  else
    gmax = -1;
  end
end

function [trueMSE, sureFullMSE, sureCalMSE, gmax] = softMSE(A, r, k, lst_c, stdev, I, x, yf, yc, gfactor, fileID)
  nc          = size(x, 3);
  trueMSE     = zeros([length(lst_c), 1]);
  sureFullMSE = zeros([length(lst_c), 1]);
  sureCalMSE  = zeros([length(lst_c), 1]);
  if (gfactor)
    gmax      = zeros([length(lst_c), 1]);
  else
    gmax      = -1;
  end

  [U, S, V] = svd(A); s = diag(S);
  weights   = svWeightsSURE(s, stdev, size(A));

  if (length(weights) < size(V, 2)) 
    weights = [weights; zeros([size(V, 2) - length(weights), 1])];
  end 

  Vt     = V * diag(weights); 
  kernel = reshape(Vt, k, k, nc, size(Vt, 2));
  [M, W] = kernelEig(kernel, [size(x, 1), size(x, 2)]);

  cdx = 0;
  for c = lst_c
    fprintf(fileID, '| | c = %0.3f\n', c);
    cdx = cdx + 1;
    if (gfactor)
      [trueMSE(cdx), sureFullMSE(cdx), sureCalMSE(cdx), gmax(cdx)] = pointMSE(M, W, r, c, stdev, I, x, yf, yc, gfactor);
    else
      [trueMSE(cdx), sureFullMSE(cdx), sureCalMSE(cdx), ~] = pointMSE(M, W, r, c, stdev, I, x, yf, yc, gfactor);
    end
  end
end

function [trueMSE, sureFullMSE, sureCalMSE, gmax] = hardMSE(A, r, k, lst_c, lst_wnsvn, stdev, I, x, yf, yc, gfactor, fileID)
  nc          = size(x, 3);
  trueMSE     = zeros([length(lst_c), length(lst_wnsvn)]);
  sureFullMSE = zeros([length(lst_c), length(lst_wnsvn)]);
  sureCalMSE  = zeros([length(lst_c), length(lst_wnsvn)]);
  if (gfactor)
    gmax      = zeros([length(lst_c), length(lst_wnsvn)]);
  else
    gmax      = -1;
  end

  [U, S, V] = svd(A); s = diag(S);

  wdx = 0;
  weights = zeros([size(V, 2), 1]);
  for w = lst_wnsvn
    fprintf(fileID, '| | w = %0.3f\n', w);
    wdx = wdx + 1;
    weights(1:round(k * k * w)) = 1;
    Vt = V * diag(weights); 

    kernel = reshape(Vt, k, k, nc, size(Vt, 2));
    [M, W] = kernelEig(kernel, [size(x, 1), size(x, 2)]);

    cdx = 0;
    for c = lst_c
      fprintf(fileID, '| | | c = %0.3f\n', c);
      cdx = cdx + 1;
      if (gfactor)
        [trueMSE(cdx, wdx), sureFullMSE(cdx, wdx), sureCalMSE(cdx, wdx), gmax(cdx, wdx)] = pointMSE(M, W, r, c, stdev, I, x, yf, yc, gfactor);
      else
        [trueMSE(cdx, wdx), sureFullMSE(cdx, wdx), sureCalMSE(cdx, wdx), ~] = pointMSE(M, W, r, c, stdev, I, x, yf, yc, gfactor);
      end
    end
  end
end
