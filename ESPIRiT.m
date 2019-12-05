function maps = ESPIRiT(X, r, k, c, w)
  % maps = ESPIRiT(X, r, k, c, w)
  %
  % Function that calculates standard ESPIRiT maps. 
  %
  % INPUTS:
  %	  X     - k-space data [kx, ky, nc]
  %	  r     - size of calibration region
  %	  k     - kernel size
  %	  c     - crop threshold
  %	  w     - window normalized number of singular values (determines subspace size)
  % 
  % OUTPUTS
  %   maps  - espirit sensitivity maps [kx, ky, nc, nc]
  %
  % Copyright 2016. The Regents of the University of California.
  %
  % Authors:
  % 2016 Siddharth Iyer <sid8795@berkeley.edu>
  nc = size(X, 3);
  C = extractCalreg(X, r); 
  imSize = [size(X, 1), size(X, 2)];
  A = calreg2calmat(C, k); [U, S, V] = svd(A); s = diag(S);
  n = max(floor(wnsvn * k * k), 1);
  vecWeight = ones([n, 1]);
  vecWeight = [vecWeight; zeros([size(V, 2) - n, 1])];
  wV = V * diag(vecWeight);
  kernel = reshape(wV, k, k, nc, size(wV, 2));
  [M, W] = kernelEig(kernel, [size(X, 1), size(X, 2)]);
  maps = cropEigvec(M, W, c);
end
