function y = F_inv(x, ax)
  % y = F_inv(x, ax) 
  %
  % Adjoint centered unitary Fourier transform along axis in ax.

  y = fftshift(ifft(ifftshift(x, ax(1)), size(x, ax(1)), ax(1)), ax(1)) .* sqrt(size(x, ax(1)));
  for k = 2:numel(ax)
    y = fftshift(ifft(ifftshift(y, ax(k)), size(y, ax(k)), ax(k)), ax(k)) .* sqrt(size(y, ax(k)));
  end

end
