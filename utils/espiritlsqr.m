function res = espiritlsqr(x, mask, maps)
  
  sx = size(maps, 1);
  sy = size(maps, 2);
  nc = size(maps, 3);
  md = size(maps, 4);

  function res = A(x, transp_flag)
    if strcmp(transp_flag,'transp')
      x = reshape(x, sx, sy, nc, 1);
      res = bsxfun(@times, x, mask);
      res = F_inv(res, [1, 2]);
      res = sum(bsxfun(@times, conj(maps), res), 3);
    elseif strcmp(transp_flag, 'notransp')
      x = reshape(x, sx, sy, 1, md);
      res = sum(bsxfun(@times, maps, x), 4);
      res = F_fwd(res, [1, 2]);
      res = bsxfun(@times, res, mask);
    end
    res = res(:);
  end

  res = reshape(lsqr(@(x, t)A(x, t), x(:), 1e-3, 100), sx, sy, 1, md);
end

