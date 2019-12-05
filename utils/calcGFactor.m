function gmap = calcGFactor(m, rx, ry)
% Reference: /https://web.stanford.edu/class/ee369c/restricted/Solutions/assignment_4_solns.pdf
  gmap = zeros(size(m, 1), size(m, 2));
  assert(mod(size(m, 1), rx) == 0);
  assert(mod(size(m, 2), ry) == 0);
  shift_x = size(m, 1)/rx;
  shift_y = size(m, 2)/ry;
  f = @(x) x(:);
  mat = Inf;
  for x=1:1:size(m, 1)
    for y=1:1:size(m, 2)
      for xx=0:(rx-1)
        for yy=0:(ry-1)
          if (xx == 0 && yy == 0)
            mat = f(m(x, y, :));
          else
            xdx = mod((x - 1) + xx * shift_x, size(m, 1)) + 1;
            ydx = mod((y - 1) + yy * shift_y, size(m, 2)) + 1;
            mat = [mat, f(m(xdx, ydx, :))];
          end
        end
      end
      scs  = mat' * mat;
      scsi = pinv(scs);
      cnd  = sqrt(scsi(1, 1) * scs(1, 1));
      gmap(x, y) = abs(cnd);
    end
  end
end
