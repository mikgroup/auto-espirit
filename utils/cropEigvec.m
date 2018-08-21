function maps = cropEigvec(M, W, c)

    nc = size(M, 3);
    maps = zeros([size(M, 1), size(M, 2), nc, nc]);

    for idx=1:1:nc
        maps(:, :, :, idx) = M(:, :, :, end-idx+1) .* repmat(W(:, :, end-idx+1) >= c, [1, 1, nc]);
    end

end
