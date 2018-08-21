function proj = projImg(im, maps)
    nc = size(  im, 3);
    m  = size(maps, 4);
    pimg = squeeze(sum(conj(maps) .* repmat(im, [1, 1, 1, m]), 3)); %proj img
    proj = zeros(size(im)); %coil projections
    for kdx=1:1:m
        proj = proj + repmat(pimg(:, :, kdx), [1, 1, nc]) .* maps(:, :, :, kdx);
    end
end
