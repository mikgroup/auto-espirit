function stdev = noistEst(C, s, k)
    E     = normsvdist(size(C), k);
    stdev = normSVNoiseEst(s, E);
end

function E = normsvdist(calibSize, k, numIters)

    if nargin < 4
        numIters = 5;
    end

    NC = normrnd(0, 1/sqrt(2), calibSize) + i * normrnd(0, 1/sqrt(2), calibSize);
    N = calib2mat(NC, k);
    [~, E, ~] = svd(N);
    E = diag(E);

    for idx=1:1:numIters-1
        NC = normrnd(0, 1/sqrt(2), calibSize) + i * normrnd(0, 1/sqrt(2), calibSize);
        N = calib2mat(NC, k);
        [~, tmpE, ~] = svd(N);
        E = E + diag(tmpE);
    end

    E = E/numIters;

end

function stdev = normSVNoiseEst(S, E)

    n = length(S);
    range = n:-1:n-n/2;

    stdev = mean(S(range)./E(range))/1.1;

end
