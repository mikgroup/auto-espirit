function weights = svWeightsSURE(s, stdev, calmatDim)

    n       = length(s);
    estMSE  = zeros([n 1]);
    Y       = prod(calmatDim) * stdev^2;

    for idx=1:1:n
        lambda      = s(idx);
        G           = sum(min(lambda, s).^2);
        div         = calcDiv(lambda, s, calmatDim);
        estMSE(idx) = -Y + G + stdev^2 * div;
    end 

    lambda  = s(min(find(min(estMSE) == estMSE)));
    weights = max(0, s - lambda)./s;

end

function x = calcDiv(lambda, s, dimA)
    % Divergence term for low-rank matrix estimate using SURE. 
    % Code modified from: http://statweb.stanford.edu/~candes/SURE/index.html 
    svThreshold = 1E-8; % threshold to determine whether two singular values are the same
    M = dimA(1);
    N = dimA(2);
    % safeguard.
    % check multiplicities of singular values in a robust manner
    z  = s(2:end);
    s  = [s(1) 1];
    Is = 1;
    while(~isempty(z)),
        idx = find(abs(z - s(Is, 1)) < svThreshold);
        if(isempty(idx))
            s    = [s; [z(1) 1]];
            z(1) = [];
            Is   = Is + 1;
        end
        z(idx)   = [];
        s(Is, 2) = s(Is, 2) + numel(idx);
    end
    clear z
    % warns the user about using SURE with a non-simple, not-full
    % rank matrix
    if(any(s(:, 1) < svThreshold)),
        fprintf('   +   [SURE_SVT] Warning: argument might be rank-deficient.\n')
    end
    if(any(s(:, 2) > 1)),
        fprintf('   +   [SURE_SVT] Warning: argument might have repeated singular values.\n')
    end
    % find singular values above the threshold
    idx_p = (s(:, 1) > lambda);
    x     = 0;
    if(any(idx_p)),
        x = x + sum( s(idx_p, 2).^2 );
        x = x + sum( (2*abs(M-N) + 1 + s(idx_p, 2).*(s(idx_p, 2) - 1)).*(max(0, s(idx_p, 1) - lambda)./s(idx_p, 1)) );
    end
    D = zeros(size(s, 1));
    for Ik = 1:size(s, 1),
        D(:, Ik) = s(Ik, 2)*s(:, 2).*s(:, 1).*max(0, s(:, 1) - lambda)./(s(:, 1).^2 - s(Ik, 1).^2);
    end
    D(isnan(D) | isinf(D) | abs(D) > 1E6) = 0;
    x = x + 4 * sum(D(:));
end
