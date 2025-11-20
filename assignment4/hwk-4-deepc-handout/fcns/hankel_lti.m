function H = hankel_lti(data, L)
% data: (n x T) matrix of time-series (n variables, T samples)
% L   : window length
%
% H   : (n*L x (T-L+1)) Hankel data matrix
%
% Column j of H is [ data(:,j); data(:,j+1); ...; data(:,j+L-1) ]

    [n, T] = size(data);
    cols = T - L + 1;
    if cols <= 0
        error('Not enough data points for the chosen window length L.');
    end

    H = zeros(n*L, cols);
    for j = 1:cols
        block = data(:, j:j+L-1);   % n x L
        H(:,j) = block(:);          % stack columns
    end