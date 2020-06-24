function A = adj_sig_level(C, alpha, Method)
% Convert a coefficient matrix into an adjacency matrix using a heuristic
% approach inspired by choosing the most significant contributors to the
% norm of the matrix
dim = length(C(:, 1));
m = C(:);
if(Method == 0)
    %% Submethod 1: The Whole Matrix
    m = C(:);
    [sort_m, idx] = sort(abs(m)/sum(abs(m)), 'descend');
    mask = cumsum(sort_m)<(1-alpha);
    
    % Obtain topology
    A = m;
    A(idx(mask)) = 1;
    A(idx(~mask)) = 0;
    A = reshape(A, dim, dim);
else
    %% Submethods 2 & 3: By Row or By Column
    A = zeros(dim, dim);
    for i = 1:dim
        if(Method == 1)
            m = C(i, :);
        elseif(Method == 2)
            m = C(:, i);
        end
        [sort_m, idx] = sort(abs(m)/sum(abs(m)), 'descend');
        mask = cumsum(sort_m)<(1-alpha);
        % Obtain topology
        a = m;
        a(idx(mask)) = 1;
        a(idx(~mask)) = 0;
        if(Method == 1)
            A(i, :) = a;
        elseif(Method == 2)
            A(:, i) = a;
        end
    end
end
end

