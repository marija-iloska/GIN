function [fs] = in_f(A, I, I0, K, A_init, C_est, mu_x, sig_x, gamma, lambda_init)

for prior = 1:length(gamma)
    
    tic
    [A_s4] = gibbs_indegree(I, I0, K, A_init, lambda_init, C_est, mu_x, sig_x, gamma(prior));
    toc 
    
    % Estimate as the mode
    A_est_4 = mode(A_s4, 3);

    % Performance
    [~, ~, fs(prior)] = adj_eval(A, A_est_4);    

end


end