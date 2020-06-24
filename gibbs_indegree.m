function [A_s, lambda_s] = gibbs_indegree(iter_num, burn_in, thin_rate, A_init, lambda_init, C, mu, sig, gamma)
% Obtain dimensions
dim_y = length(C(:, 1));
% Initialize arrays to store samples
A_s = zeros(dim_y, dim_y, iter_num+1);
lambda_s = zeros(iter_num+1, 1);
% Set the initial states of the Markov chains
A_s(:, :, 1) = A_init;
lambda_s(1) = lambda_init;
A_old = A_s(:, :, 1);
% For loop for the Gibbs
for i = 2:iter_num+1
    for j = 1:dim_y
        for k = 1:dim_y
            
            % Conside that A(j, k)=1
            A_old(j,k) = 1;
            C_temp = C.*A_old; 
            log_pa1 = -lambda_s(i-1)*sum(A_old(j,:)) + logmvnpdf(C_temp(:)', mu', sig);
           
            % Conside that A(j, k)=0
            A_old(j,k) = 0;
            C_temp = C.*A_old; 
            log_pa0 = -lambda_s(i-1)*sum(A_old(j,:)) + logmvnpdf(C_temp(:)', mu', sig);
           
            % Sample the topology
            pa0 = exp(log_pa0-max([log_pa0, log_pa1]));
            pa1 = exp(log_pa1-max([log_pa0, log_pa1]));
            prob_1 = pa1/(pa1+pa0);
            A_old(j,k) = rand<prob_1;
           
            %lambda_s(j,k,i) = exprnd(gamma + sum(A_old(j,:) + A_old(:,k)'));
        end
    end
    lambda_s(i) = gamrnd(dim_y+1, gamma + dim_y*sum(A_old(j,:)));
    A_s(:, :, i) = A_old; 
    % Sample a new value of lambda
    %lambda_s(i) = exprnd(gamma+sum(sum(A_old)));

end
% Apply burn-in 
A_s = A_s(:, :, burn_in+1:thin_rate:iter_num+1);
lambda_s = lambda_s(burn_in+1:thin_rate:iter_num+1);

end

