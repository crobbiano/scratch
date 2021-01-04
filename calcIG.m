function [ig, entropy_alphas] = calcIG(old_alphas, class_label)
%calcIGCone calculates the information gain given the current state (old_alphas)
% and a new measurement
% This requires that the alphas related to class 1 is in index 1 of
% old_alphas, class 2 in index 2, etc.  
% This also requires that the class label is the numeric version of the
% class, so if cars are class 2, then class label would be 2


% For a single measurement, we update the pseudocount associated with the
% new class label
entropy_alphas = old_alphas;
entropy_alphas(class_label) = entropy_alphas(class_label) + 1;

alpha0 = sum(entropy_alphas);
entropy_lambdas = entropy_alphas./alpha0;

h_lambda = logBeta(entropy_alphas) + (alpha0 - length(entropy_alphas))*psi(alpha0)...
    - sum(arrayfun(@(x) (x-1)*psi(x), entropy_alphas));

% calculate the conditional entropy term h(\lambda | z)
% \sum_{z=1}^K \frac{\alpha_z}{\alpha_0} \left( \log{B(\ba')} + (\alpha_0' - K)\psi(\alpha_0') - \sum_{k=1}^K (\alpha_k' - 1)\psi(\alpha_k')\right) \right].
temp_alphas = entropy_alphas + 1;
delta_alphas = temp_alphas - entropy_alphas;

h_lambda_z = 0;
for z=1:length(delta_alphas)
    cond_entropy_alphas = entropy_alphas;
    cond_entropy_alphas(z) = cond_entropy_alphas(z) + delta_alphas(z);
    alpha0 = sum(cond_entropy_alphas);
    h_lambda_z = h_lambda_z + entropy_lambdas(z) * (logBeta(cond_entropy_alphas) ...
        + (alpha0 - length(cond_entropy_alphas))*psi(alpha0) ...
        - sum(arrayfun(@(x) (x-1)*psi(x), cond_entropy_alphas)));
end
ig = h_lambda - h_lambda_z;
end

function [outs] = logBeta(alpha)
%LOGBETA calculates the log of multivariate beta function
%   This recursively, and quickly, calculates the log beta
%   function.  If you tried to calculate from just the
%   definition, then everything breaks due to HUGE AF factors
%   from the factorial function for large values of alpha.

K = length(alpha);
alpha0 = sum(alpha);
logalpha0 = 0;
for i=1:alpha0-1
    logalpha0 = logalpha0 + log(i);
end

logalpha = zeros(1,K);
for k=1:K
    for j=1:alpha(k)-1
        logalpha(k) = logalpha(k) + log(j);
    end
end

outs = sum(logalpha) - logalpha0;
end
