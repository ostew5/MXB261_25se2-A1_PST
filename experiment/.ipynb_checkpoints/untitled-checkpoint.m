% Define the true pmf
lambda = 4;

k = 0:15;

Prob_k = (lambda .^ k * exp(-lambda)) ./ factorial(k);

true_pmf = Prob_k / sum(Prob_k);

% Sampling Implementation
true_cdf = cumsum(true_pmf);

function sampled = sample(N, true_cdf)
    sampled = zeros(N, 1);
    for i = 1:N
        % use a uniform random number to sample PMF
        u = rand;
        sampled(i, 1) = find(true_cdf >= u, 1, 'first'); % -1 cause MATLAB indexes from 1, not 0 like every other language
    end
end

sample_sizes = [10, 25, 60, 100, 175, 250];
kld_mean = zeros(6, 1);
kld_se = zeros(6, 1);

% KL Divergence Calculation
for i = 1:6
    kl_divergences = zeros(100, 1);
    for ii = 1:100
        curr_sample = sample(sample_sizes(i, 1), true_cdf);
        counts = histcounts(curr_sample, 0:16);
        empirical_pmf = max(counts ./ sum(counts), 1e-10); % avoid numerical issues with log(0)
        kl_divergences(ii, 1) = sum(true_pmf .* log(true_pmf ./ empirical_pmf));
    end

    kld_mean(i, 1) = mean(kl_divergences);
    kld_se(i, 1) = std(kl_divergences) / sqrt(100);
end