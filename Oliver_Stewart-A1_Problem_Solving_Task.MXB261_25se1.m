function xys = simulate(N, P, sew)
% simulate : Function to simulate a biased random walk
%
% N   : arg, the number of particles to simulate
% P   : arg, the starting postion, can be an enum of 1 or 'rand'
% sew : arg, array of size 3, indicating the probabilities of any 
%      particle navigating South, East, and West (no North)
%
% particles : return value, an array of size N containing x, y 
%      coordinates for the final position of the simulated particles

    % change in x and y with each direction 
    %          S;       E;       W
    dxdys = [[0, -1]; [1, 0]; [-1, 0]];

    % grid size is maxxy x maxxy, e.g.: maxxy = 99 would result in a 99x99 grid
    maxxy = 99;
    
    % calculate the cumulative distribution function
    cdf = cumsum(sew);
    
    % initialise
    xys = zeros(N, 2); % rows are particles, columns are x-posn, y-posn
    stopped = false;
    
    for i = 1:N
        if P == 1
            xys(i, :) = [ceil(maxxy / 2), maxxy];
        end
        if P == 'rand'
            xys(i, :) = [randi([0, maxxy]), maxxy];
        end
        while stopped == false
            % use a uniform random number to sample PMF
            u = rand;
            dirn_idx = find(cdf >= u, 1);
        
            % find the to-be x and y position
            dxdy = dxdys(dirn_idx, :);
            tobexy = xys(i, :) + dxdy;
            
            % colision logic
            if tobexy(1) == -1 % colide with West border, wrap to East
                tobexy(1) = maxxy;
            end

            if tobexy(1) == maxxy + 1 % colide with East border, wrap to West
                tobexy(1) = 0;
            end
            
            if tobexy(2) == -1 % colide with South border
                stopped = true;
                continue;
            end

            for ii = 1:N
                if i == ii | ~isequal(tobexy, xys(ii, :))
                    continue;
                end
                if dxdy == dxdys(1, :) % colide vertically
                    tobexy = xys(i, :); % reset move
                    stopped = true; % stop outer while loop
                    break;
                end
                % default will trigger on coliding and not coliding vertically, coliding horizontally
                tobexy = xys(i, :); % reset move without stopping, regen. move
                break;
            end
        
            % update stored of results
            xys(i, :) = tobexy;
        end
        stopped = false;
    end

    for i = 1:N
        xys(i, :) = xys(i, :) - ceil(maxxy / 2);
    end
end

function graph(N, P, figurePrefix)
% graph : Function to graph simulated random walks performed by the simulate function, 
%     will plot the result of simulate N, P, for each of sews:
%        1/3, 1/3, 1/3;
%        2/3, 1/6, 1/6;
%        3/5, 1/10, 3/10;
%        3/5, 3/10, 1/10;
%
% N            : arg, the number of particles to simulate
% P            : arg, the starting postion, can be an enum of 1 or 'rand'
% figurePrefix : arg, the prefix for the figure title

    sews = [
        1/3, 1/3, 1/3;
        2/3, 1/6, 1/6;
        3/5, 1/10, 3/10;
        3/5, 3/10, 1/10;
    ];
    
    figure('Units','pixels','Position',[100 100 1200 800]);
    for i = 1:size(sews,1)
        sew = sews(i,:);
        final_positions = simulate(N, P, sew);
        x_coords = final_positions(:,1);
        
        subplot(2,2,i);
        h = histogram(x_coords, 'BinEdges', -50:50);
        
        ax = gca;
        ax.XAxisLocation = 'origin';
        
        title(sprintf('SEW %d', i));
        ylabel('Frequency');
        
        xlim([-50 49]);
        ylim([0 max(h.Values)+1]);
    end
    sgtitle(sprintf('%s: Random Walk Distribution for N: %d, P: %s', figurePrefix, N, mat2str(P)));
end

% Figures 1.a-b: Random Walk Distribution for N, P
graph(100, 1, 'Figure 1.a');
graph(200, 1, 'Figure 1.b');
graph(100, 'rand', 'Figure 1.c');
graph(200, 'rand', 'Figure 1.d');

% Define the true pmf
lambda = 4;
k = 0:15;
Prob_k = (lambda .^ k * exp(-lambda)) ./ factorial(k);
true_pmf = Prob_k / sum(Prob_k);

% Sampling Implementation
true_cdf = cumsum(true_pmf);

function sampled = sample(N, true_cdf)
% sampled : 
%
% N        : arg, the number of samples to sample from the cdf
% true_cdf : arg, the cdf to sample from
    sampled = zeros(N, 1);
    for i = 1:N
        % use a uniform random number to sample PMF
        u = rand;
        sampled(i, 1) = find(true_cdf >= u, 1, 'first') - 1; % -1 cause MATLAB indexes from 1, not 0 like every other language
    end
end

% KL Divergence Calculation
sample_sizes = [10, 25, 60, 100, 175, 250];
kld_mean = zeros(6, 1);
kld_se = zeros(6, 1);
example_pmfs = zeros(6, 16);

for i = 1:6
    kl_divergences = zeros(100, 1);
    for ii = 1:100
        curr_sample = sample(sample_sizes(i), true_cdf);
        counts = histcounts(curr_sample, 0:16);
        empirical_pmf = counts ./ sum(counts);
        kl_divergences(ii, 1) = sum(true_pmf .* log(true_pmf ./  max(empirical_pmf, 1e-10)));% avoid numerical issues with log(0)
        if ii == 1
            example_pmfs(i, :) = empirical_pmf; % sample pmf for plotting
        end
    end

    kld_mean(i, 1) = mean(kl_divergences);
    kld_se(i, 1) = std(kl_divergences) / sqrt(100);
end

% Figure 2.a: KL Divergence by Sample Size
figure;
errorbar(sample_sizes, kld_mean(:, 1), kld_se(:, 1));
sgtitle('Figure 2.a: KL Divergence by Sample Size')
xlabel('Sample Size');
ylabel('KL Divergence Mean');

% Figure 2.b: True PMF vs. Empirical PMF per. sample size
figure('Units','pixels','Position',[100 100 1200 800]);
for i = 1:6  
    subplot(2,3,i);
    bar(0:15, [true_pmf(:), example_pmfs(i,:).'], 'grouped');
    
    title(sprintf('No. Samples: %d', sample_sizes(i)));
    xlabel('x')
    ylabel('Probability')

    if i == 3
        legend('True PMF','Empirical PMF');
    end
end
sgtitle('Figure 2.b: True PMF vs. Empirical PMF per. sample size');