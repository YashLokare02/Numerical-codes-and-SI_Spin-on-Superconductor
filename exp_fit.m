% Exponential fits -- to estimate the superconducting coherence length along different line cuts

function [coherence_len, para_vals] = exp_fit(data1, data2, guess_1, guess_2)

% Temporary assignment
data_real = data1;

% Normalization of input datasets
data1 = (data1 - min(data1))./(max(data1 - min(data1))); % scan range values (1D)
data2 = (data2 - min(data2))./(max(data2 - min(data2))); % autocorrelation values

% Initial guesses
beta0 = [guess_1; guess_2]; % guess_1 -- for coefficient; guess_2 -- for length scale

% Statistical measures
mu_distance = mean(data_real);
std_distance = std(data_real);

% Fitting to an exponential model
f = @(b, data1) b(1).*exp(-data1./b(2)); % Exponential model

% Parametric fitting using fminsearch()
para_vals = fminsearch(@(b) norm(data2 - f(b, data1)), beta0); % Estimating fitting parameters

coherence_len = para_vals(2)*std_distance + mu_distance; % estimated superconducting coherence length

figure(1);
plot(data1, data2, 'pg');
hold on
plot(data1, f(para_vals, data1), 'r');
grid
xlabel('Distance (arb. units)', 'FontSize', 18);
ylabel('ACF (arb. units)', 'FontSize', 18);
legend('Autocorrelation data', 'Exponential fit', 'FontSize', 16, 'Orientation', 'vertical', 'Location', 'best');

end