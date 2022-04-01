% To estimate the MZM localization length for our topological superconducting system 

% Extracting bias voltage data
[volt] = Nanonis_reader('Line Spectroscopy006.3ds', 1);
volt_data = volt(:, 3); %% Bias voltage values [V]
volt_data = volt_data.*1e3; % in mV

% Parameters
n = 64; % no. of grid points
rows = 128;
cols = n*3; % to store all the STM tunneling spectral data

% Initialization of the data matrix
data = zeros(rows, cols);
data(:, 1:3) = Nanonis_reader('Line Spectroscopy006.3ds', 0); % extracting data for grid point 0

% Data acquisition for the rest of the grid points
for i = 1:n-1
    data(:, 3*i+1:3*i+3) = Nanonis_reader('Line Spectroscopy006.3ds', i);
end

% Acquiring the STM data exclusively
STM_data = zeros(rows, n);

STM_data(:, 1) = data(:, 2);

for i = 2:n
    STM_data(:, i) = data(:, 3*i-1);
end

% Normalization of the acquired STM data
for i = 1:n
    STM_data(:, i) = (STM_data(:, i) - min(STM_data(:, i)))./(max(STM_data(:, i)) - min(STM_data(:, i))); % Normalization
end

% Defining line scan range
x = linspace(0, 3.42*1e-9, n); % 3.42 nm line scan considered

% Extracting zero bias STM data (to visuualize the spatially-resolved profile of the ZBCPs)
ZBCP_data = zeros(1, n);

for i = 1:n
    ZBCP_data(i) = STM_data(rows/2, i);
end

% Locating the index positions of max(ZBCP_data)
k = max(ZBCP_data); % to locate all the index positions of max(ZBCP_data)
idx = find(ZBCP_data == k);

% Relevant values of ZBCP_data and x to be considered
ZBCP_actual = ZBCP_data(idx(end):end);
x_real = x; % temporary storage of actual scan range values
x_actual = x(idx(end):end); 

% Fitting to an exponential model
x_actual = (x_actual - min(x_actual))./(max(x_actual) - min(x_actual)); % normalizing the scan range values

beta0 = [1.1; 0.04; 0.02]; % initial gueses for the parametric fit

% Statistical measures
mu_x = mean(x_real);
std_x = std(x_real);

% Fitting to an exponential model
f = @(b, x_actual) b(1).*exp(-x_actual./b(2)) + b(3); % Exponential model

% Parametric fitting using fminsearch()
para_vals = fminsearch(@(b) norm(ZBCP_actual - f(b, x_actual)), beta0); % Estimating fitting parameters

% Exponential fit function
yfit = f(para_vals, x_actual);

% Determining MZM localization length
MZM_len = para_vals(2)*std_x + mu_x; % estimated MZM localization length

% Visualization
figure(1); % for the ZBCP profile (full profile)
plot(x_real, smooth(ZBCP_data), 'k', 'LineWidth', 2);
grid
xlabel('Distance [m]', 'FontSize', 19);
ylabel('$dI/dV$ (arb. units)', 'Interpreter', 'latex', 'FontSize', 19);

figure(2); % for the exponential fit
plot(x_actual, ZBCP_actual, 'ko');
hold on
plot(x_actual, yfit, 'r', 'LineWidth', 2);
grid
xlabel('Distance (arb. units)', 'FontSize', 20);
ylabel('$dI/dV$ (arb. units)', 'Interpreter', 'latex', 'FontSize', 20);
legend('Differential conductance data', 'Exponential fit', 'FontSize', 19, 'Orientation', 'vertical', 'Location', 'best');

figure(3); % exponential fit with error bars included 
plot(x_actual, yfit, 'ro-', 'LineWidth', 2);
hold on
err = abs(ZBCP_actual - yfit); % error computation
eb = errorbar(x_actual, yfit, err, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1.5);
grid
xlabel('Distance (arb. units)', 'FontSize', 20);
ylabel('$dI/dV$ (arb. units)', 'Interpreter', 'latex', 'FontSize', 20);
legend('Exponential fit', 'Error', 'FontSize', 19, 'Orientation', 'vertical', 'Location', 'best');
hold off

















