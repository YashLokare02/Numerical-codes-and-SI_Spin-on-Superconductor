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

% Trial solution (to facilitate a better fit of the full data)
STM_trial = STM_data;
data_trial = STM_data(rows/2, 1:end);

% Defining line scan range
x = linspace(0, 3.42*1e-9, n); % 3.42 nm line scan considered
x_real = x; % temporary storage of actual scan range values

% Normalization (for the full dataset)
x_real = (x_real - min(x_real))./(max(x_real) - min(x_real)); % normalization of the full line scan range

% Additional calculations (Fermi wavevector and Fermi wavelength; also period of oscillations of the LDOS data)
esp = (x(end) - x(1))/length(x); 
FFT_vals = fft(STM_data); % Fourier transform of the conductance map
ZBCP_FFT = real(FFT_vals(rows/2, 1:end)); % FFT of zero bias data
[~, locs] = findpeaks(ZBCP_FFT); % locate peaks 
T = mean(diff(locs)*esp); % period of oscillations
k_Fermi = pi/T; % Fermi wavevector
lambda_Fermi = 2*pi/k_Fermi; % Fermi wavelength

% Extracting zero bias STM data (to visuualize the spatially-resolved profile of the ZBCPs)
ZBCP_data = zeros(1, n);

for i = 1:n
    ZBCP_data(i) = STM_data(rows/2, i);
end

% Fitting the theoretical MBS probability density to the ZBCP data
[yfit_complete, gof_2] = MZM_prob_density(x_real, data_trial); % for probability density fit (full dataset)

% Computng parametric fitting errors
y_vals = yfit_complete(x_real);
err = abs(y_vals' - data_trial);
err(1:3) = 0; % deleting irrelevant error values

% Visualization of the fit
figure(2); % for the ZBCP profile (full profile)
plot(x_real, smooth(ZBCP_data), 'k', 'LineWidth', 2);
grid
xlabel('Distance [m]', 'FontSize', 19);
ylabel('$dI/dV$ (arb. units)', 'Interpreter', 'latex', 'FontSize', 19);

% Color map for dI/dV vs. distance and bias voltage (Supplementary Information; not required as such)
% Symmetrization of the STM data
for i = 1:n
    STM_data((rows/2)+1:rows, i) = STM_data(rows/2:-1:1, i); % symmetrization
end

% Rescaling
x = x./1e-9; % in nm

% Visualization
figure(3);
surf(x, volt_data, STM_data)
brighten(0.6);
view([0 0 90])
hold on % holds the opened figure for further modification
shading interp % interpolates datapoints, ideal for 3D graphs
grid off % removes graph grid
box on % shows all the borders of the graph instead of just a few
set(gca, 'layer', 'top') % puts all labels and axes (borders) on top of the plotted data
set(gca, 'linewidth', 2) % controls the linewidth of the graph borders
set(gca, 'fontsize', 18) % sets the font size for all elements of the plot
set(gca, 'fontname', 'times new roman') % sets the font type for all elements of the plot
set(gcf, 'windowstate', 'maximized') % maximizes the figure window
set(gcf, 'color', 'w'); % sets the background color behind the plot ('w': white)
pbaspect([1 1 1]) % ratio of x, y, z axes. [1 1 1] for square plots
% to include and modify a colorbar, uncomment the following line:
colormap(flipud(jet));
xlabel('Distance [nm]', 'FontSize', 21);
ylabel('Bias voltage [mV]', 'FontSize', 21);
colorbar
axis tight

% Visualization of the fit with error bars included
figure(4);
plot(yfit_complete);
hold on
errorbar(x_real, y_vals, err, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 2);
xlabel('Distance (arb. units)', 'FontSize', 20);
ylabel('$dI/dV$ (arb. units)', 'Interpreter', 'latex', 'FontSize', 20);
legend('MBS probability density fit', 'Error', 'FontSize', 19, 'Orientation', 'vertical', 'Location', 'best');
grid on
hold off







































