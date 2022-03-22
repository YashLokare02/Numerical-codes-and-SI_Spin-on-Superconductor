%% Generation of the superconducting energy gap map

% Bias voltage values
[volt] = Nanonis_reader('Gap map1.3ds', 1); 
volt_data = volt(:, 3); %% Bias voltage values [V]
volt_data = volt_data.*1e3;

% Paramters (no. of data points)
n = 3600; %% Number of grid points in total
rows = 192; %% No. of bias voltage values (sweep range)
cols = n*3; %% Number of columns in the data matrix
k = 60; %% Number of grid points along each row and column of the scan range

% Initialization of the data matrix
data = zeros(rows, cols);

data(:, 1:3) = Nanonis_reader('Gap map1.3ds', 0);

% Data acquisition
for i = 1:n-1
    data(:, 3*i+1:3*i+3) = Nanonis_reader('Gap map1.3ds', i);
end

% STM spectral data symmetrization and normalization
for i = 2:3:cols-1
    data((rows/2)+1:rows, i) = data(rows/2:-1:1, i); %% Symmetrization
    data(:, i) = (data(:, i) - min(data(:, i)))./(max(data(:, i)) - min(data(:, i))); %% Normalization
end

% Storing the symmetrized STM data in a separate matrix
data_sym = zeros(rows, n);

data_sym(:, 1) = data(:, 2);

for i = 2:n
    data_sym(:, i) = data(:, 3*i-1);
end

% Scan area dimensions
x = linspace(0, 20, k);
y = linspace(0, 20, k);

% Sorting volt_data to locate peaks and compute their widths using findpeaks()
volt_val = sort(volt_data, 'ascend');

% Extracting the relevant data
volt_val1 = volt_val(71:122);

no_vals = 122 - 71 + 1;

data_gap = zeros(no_vals, n);

for i = 1:n
    data_gap(:, i) = data_sym(71:122, i);
end

% Smoothening of data
for i = 1:n
    data_gap(:, i) = smooth(data_gap(:, i));
end

% Initializing peak_widths
peak_widths = zeros(5, n); % each row contains superconducting gap values corresponding to a certain normalization % value
norm_vals = [0.77, 0.7, 0.65, 0.6, 0.5]; % normalization % values

% Superconducting energy gap estimation
for j = 1:length(norm_vals)
    for i = 1:n
        [x_int, y_int, iout, jout] = intersections(volt_val1', (data_gap(:, i))', [volt_val1(1) volt_val1(end)], [norm_vals(j)*max(data_gap(:, i)) norm_vals(j)*max(data_gap(:, i))]);
    if length(x_int) == 8
        peak_widths(j, i) = x_int(5) - x_int(4);
    elseif length(x_int) == 6
        peak_widths(j, i) = x_int(4) - x_int(3);
    elseif length(x_int) == 4
        peak_widths(j, i) = x_int(3) - x_int(2);
    else
        peak_widths(j, i) = x_int(end) - x_int(1); % for length(x_int) == 2
    end
    peak_widths(j, i) = peak_widths(j, i)/2;
    end
end

% Calibration of \Delta values (due to imperfections in the STM measurements)
peak_widths = 0.7*peak_widths; % in meV

% Visualization (2D color maps for \Delta, corresponding to each normalization slice)
% For 77% normalization slice
gap_matrix_1 = reshape(peak_widths(1, :), k, k);
gap_vals_1 = gap_matrix_1'; % aligning along the corresponding grid points

% For 70% normalization slice
gap_matrix_2 = reshape(peak_widths(2, :), k, k);
gap_vals_2 = gap_matrix_2'; % aligning along the corresponding grid points

% For 65% normalization slice
gap_matrix_3 = reshape(peak_widths(3, :), k, k);
gap_vals_3 = gap_matrix_3'; % aligning along the corresponding grid points

% For 60% normalization slice
gap_matrix_4 = reshape(peak_widths(4, :), k, k);
gap_vals_4 = gap_matrix_4'; % aligning along the corresponding grid points

% For 50% normalization slice
gap_matrix_5 = reshape(peak_widths(5, :), k, k);
gap_vals_5 = gap_matrix_5'; % aligning along the corresponding grid points

% Figures
figure(1);
surf(x, y, gap_vals_1)
brighten(0.6);
view([0 0 90])
hold on % holds the opened figure for further modification
shading interp % interpolates datapoints, ideal for 3D graphs
grid off % removes graph grid
box on % shows all the borders of the graph instead of just a few
set(gca, 'layer', 'top') % puts all labels and axes (borders) on top of the plotted data
set(gca, 'linewidth', 2) % controls the linewidth of the graph borders
set(gca, 'fontsize', 15) % sets the font size for all elements of the plot
set(gca, 'fontname', 'times new roman') % sets the font type for all elements of the plot
set(gcf, 'windowstate', 'maximized') % maximizes the figure window
set(gcf, 'color', 'w'); % sets the background color behind the plot ('w': white)
pbaspect([1 1 1]) % ratio of x, y, z axes. [1 1 1] for square plots
% to include and modify a colorbar uncomment the following line:
colormap(jet(4096));
colorbar
xlabel('Distance along x [nm]', 'FontSize', 18);
ylabel('Distance along y [nm]', 'FontSize', 18);
title('2D color map for $\Delta$ (meV) [at 77$\%$ normalization slice]', 'Interpreter', 'latex', 'FontSize', 17);

figure(2);
surf(x, y, gap_vals_2)
brighten(0.6);
view([0 0 90])
hold on % holds the opened figure for further modification
shading interp % interpolates datapoints, ideal for 3D graphs
grid off % removes graph grid
box on % shows all the borders of the graph instead of just a few
set(gca, 'layer', 'top') % puts all labels and axes (borders) on top of the plotted data
set(gca, 'linewidth', 2) % controls the linewidth of the graph borders
set(gca, 'fontsize', 15) % sets the font size for all elements of the plot
set(gca, 'fontname', 'times new roman') % sets the font type for all elements of the plot
set(gcf, 'windowstate', 'maximized') % maximizes the figure window
set(gcf, 'color', 'w'); % sets the background color behind the plot ('w': white)
pbaspect([1 1 1]) % ratio of x, y, z axes. [1 1 1] for square plots
% to include and modify a colorbar uncomment the following line:
colormap(jet(4096));
colorbar
xlabel('Distance along x [nm]', 'FontSize', 18);
ylabel('Distance along y [nm]', 'FontSize', 18);
title('2D color map for $\Delta$ (meV) [at 70$\%$ normalization slice]', 'Interpreter', 'latex', 'FontSize', 17);

figure(3);
surf(x, y, gap_vals_3)
brighten(0.6);
view([0 0 90])
hold on % holds the opened figure for further modification
shading interp % interpolates datapoints, ideal for 3D graphs
grid off % removes graph grid
box on % shows all the borders of the graph instead of just a few
set(gca, 'layer', 'top') % puts all labels and axes (borders) on top of the plotted data
set(gca, 'linewidth', 2) % controls the linewidth of the graph borders
set(gca, 'fontsize', 15) % sets the font size for all elements of the plot
set(gca, 'fontname', 'times new roman') % sets the font type for all elements of the plot
set(gcf, 'windowstate', 'maximized') % maximizes the figure window
set(gcf, 'color', 'w'); % sets the background color behind the plot ('w': white)
pbaspect([1 1 1]) % ratio of x, y, z axes. [1 1 1] for square plots
% to include and modify a colorbar uncomment the following line:
colormap(jet(4096));
colorbar
xlabel('Distance along x [nm]', 'FontSize', 18);
ylabel('Distance along y [nm]', 'FontSize', 18);
title('2D color map for $\Delta$ (meV) [at 65$\%$ normalization slice]', 'Interpreter', 'latex', 'FontSize', 17);

figure(4);
surf(x, y, gap_vals_4)
brighten(0.6);
view([0 0 90])
hold on % holds the opened figure for further modification
shading interp % interpolates datapoints, ideal for 3D graphs
grid off % removes graph grid
box on % shows all the borders of the graph instead of just a few
set(gca, 'layer', 'top') % puts all labels and axes (borders) on top of the plotted data
set(gca, 'linewidth', 2) % controls the linewidth of the graph borders
set(gca, 'fontsize', 15) % sets the font size for all elements of the plot
set(gca, 'fontname', 'times new roman') % sets the font type for all elements of the plot
set(gcf, 'windowstate', 'maximized') % maximizes the figure window
set(gcf, 'color', 'w'); % sets the background color behind the plot ('w': white)
pbaspect([1 1 1]) % ratio of x, y, z axes. [1 1 1] for square plots
% to include and modify a colorbar uncomment the following line:
colormap(jet(4096));
colorbar
xlabel('Distance along x [nm]', 'FontSize', 18);
ylabel('Distance along y [nm]', 'FontSize', 18);
title('2D color map for $\Delta$ (meV) [at 60$\%$ normalization slice]', 'Interpreter', 'latex', 'FontSize', 17);

figure(5);
surf(x, y, gap_vals_5)
brighten(0.6);
view([0 0 90])
hold on % holds the opened figure for further modification
shading interp % interpolates datapoints, ideal for 3D graphs
grid off % removes graph grid
box on % shows all the borders of the graph instead of just a few
set(gca, 'layer', 'top') % puts all labels and axes (borders) on top of the plotted data
set(gca, 'linewidth', 2) % controls the linewidth of the graph borders
set(gca, 'fontsize', 15) % sets the font size for all elements of the plot
set(gca, 'fontname', 'times new roman') % sets the font type for all elements of the plot
set(gcf, 'windowstate', 'maximized') % maximizes the figure window
set(gcf, 'color', 'w'); % sets the background color behind the plot ('w': white)
pbaspect([1 1 1]) % ratio of x, y, z axes. [1 1 1] for square plots
% to include and modify a colorbar uncomment the following line:
colormap(jet(4096));
colorbar
xlabel('Distance along x [nm]', 'FontSize', 18);
ylabel('Distance along y [nm]', 'FontSize', 18);
title('2D color map for $\Delta$ (meV) [at 50$\%$ normalization slice]', 'Interpreter', 'latex', 'FontSize', 17);

% Histogram distributions (2D) of \Delta at different normalization slices
% Computing bin widths for each histogram distribution
no_bins = sqrt(n);
rounded_no_bins = round(no_bins, 0);
bin_width = zeros(1, length(norm_vals));

for i = 1:length(norm_vals)
    range = max(peak_widths(i, :)) - min(peak_widths(i, :));
    bin_width(i) = range/rounded_no_bins;
end

% Visualization
figure(6); % at 77% normalization slice
histogram(peak_widths(1, :), 'FaceAlpha', 0.1, 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', bin_width(1))
title('2D histogram distribution of $\Delta$ (at 77$\%$ normalization slice)', 'Interpreter', 'latex', 'FontSize', 17);
xlabel('$\Delta$ [meV]', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Probability', 'FontSize', 18);

figure(7); % at 70% normalization slice
histogram(peak_widths(2, :), 'FaceAlpha', 0.1, 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', bin_width(2))
title('2D histogram distribution of $\Delta$ (at 70$\%$ normalization slice)', 'Interpreter', 'latex', 'FontSize', 17);
xlabel('$\Delta$ [meV]', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Probability', 'FontSize', 18);

figure(8); % at 65% normalization slice
histogram(peak_widths(3, :), 'FaceAlpha', 0.1, 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', bin_width(3))
title('2D histogram distribution of $\Delta$ (at 65$\%$ normalization slice)', 'Interpreter', 'latex', 'FontSize', 17);
xlabel('$\Delta$ [meV]', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Probability', 'FontSize', 18);

figure(9); % at 60% normalization slice
histogram(peak_widths(4, :), 'FaceAlpha', 0.1, 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', bin_width(4))
title('2D histogram distribution of $\Delta$ (at 60$\%$ normalization slice)', 'Interpreter', 'latex', 'FontSize', 17);
xlabel('$\Delta$ [meV]', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Probability', 'FontSize', 18);

figure(10); % at 50% normalization slice
histogram(peak_widths(5, :), 'FaceAlpha', 0.1, 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', bin_width(5))
title('2D histogram distribution of $\Delta$ (at 50$\%$ normalization slice)', 'Interpreter', 'latex', 'FontSize', 17);
xlabel('$\Delta$ [meV]', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Probability', 'FontSize', 18);

% Numerical estimation of 2D autocorrelation functions for \Delta values at different normalization slices
% For 77% normalization slice values
auto2d_77 = autocorr2d(gap_vals_1);

% For 70% normalization slice values
auto2d_70 = autocorr2d(gap_vals_2);

% For 65% normalization slice values
auto2d_65 = autocorr2d(gap_vals_3);

% For 60% normalization slice values
auto2d_60 = autocorr2d(gap_vals_4);

% For 50% normalization slice values
auto2d_50 = autocorr2d(gap_vals_5);

% Visualization of the autocorrelation maps for \Delta values at all normalization slices
figure(11); % at 77% normalization slice
surf(auto2d_77)
brighten(0.6);
view([0 0 90])
hold on % holds the opened figure for further modification
shading interp % interpolates datapoints, ideal for 3D graphs
grid off % removes graph grid
box on % shows all the borders of the graph instead of just a few
set(gca, 'layer', 'top') % puts all labels and axes (borders) on top of the plotted data
set(gca, 'linewidth', 2) % controls the linewidth of the graph borders
set(gca, 'fontsize', 16) % sets the font size for all elements of the plot
set(gca, 'fontname', 'times new roman') % sets the font type for all elements of the plot
set(gcf, 'windowstate', 'maximized') % maximizes the figure window
set(gcf, 'color', 'w'); % sets the background color behind the plot ('w': white)
pbaspect([1 1 1]) % ratio of x, y, z axes. [1 1 1] for square plots
% to include and modify a colorbar uncomment the following line:
colormap(flipud(jet));
colorbar
xlabel('Grid points along x', 'FontSize', 23);
ylabel('Grid points along y', 'FontSize', 23);
title('2D autocorrelation function (77$\%$ normalization slice)', 'Interpreter', 'latex', 'FontSize', 17);
axis tight

figure(12); % at 70% normalization slice
surf(auto2d_70)
brighten(0.6);
view([0 0 90])
hold on % holds the opened figure for further modification
shading interp % interpolates datapoints, ideal for 3D graphs
grid off % removes graph grid
box on % shows all the borders of the graph instead of just a few
set(gca, 'layer', 'top') % puts all labels and axes (borders) on top of the plotted data
set(gca, 'linewidth', 2) % controls the linewidth of the graph borders
set(gca, 'fontsize', 16) % sets the font size for all elements of the plot
set(gca, 'fontname', 'times new roman') % sets the font type for all elements of the plot
set(gcf, 'windowstate', 'maximized') % maximizes the figure window
set(gcf, 'color', 'w'); % sets the background color behind the plot ('w': white)
pbaspect([1 1 1]) % ratio of x, y, z axes. [1 1 1] for square plots
% to include and modify a colorbar uncomment the following line:
colormap(flipud(jet));
colorbar
xlabel('Grid points along x', 'FontSize', 23);
ylabel('Grid points along y', 'FontSize', 23);
title('2D autocorrelation function (70$\%$ normalization slice)', 'Interpreter', 'latex', 'FontSize', 17);
axis tight

figure(13); % at 65% normalization slice
surf(auto2d_65)
brighten(0.6);
view([0 0 90])
hold on % holds the opened figure for further modification
shading interp % interpolates datapoints, ideal for 3D graphs
grid off % removes graph grid
box on % shows all the borders of the graph instead of just a few
set(gca, 'layer', 'top') % puts all labels and axes (borders) on top of the plotted data
set(gca, 'linewidth', 2) % controls the linewidth of the graph borders
set(gca, 'fontsize', 16) % sets the font size for all elements of the plot
set(gca, 'fontname', 'times new roman') % sets the font type for all elements of the plot
set(gcf, 'windowstate', 'maximized') % maximizes the figure window
set(gcf, 'color', 'w'); % sets the background color behind the plot ('w': white)
pbaspect([1 1 1]) % ratio of x, y, z axes. [1 1 1] for square plots
% to include and modify a colorbar uncomment the following line:
colormap(flipud(jet));
colorbar
xlabel('Grid points along x', 'FontSize', 23);
ylabel('Grid points along y', 'FontSize', 23);
title('2D autocorrelation function (65$\%$ normalization slice)', 'Interpreter', 'latex', 'FontSize', 17);
axis tight

figure(14); % at 77% normalization slice
surf(auto2d_60)
brighten(0.6);
view([0 0 90])
hold on % holds the opened figure for further modification
shading interp % interpolates datapoints, ideal for 3D graphs
grid off % removes graph grid
box on % shows all the borders of the graph instead of just a few
set(gca, 'layer', 'top') % puts all labels and axes (borders) on top of the plotted data
set(gca, 'linewidth', 2) % controls the linewidth of the graph borders
set(gca, 'fontsize', 16) % sets the font size for all elements of the plot
set(gca, 'fontname', 'times new roman') % sets the font type for all elements of the plot
set(gcf, 'windowstate', 'maximized') % maximizes the figure window
set(gcf, 'color', 'w'); % sets the background color behind the plot ('w': white)
pbaspect([1 1 1]) % ratio of x, y, z axes. [1 1 1] for square plots
% to include and modify a colorbar uncomment the following line:
colormap(flipud(jet));
colorbar
xlabel('Grid points along x', 'FontSize', 23);
ylabel('Grid points along y', 'FontSize', 23);
title('2D autocorrelation function (60$\%$ normalization slice)', 'Interpreter', 'latex', 'FontSize', 17);
axis tight

figure(15); % at 77% normalization slice
surf(auto2d_50)
brighten(0.6);
view([0 0 90])
hold on % holds the opened figure for further modification
shading interp % interpolates datapoints, ideal for 3D graphs
grid off % removes graph grid
box on % shows all the borders of the graph instead of just a few
set(gca, 'layer', 'top') % puts all labels and axes (borders) on top of the plotted data
set(gca, 'linewidth', 2) % controls the linewidth of the graph borders
set(gca, 'fontsize', 16) % sets the font size for all elements of the plot
set(gca, 'fontname', 'times new roman') % sets the font type for all elements of the plot
set(gcf, 'windowstate', 'maximized') % maximizes the figure window
set(gcf, 'color', 'w'); % sets the background color behind the plot ('w': white)
pbaspect([1 1 1]) % ratio of x, y, z axes. [1 1 1] for square plots
% to include and modify a colorbar uncomment the following line:
colormap(flipud(jet));
colorbar
xlabel('Grid points along x', 'FontSize', 23);
ylabel('Grid points along y', 'FontSize', 23);
title('2D autocorrelation function (50$\%$ normalization slice)', 'Interpreter', 'latex', 'FontSize', 17);
axis tight

% Superconducting coherence length estimation
% Initializing line scan ranges
x_half = x((k/2)+1:k); % along x
y_half = y((k/2)+1:k); % along y

% For 77% normalization slice values
% Along x
x_correlation_77 = auto2d_77(k/2, (k/2)+1:k);

[xi_x_77, para_vals_x_77] = exp_fit(x_half, x_correlation_77, 0.17, 0.9);

Along y
y_correlation_77 = auto2d_77(1:k/2, k/2);
y_correlation_77 = flip(y_correlation_77);
y_correlation_77 = y_correlation_77';

[xi_y_77, para_vals_y_77] = exp_fit(y_half, y_correlation_77, 0.17, 0.9);

% For 70% normalization slice values
Along x
x_correlation_70 = auto2d_70(k/2, (k/2)+1:k);

[xi_x_70, para_vals_x_70] = exp_fit(x_half, x_correlation_70, 0.17, 0.9);

% Along y
y_correlation_70 = auto2d_70(1:k/2, k/2);
y_correlation_70 = flip(y_correlation_70);
y_correlation_70 = y_correlation_70';

[xi_y_70, para_vals_y_70] = exp_fit(y_half, y_correlation_70, 0.17, 0.9);

% For 65% normalization slice values
% Along x
x_correlation_65 = auto2d_65(k/2, (k/2)+1:k);

[xi_x_65, para_vals_x_65] = exp_fit(x_half, x_correlation_65, 0.17, 0.9);

% Along y
y_correlation_65 = auto2d_65(1:k/2, k/2);
y_correlation_65 = flip(y_correlation_65);
y_correlation_65 = y_correlation_65';

[xi_y_65, para_vals_y_65] = exp_fit(y_half, y_correlation_65, 0.17, 0.9);

% For 60% normalization slice values
% Along x
x_correlation_60 = auto2d_60(k/2, (k/2)+1:k);

[xi_x_60, para_vals_x_60] = exp_fit(x_half, x_correlation_60, 0.17, 0.9);

% Along y
y_correlation_60 = auto2d_60(1:k/2, k/2);
y_correlation_60 = flip(y_correlation_60);
y_correlation_60 = y_correlation_60';

[xi_y_60, para_vals_y_60] = exp_fit(y_half, y_correlation_60, 0.17, 0.9);

% For 50% normalization slice values
% Along x
x_correlation_50 = auto2d_50(k/2, (k/2)+1:k);

[xi_x_50, para_vals_x_50] = exp_fit(x_half, x_correlation_50, 0.17, 0.9);

% Along y
y_correlation_50 = auto2d_50(1:k/2, k/2);
y_correlation_50 = flip(y_correlation_50);
y_correlation_50 = y_correlation_50';

[xi_y_50, para_vals_y_50] = exp_fit(y_half, y_correlation_50, 0.17, 0.9);

% Visualization of the dI/dV spectra (line scan plots)
% Considered scan area -- 10 x 10 nm^2
% Scan length scales
y1 = linspace(0, 10, 10);
x1 = y1;
% Extracting spectral datasets along y (taken in the middle of the 10 x 10 nm^2 scan range)
% Normalizing bias voltage values
volt_data = volt_data.*0.7; % in mV

% Initialization of the reduced data matrices
data_along_y = data_sym(:, 15:61:564); % along y
data_along_x = data_sym(:, 854:61:1403); % along x

figure(16); % along y
surf(y1, volt_data, data_along_y)
view([0 0 90])
ylim([-2 2])
hold on % holds the opened figure for further modification
shading interp % interpolates datapoints, ideal for 3D graphs
grid off % removes graph grid
box on % shows all the borders of the graph instead of just a few
set(gca, 'layer', 'top') % puts all labels and axes (borders) on top of the plotted data
set(gca, 'linewidth', 2) % controls the linewidth of the graph borders
set(gca, 'fontsize', 18) % sets the font size for all elements of the plot
set(gca, 'fontname', 'times new roman') % sets the font type for all elements of the plot
set(gcf, 'color', 'w'); % sets the background color behind the plot ('w': white)
% to include and modify a colorbar uncomment the following line:
colormap(flipud(jet));
colorbar
title('10 nm vertical line scan', 'FontSize', 19);

figure(17); % along x
surf(x1, volt_data, data_along_x)
view([0 0 90])
ylim([-2 2])
hold on % holds the opened figure for further modification
shading interp % interpolates datapoints, ideal for 3D graphs
grid off % removes graph grid
box on % shows all the borders of the graph instead of just a few
set(gca, 'layer', 'top') % puts all labels and axes (borders) on top of the plotted data
set(gca, 'linewidth', 2) % controls the linewidth of the graph borders
set(gca, 'fontsize', 18) % sets the font size for all elements of the plot
set(gca, 'fontname', 'times new roman') % sets the font type for all elements of the plot
set(gcf, 'color', 'w'); % sets the background color behind the plot ('w': white)
% to include and modify a colorbar uncomment the following line:
colormap(flipud(jet));
colorbar
title('10 nm horizontal line scan', 'FontSize', 19);

% Visualization of the same color maps in 3D (waterfall plots)
figure(18); % along y
surf(y1, volt_data, data_along_y)
ylim([-5 5])
hold on % holds the opened figure for further modification
shading interp % interpolates datapoints, ideal for 3D graphs
grid off % removes graph grid
box on % shows all the borders of the graph instead of just a few
set(gca, 'linewidth', 2) % controls the linewidth of the graph borders
set(gca, 'layer', 'top') % puts all labels and axes (borders) on top of the plotted data
set(gca, 'fontsize', 18) % sets the font size for all elements of the plot
set(gca, 'fontname', 'times new roman') % sets the font type for all elements of the plot
set(gcf, 'color', 'w'); % sets the background color behind the plot ('w': white)
set(gca, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none')
% to include and modify a colorbar uncomment the following line:
colormap(flipud(jet));
colorbar
title('10 nm vertical line scan', 'FontSize', 19);

figure(19); % along x
surf(x1, volt_data, data_along_x)
ylim([-5 5])
hold on % holds the opened figure for further modification
shading interp % interpolates datapoints, ideal for 3D graphs
grid off % removes graph grid
box on % shows all the borders of the graph instead of just a few
set(gca, 'linewidth', 2) % controls the linewidth of the graph borders
set(gca, 'layer', 'top') % puts all labels and axes (borders) on top of the plotted data
set(gca, 'fontsize', 18) % sets the font size for all elements of the plot
set(gca, 'fontname', 'times new roman') % sets the font type for all elements of the plot
set(gcf, 'color', 'w'); % sets the background color behind the plot ('w': white)
set(gca, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none')
% to include and modify a colorbar uncomment the following line:
colormap(flipud(jet));
colorbar
title('10 nm horizontal line scan', 'FontSize', 19);













    



