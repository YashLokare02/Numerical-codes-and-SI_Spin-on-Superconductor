%% Fitting the BCS Dynes density of states to the tunneling spectra data (obtained using STM)
% Aim: To estimate the superconducting energy gap and the thermal broadening parameter from the generated fit
% Note: All values have been considered and/or computed in units of eV

% Importing the STM tunneling spectra into MATLAB
spectra_data1 = load('Dataset I-V curve 1.08 K.txt');
spectra_data2 = load('Dataset I-V curve 1.25 K.txt');
spectra_data3 = load('Dataset I-V curve 1.5 K.txt');
spectra_data4 = load('Dataset I-V curve 1.75 K.txt');
spectra_data5 = load('Dataset I-V curve 2 K.txt');
spectra_data6 = load('Dataset I-V curve 2.25 K.txt');
spectra_data7 = load('Dataset I-V curve 2.5 K.txt');
spectra_data8 = load('Dataset I-V curve 2.75 K.txt');
spectra_data9 = load('Dataset I-V curve 3 K.txt');
spectra_data10 = load('Dataset I-V curve 3.25 K.txt');
spectra_data11 = load('Dataset I-V curve 3.5 K.txt');
spectra_data12 = load('Dataset I-V curve 3.75 K.txt');
spectra_data13 = load('Dataset I-V curve 4 K.txt');

% Initializing data storage vectors + storing the STM tunneling spectra data
volt_data = zeros(13, length(spectra_data1));
sigma_val = zeros(13, length(spectra_data1));

% For $T = 1.08 K$
volt_data(1, :) = spectra_data1(:, 1)';
sigma_val(1, :) = spectra_data1(:, 2)';

% For $T = 1.25 K$
volt_data(2, :) = spectra_data2(:, 1)';
sigma_val(2, :) = spectra_data2(:, 2)';

% For $T = 1.5 K$
volt_data(3, :) = spectra_data3(:, 1)';
sigma_val(3, :) = spectra_data3(:, 2)';

% For $T = 1.75 K$
volt_data(4, :) = spectra_data4(:, 1)';
sigma_val(4, :) = spectra_data4(:, 2)';

% For $T = 2 K$
volt_data(5, :) = spectra_data5(:, 1)';
sigma_val(5, :) = spectra_data5(:, 2)';

% For %T = 2.25 K$
volt_data(6, :) = spectra_data6(:, 1)';
sigma_val(6, :) = spectra_data6(:, 2)';

% For $T = 2.5 K$
volt_data(7, :) = spectra_data7(:, 1)';
sigma_val(7, :) = spectra_data7(:, 2)';

% For $T = 2.75 K$
volt_data(8, :) = spectra_data8(:, 1)';
sigma_val(8, :) = spectra_data8(:, 2)';

% For %T = 3 K$
volt_data(9, :) = spectra_data9(:, 1)';
sigma_val(9, :) = spectra_data9(:, 2)';

% For $T = 3.25 K$
volt_data(10, :) = spectra_data10(:, 1)';
sigma_val(10, :) = spectra_data10(:, 2)';

% For $T = 3.5 K$
volt_data(11, :) = spectra_data11(:, 1)';
sigma_val(11, :) = spectra_data11(:, 2)';

% For $T = 3.75 K$
volt_data(12, :) = spectra_data12(:, 1)';
sigma_val(12, :) = spectra_data12(:, 2)';

% For $T = 4 K$
volt_data(13, :) = spectra_data13(:, 1)';
sigma_val(13, :) = spectra_data13(:, 2)';

% Defining parameters
T = [1.08 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4]; % Vector containing the temperatures at which the STM data was collected 
k = (1.38/1.6)*1.0e-4; % Boltzmann's constant (in units of eV)
e = 1; % Electronic charge (in units of eV)
Ef = 0; % We take the Fermi energy level (reference point) as zero

%% Generating the BCS Dynes density of states fit
% Note: For this, we make use of the BCS Dynes density of states expression

options = optimoptions('lsqcurvefit', 'Algorithm', 'trust-region-reflective', 'OptimalityTolerance', 1e-6, 'FunctionTolerance', 1e-6); % Setting optimization options
x0 = [0.45*1.0e-4 1e-5]; % Initial guesses for the parametric fit
BCS_sol = zeros(13, 2); % Vector containing the fitting parameter values (order: [\Gamma \Delta])
gamma = zeros(13, 1); % Vector containing the thermal broadening parameter values
delta = zeros(13, 1); % Vector containing the superconducting energy gap values

% Numerically estimating the fitting parameters

for i = 1:length(T)
    
    BCS_fit = @(gamma, delta, tt) arrayfun(@(V) integral(@(E) real((E - 1i*gamma)./sqrt((E - 1i*gamma).^2 - delta.^2)).*(cosh((E + e*V)./2*k.*T(i))).^-2, Ef, Ef + e*V), tt);
    BCS_sol(i, :) = lsqcurvefit(@(para_val, V) BCS_fit(para_val(1), para_val(2), V), x0, volt_data(i, :), sigma_val(i, :), [], [], options);
    gamma(i) = BCS_sol(i, 1);
    delta(i) = BCS_sol(i, 2);
        
end

%% Visualization of the BCS generated fits
figure(1)
axes();
hold on
plot(volt_data(1, :), sigma_val(1, :), 'bo', 'DisplayName', 'Experimental Data');
hold on
plot(volt_data(1, :), BCS_fit(gamma(1), delta(1), sigma_val(1, :)), 'r', 'DisplayName', 'Estimated Model (BCS fit)');
hold on
plot(volt_data(2, :), sigma_val(2, :), 'bo', 'DisplayName', 'Experimental Data');
hold on
plot(volt_data(2, :), BCS_fit(gamma(2), delta(2), sigma_val(2, :)), 'r', 'DisplayName', 'Estimated Model (BCS fit)');
hold on
plot(volt_data(3, :), sigma_val(3, :), 'bo', 'DisplayName', 'Experimental Data');
hold on
plot(volt_data(3, :), BCS_fit(gamma(3), delta(3), sigma_val(3, :)), 'r', 'DisplayName', 'Estimated Model (BCS fit)');
hold on
plot(volt_data(4, :), sigma_val(4, :), 'bo', 'DisplayName', 'Experimental Data');
hold on
plot(volt_data(4, :), BCS_fit(gamma(4), delta(4), sigma_val(4, :)), 'r', 'DisplayName', 'Estimated Model (BCS fit)');
hold on
plot(volt_data(5, :), sigma_val(5, :), 'bo', 'DisplayName', 'Experimental Data');
hold on
plot(volt_data(5, :), BCS_fit(gamma(5), delta(5), sigma_val(5, :)), 'r', 'DisplayName', 'Estimated Model (BCS fit)');
hold on
plot(volt_data(6, :), sigma_val(6, :), 'bo', 'DisplayName', 'Experimental Data');
hold on
plot(volt_data(6, :), BCS_fit(gamma(6), delta(6), sigma_val(6, :)), 'r', 'DisplayName', 'Estimated Model (BCS fit)');
hold on
plot(volt_data(7, :), sigma_val(7, :), 'bo', 'DisplayName', 'Experimental Data');
hold on
plot(volt_data(7, :), BCS_fit(gamma(7), delta(7), sigma_val(7, :)), 'r', 'DisplayName', 'Estimated Model (BCS fit)');
hold on
plot(volt_data(8, :), sigma_val(8, :), 'bo', 'DisplayName', 'Experimental Data');
hold on
plot(volt_data(8, :), BCS_fit(gamma(8), delta(8), sigma_val(8, :)), 'r', 'DisplayName', 'Estimated Model (BCS fit)');
hold on
plot(volt_data(9, :), sigma_val(9, :), 'bo', 'DisplayName', 'Experimental Data');
hold on
plot(volt_data(9, :), BCS_fit(gamma(9), delta(9), sigma_val(9, :)), 'r', 'DisplayName', 'Estimated Model (BCS fit)');
hold on
plot(volt_data(10, :), sigma_val(10, :), 'bo', 'DisplayName', 'Experimental Data');
hold on
plot(volt_data(10, :), BCS_fit(gamma(10), delta(10), sigma_val(10, :)), 'r', 'DisplayName', 'Estimated Model (BCS fit)');
hold on
plot(volt_data(11, :), sigma_val(11, :), 'bo', 'DisplayName', 'Experimental Data');
hold on
plot(volt_data(11, :), BCS_fit(gamma(11), delta(11), sigma_val(11, :)), 'r', 'DisplayName', 'Estimated Model (BCS fit)');
hold on
plot(volt_data(12, :), sigma_val(12, :), 'bo', 'DisplayName', 'Experimental Data');
hold on
plot(volt_data(12, :), BCS_fit(gamma(12), delta(12), sigma_val(12, :)), 'r', 'DisplayName', 'Estimated Model (BCS fit)');
hold on
plot(volt_data(13, :), sigma_val(13, :), 'bo', 'DisplayName', 'Experimental Data');
hold on
plot(volt_data(13, :), BCS_fit(gamma(13), delta(13), sigma_val(13, :)), 'r', 'DisplayName', 'Estimated Model (BCS fit)');
title('BCS Dynes DOS fit to the STM tunneling spectra', 'FontSize', 20);
xlabel('Voltage [V]', 'FontSize', 18);
ylabel('$\sigma(V) \equiv dI/dV$', 'Interpreter', 'latex', 'FontSize', 18);
legend('FontSize', 16, 'Location', 'best');
hold off

%% Visualization of the fitted parameters vs. temperature
figure(2)
hold on
plot(T, gamma, 'bo');
title('Plot for the thermal broadening parameter vs. temperature', 'FontSize', 20);
xlabel('Temperature [K]', 'FontSize', 18);
ylabel('$\Gamma$', 'Interpreter', 'latex', 'FontSize', 18);
hold off

figure(3)
hold on
plot(T, abs(delta), 'bo');
title('Plot for the superconducting energy gap vs. temperature', 'FontSize', 20);
xlabel('Temperature [K]', 'FontSize', 18);
ylabel('$\Delta$', 'Interpreter', 'latex', 'FontSize', 18);
hold off
