%% Fitting the BCS Dynes density of states to the tunneling spectra data (obtained using STM)
% Aim: To estimate the superconducting energy gap and the thermal broadening parameter from the generated fit

spectra_data = load('Dataset I-V curve 1.08 K.txt');
sigma_val = spectra_data(:, 1); % Vector containing the differential oonductance values (normalized)
volt_data = spectra_data(:, 2); % Vector containing the applied bias voltage values

%% Visualization of the imported dataset
figure(1);
hold on
plot(sigma_val, volt_data, 'bo');
title('Plot for the normalized tunneling spectra data vs. applied bias voltage', 'FontSize', 20);
xlabel('Voltage [V]', 'FontSize', 18);
ylabel('$\sigma(V) \equiv dI/dV$', 'Interpreter', 'latex', 'FontSize', 18);
hold off

%% Fitting the BCS Dynes density of states to the STM tunneling spectra
% Generating the numerical fit
para_fit = lsqcurvefit(@(c, xdata) BCS_fit(c(1), c(2), xdata), [0.65 0.15], volt_data, sigma_val);

%% Visualization of the fit
figure(2)
plot(sigma_val, volt_data, 'bo');
hold on
plot(volt_data, BCS_fit(para_fit(1), para_fit(2), volt_data), 'b') % Plotting the BCS fit as a smooth curve
title('Plot for the normalized tunneling spectra data vs. applied bias voltage', 'FontSize', 20);
xlabel('Voltage [V]', 'FontSize', 18);
ylabel('$\sigma(V) \equiv dI/dV$', 'Interpreter', 'latex', 'FontSize', 18);
hold off



