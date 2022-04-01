function [fitresult, gof] = createFit(x_actual, ZBCP_actual)

%      X input: x_actual
%      Y output: ZBCP_actual
%  Output:
%      fitresult: a fit object representing the fit
%      gof: goodness-of-fit info
%

%% Fit: MBS wavefunction fit to the ZBCP data
[xData, yData] = prepareCurveData(x_actual, ZBCP_actual);

% Set up fittype and options
ft = fittype('1./a.*(exp(-2.*x./b) + exp(-2.*x./c) - 2.*exp(-x.*(1./b + 1./c)).*cos(2.*x.*d));', 'independent', 'x', 'dependent', 'y');
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.StartPoint = [0.575045998630839 0.373324874708725 0.0342058755847451 0.913375856139019];

% Fit model to data.
[fitresult, gof] = fit(xData, yData, ft, opts);

% Plot fit with data
figure('Name', 'MBS wavefuntion fit to the ZBCP data');
h = plot(fitresult, xData, yData);
legend(h, 'Differential conductance data', 'MBS wavefuntion fit', 'FontSize', 20, 'Location', 'best', 'Interpreter', 'none');
% Label axes
xlabel('Distance (arb. units)', 'Interpreter', 'none', 'FontSize', 20);
ylabel( '$dI/dV$ (arb. units)', 'Interpreter', 'latex', 'FontSize', 20);
grid on


