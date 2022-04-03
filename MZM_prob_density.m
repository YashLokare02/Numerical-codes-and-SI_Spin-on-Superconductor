function [fitresult, gof] = MZM_prob_density(x_real, data_trial)

%      X input: x_actual
%      Y output: ZBCP_actual
%  Output:
%      fitresult: a fit object representing the fit
%      gof: goodness-of-fit information
%

%% Fit: MBS probability density fit to the ZBCP data
[xData, yData] = prepareCurveData(x_real, data_trial);

% Set up fittype and options
ft = fittype('(1./a.*(exp(-2.*x./b) + exp(-2.*x./c) - 2.*exp(-x.*(1./b + 1./c)).*cos(2.*x.*d)))', 'independent', 'x', 'dependent', 'y');
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.StartPoint = [0.692202047570882 0.226700063401286 0.766536307847191 60];

% Fit model to data
[fitresult, gof] = fit(xData, yData, ft, opts);

% % Plot fit with data
figure('Name', 'MBS probability density fit to the ZBCP data');
h = plot(fitresult, xData, yData);
ylim([-0.1 3.5])
hold on
legend(h, 'Differential conductance data', 'MBS probability density fit', 'FontSize', 17, 'Location', 'best', 'Interpreter', 'none');
% Label axes
xlabel('Distance (arb. units)', 'Interpreter', 'none', 'FontSize', 20);
ylabel('$dI/dV$ (arb. units)', 'Interpreter', 'latex', 'FontSize', 20);
grid on
