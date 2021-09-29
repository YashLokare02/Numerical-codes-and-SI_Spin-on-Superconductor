clc;
clear;
close all;

load('2K.mat');
load('Delta.mat');


spectra_data1 = load('1.08K.txt');
spectra_data2 = load('1.25K.txt');
spectra_data3 = load('1.5K.txt');
spectra_data4 = load('1.75K.txt');
spectra_data5 = load('2K.txt');
spectra_data6 = load('2.25K.txt');
spectra_data7 = load('2.5K.txt');
spectra_data8 = load('2.75K.txt');
spectra_data9 = load('3K.txt');
spectra_data10 = load('3.25K.txt');
spectra_data11 = load('3.5K.txt');
spectra_data12 = load('3.75K.txt');
spectra_data13 = load('4K.txt');

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

%%%====================================== 
T = [1.08 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4]; % Vector containing the temperatures at which the STM data was collected 
kb = 1.3807*1.0e-4; % Boltzmann's constant (in units of eV)
e = 1; 



%  % Vector containing the fitting parameter values (order: [\Gamma \Delta])
cond = zeros(13, 193);
% 

for i = 1:length(T)
    
    %BCS_sol(i, :) = (@(para_val, V) BCS_fit(para_val(1), para_val(2), V), x0, volt_data(i, :), sigma_val(i, :), [], [], options);
    delta = delta_mv(i)*0.001;
    tao = Gamma_mv(i)*0.001;
    val = val_mv(i)*0.001;
    cond(i, :) = abs(real(complex(x,tao)./sqrt(complex(x,tao).^2-delta^2))).*(1+cosh((e*x)./2*kb*T(i))).^-1;
    z = cond(i,:);
    tmp=x-val;
    [idx idx]= min(abs(tmp));                       %Normalization
    norm = y(idx);
    cond(i,:)=cond(i,:)/z(idx)*norm;
        
end

figure(1)
axes()
hold on
plot(x,cond(1,:),'b');
hold on
%plot(volt_data(1, :), sigma_val(1, :), '*r');
%hold on
plot(x,cond(2,:),'b');
hold on
%plot(volt_data(2, :), sigma_val(2, :), '*r');
%hold on
plot(x,cond(3,:),'b');
hold on
%plot(volt_data(3, :), sigma_val(3, :), '*r');
%hold on
plot(x,cond(4,:),'b');
hold on
%plot(volt_data(4, :), sigma_val(4, :), '*r');
%hold on
plot(x,cond(5,:),'b');
hold on
%plot(volt_data(5, :), sigma_val(5, :), '*r');
%hold on
plot(x,cond(6,:),'b');
hold on
%plot(volt_data(6, :), sigma_val(6, :), '*r');
%hold on
plot(x,cond(7,:),'b');
hold on
%plot(volt_data(7, :), sigma_val(7, :), '*r');
%hold on
plot(x,cond(8,:),'b');
hold on
%plot(volt_data(8, :), sigma_val(8, :), '*r');
%hold on
plot(x,cond(9,:),'b');
hold on
%plot(volt_data(9, :), sigma_val(9, :), '*r');
%hold on
plot(x,cond(10,:),'b');
hold on
%plot(volt_data(10, :), sigma_val(10, :), '*r');
%hold on
plot(x,cond(11,:),'b');
hold on
%plot(volt_data(11, :), sigma_val(11, :), '*r');
%hold on
plot(x,cond(12,:),'b');
hold on
%plot(volt_data(12, :), sigma_val(12, :), '*r');
%hold on
plot(x,cond(13,:),'b');
hold on
%plot(volt_data(13, :), sigma_val(13, :), '*r');

title('BCS Dynes DOS fit to the STM tunneling spectra', 'FontSize', 20);
xlabel('Voltage [V]', 'FontSize', 18);
ylabel('$\sigma(V) \equiv dI/dV$', 'Interpreter', 'latex', 'FontSize', 18);
legend('FontSize', 16, 'Location', 'best');
hold off
