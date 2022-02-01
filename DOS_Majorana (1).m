%% Majorana localization -- data analysis
data = zeros(128, 189); %% Matrix containing all the spectroscopic data
k = 0.15; %% Constant vertical shift of all plots (for the full spectra)
k1 = 0.1; %% For the full gapped states

%% Extracting the spectroscopic data at various grid points
[volt] = Nanonis_reader('Majorana localization-01.3ds', 1); 
volt_val = volt(:, 3); %% Bias voltage values [mV]

%% Importing spectroscopic data
[data1] = Nanonis_reader('Majorana localization-01.3ds', 1);
[data2] = Nanonis_reader('Majorana localization-01.3ds', 2);
[data3] = Nanonis_reader('Majorana localization-01.3ds', 3);
[data4] = Nanonis_reader('Majorana localization-01.3ds', 4);
[data5] = Nanonis_reader('Majorana localization-01.3ds', 5);
[data6] = Nanonis_reader('Majorana localization-01.3ds', 6);
[data7] = Nanonis_reader('Majorana localization-01.3ds', 7);
[data8] = Nanonis_reader('Majorana localization-01.3ds', 8);
[data9] = Nanonis_reader('Majorana localization-01.3ds', 9);
[data10] = Nanonis_reader('Majorana localization-01.3ds', 10);
[data11] = Nanonis_reader('Majorana localization-01.3ds', 11);
[data12] = Nanonis_reader('Majorana localization-01.3ds', 12);
[data13] = Nanonis_reader('Majorana localization-01.3ds', 13);
[data14] = Nanonis_reader('Majorana localization-01.3ds', 14);
[data15] = Nanonis_reader('Majorana localization-01.3ds', 15);
[data16] = Nanonis_reader('Majorana localization-01.3ds', 16);
[data17] = Nanonis_reader('Majorana localization-01.3ds', 17);
[data18] = Nanonis_reader('Majorana localization-01.3ds', 18);
[data19] = Nanonis_reader('Majorana localization-01.3ds', 19);
[data20] = Nanonis_reader('Majorana localization-01.3ds', 20);
[data21] = Nanonis_reader('Majorana localization-01.3ds', 21);
[data22] = Nanonis_reader('Majorana localization-01.3ds', 22);
[data23] = Nanonis_reader('Majorana localization-01.3ds', 23);
[data24] = Nanonis_reader('Majorana localization-01.3ds', 24);
[data25] = Nanonis_reader('Majorana localization-01.3ds', 25);
[data26] = Nanonis_reader('Majorana localization-01.3ds', 26);
[data27] = Nanonis_reader('Majorana localization-01.3ds', 27);
[data28] = Nanonis_reader('Majorana localization-01.3ds', 28);
[data29] = Nanonis_reader('Majorana localization-01.3ds', 29);
[data30] = Nanonis_reader('Majorana localization-01.3ds', 30);
[data31] = Nanonis_reader('Majorana localization-01.3ds', 31);
[data32] = Nanonis_reader('Majorana localization-01.3ds', 32);
[data33] = Nanonis_reader('Majorana localization-01.3ds', 33);
[data34] = Nanonis_reader('Majorana localization-01.3ds', 34);
[data35] = Nanonis_reader('Majorana localization-01.3ds', 35);
[data36] = Nanonis_reader('Majorana localization-01.3ds', 36);
[data37] = Nanonis_reader('Majorana localization-01.3ds', 37);
[data38] = Nanonis_reader('Majorana localization-01.3ds', 38);
[data39] = Nanonis_reader('Majorana localization-01.3ds', 39);
[data40] = Nanonis_reader('Majorana localization-01.3ds', 40);
[data41] = Nanonis_reader('Majorana localization-01.3ds', 41);
[data42] = Nanonis_reader('Majorana localization-01.3ds', 42);
[data43] = Nanonis_reader('Majorana localization-01.3ds', 43);
[data44] = Nanonis_reader('Majorana localization-01.3ds', 44);
[data45] = Nanonis_reader('Majorana localization-01.3ds', 45);
[data46] = Nanonis_reader('Majorana localization-01.3ds', 46);
[data47] = Nanonis_reader('Majorana localization-01.3ds', 47);
[data48] = Nanonis_reader('Majorana localization-01.3ds', 48);
[data49] = Nanonis_reader('Majorana localization-01.3ds', 49);
[data50] = Nanonis_reader('Majorana localization-01.3ds', 50);
[data51] = Nanonis_reader('Majorana localization-01.3ds', 51);
[data52] = Nanonis_reader('Majorana localization-01.3ds', 52);
[data53] = Nanonis_reader('Majorana localization-01.3ds', 53);
[data54] = Nanonis_reader('Majorana localization-01.3ds', 54);
[data55] = Nanonis_reader('Majorana localization-01.3ds', 55);
[data56] = Nanonis_reader('Majorana localization-01.3ds', 56);
[data57] = Nanonis_reader('Majorana localization-01.3ds', 57);
[data58] = Nanonis_reader('Majorana localization-01.3ds', 58);
[data59] = Nanonis_reader('Majorana localization-01.3ds', 59);
[data60] = Nanonis_reader('Majorana localization-01.3ds', 60);
[data61] = Nanonis_reader('Majorana localization-01.3ds', 61);
[data62] = Nanonis_reader('Majorana localization-01.3ds', 62);
[data63] = Nanonis_reader('Majorana localization-01.3ds', 63);

%% Sorting spectroscopic data into the main data matrix
data(:, 1:3) = data1;
data(:, 4:6) = data2;
data(:, 7:9) = data3;
data(:, 10:12) = data4;
data(:, 13:15) = data5;
data(:, 16:18) = data6;
data(:, 19:21) = data7;
data(:, 22:24) = data8;
data(:, 25:27) = data9;
data(:, 28:30) = data10;
data(:, 31:33) = data11;
data(:, 34:36) = data12;
data(:, 37:39) = data13;
data(:, 40:42) = data14;
data(:, 43:45) = data15;
data(:, 46:48) = data16;
data(:, 49:51) = data17;
data(:, 52:54) = data18;
data(:, 55:57) = data19;
data(:, 58:60) = data20;
data(:, 61:63) = data21;
data(:, 64:66) = data22;
data(:, 67:69) = data23;
data(:, 70:72) = data24;
data(:, 73:75) = data25;
data(:, 76:78) = data26;
data(:, 79:81) = data27;
data(:, 82:84) = data28;
data(:, 85:87) = data29;
data(:, 88:90) = data30;
data(:, 91:93) = data31;
data(:, 94:96) = data32;
data(:, 97:99) = data33;
data(:, 100:102) = data34;
data(:, 103:105) = data35;
data(:, 106:108) = data36;
data(:, 109:111) = data37;
data(:, 112:114) = data38;
data(:, 115:117) = data39;
data(:, 118:120) = data40;
data(:, 121:123) = data41;
data(:, 124:126) = data42;
data(:, 127:129) = data43;
data(:, 130:132) = data44;
data(:, 133:135) = data45;
data(:, 136:138) = data46;
data(:, 139:141) = data47;
data(:, 142:144) = data48;
data(:, 145:147) = data49;
data(:, 148:150) = data50;
data(:, 151:153) = data51;
data(:, 154:156) = data52;
data(:, 157:159) = data53;
data(:, 160:162) = data54;
data(:, 163:165) = data55;
data(:, 166:168) = data56;
data(:, 169:171) = data57;
data(:, 172:174) = data58;
data(:, 175:177) = data59;
data(:, 178:180) = data60;
data(:, 181:183) = data61;
data(:, 184:186) = data62;
data(:, 187:189) = data63;

%% Locating the grid point where the ZBP appears in the spectroscopic data
ZBP_spectra = data(:, 1:11); % In retrospect, the 4th grid point is where the ZBP appears (confirmed)
[maxVal, max_idx] = max(ZBP_spectra, [], 2); %% Returns the max of all columns, along with the index positions of the max elements
max_data = max(maxVal); % Max of all elements in maxVal
[row, col] = find(data == max_data);
ZBP_data = data(:, col);

%% Implementing the scan range
dist = zeros(1, length(ZBP_data));
dist(1) = 0;
dist_0 = dist(1);
common_diff = 0.0267e-9; 

%% Scan range initialization
for i = 1:length(dist)
    dist(i) = dist_0 + (i - 1)*common_diff;
end

%% Computing the FWHM of the zero bias peak (to estimate the superconducting coherence length)
halfMax = (min(ZBP_data) + max(ZBP_data))/2; % Computes the position where the FWHM is located
idx_1 = find(ZBP_data >= halfMax, 1, 'first'); % Returns the position where ZBP_data first rises above halfMax
idx_2 = find(ZBP_data >= halfMax, 1, 'last'); % Returns the position where ZBP_data first falls below halfMax
FWHM = dist(idx_2) - dist(idx_1); % FWHM of ZBP_data
super_len = FWHM/2; % Rough estimate for the superconducting coherence length of the sample

%% Visualization
% Introducing a constant vertical shift between the plots
m = 62;
for i = 2:3:188
        data(:, i) = normalize(data(:, i), 'range') + m*k;
        m = m - 1;
end

% Plotting
figure(1);
for i = 2:3:188
    plot(volt_val, smooth(data(:, i)), 'b');
    xlabel('Voltage [mV]', 'FontSize', 18);
    ylabel('$\sigma \equiv dI/dV$', 'Interpreter', 'latex', 'FontSize', 18);
    title('STM tunneling spectra at various grid points', 'FontSize', 20);
    hold on
end

%% Symmetrization of the spectroscopic data exhibiting an SC gap-like feature
data_sym = zeros(128, 54);

% Symmetrization
data_trial1 = data10(:, 2);
data_trial1(65:128) = data_trial1(64:-1:1);

data_trial2 = data11(:, 2);
data_trial2(65:128) = data_trial2(64:-1:1);

data_trial3 = data12(:, 2);
data_trial3(65:128) = data_trial3(64:-1:1);

data_trial4 = data13(:, 2);
data_trial4(65:128) = data_trial4(64:-1:1);

data_trial5 = data14(:, 2);
data_trial5(65:128) = data_trial5(64:-1:1);

data_trial6 = data15(:, 2);
data_trial6(65:128) = data_trial6(64:-1:1);

data_trial7 = data16(:, 2);
data_trial7(65:128) = data_trial7(64:-1:1);

data_trial8 = data17(:, 2);
data_trial8(65:128) = data_trial8(64:-1:1);

data_trial9 = data18(:, 2);
data_trial9(65:128) = data_trial9(64:-1:1);

data_trial10 = data19(:, 2);
data_trial10(65:128) = data_trial10(64:-1:1);

data_trial11 = data20(:, 2);
data_trial11(65:128) = data_trial11(64:-1:1);

data_trial12 = data21(:, 2);
data_trial12(65:128) = data_trial12(64:-1:1);

data_trial13 = data22(:, 2);
data_trial13(65:128) = data_trial13(64:-1:1);

data_trial14 = data23(:, 2);
data_trial14(65:128) = data_trial14(64:-1:1);

data_trial15 = data24(:, 2);
data_trial15(65:128) = data_trial15(64:-1:1);

data_trial16 = data25(:, 2);
data_trial16(65:128) = data_trial16(64:-1:1);

data_trial17 = data26(:, 2);
data_trial17(65:128) = data_trial17(64:-1:1);

data_trial18 = data27(:, 2);
data_trial18(65:128) = data_trial18(64:-1:1);

data_trial19 = data28(:, 2);
data_trial19(65:128) = data_trial19(64:-1:1);

data_trial20 = data29(:, 2);
data_trial20(65:128) = data_trial20(64:-1:1);

data_trial21 = data30(:, 2);
data_trial21(65:128) = data_trial21(64:-1:1);

data_trial22 = data31(:, 2);
data_trial22(65:128) = data_trial22(64:-1:1);

data_trial23 = data32(:, 2);
data_trial23(65:128) = data_trial23(64:-1:1);

data_trial24 = data33(:, 2);
data_trial24(65:128) = data_trial24(64:-1:1);

data_trial25 = data34(:, 2);
data_trial25(65:128) = data_trial25(64:-1:1);

data_trial26 = data35(:, 2);
data_trial26(65:128) = data_trial26(64:-1:1);

data_trial27 = data36(:, 2);
data_trial26(65:128) = data_trial26(64:-1:1);

data_trial28 = data37(:, 2);
data_trial28(65:128) = data_trial28(64:-1:1);

data_trial29 = data38(:, 2);
data_trial29(65:128) = data_trial29(64:-1:1);

data_trial30 = data39(:, 2);
data_trial30(65:128) = data_trial30(64:-1:1);

data_trial31 = data40(:, 2);
data_trial31(65:128) = data_trial31(64:-1:1);

data_trial32 = data41(:, 2);
data_trial32(65:128) = data_trial32(64:-1:1);

data_trial33 = data42(:, 2);
data_trial33(65:128) = data_trial33(64:-1:1);

data_trial34 = data43(:, 2);
data_trial34(65:128) = data_trial34(64:-1:1);

data_trial35 = data44(:, 2);
data_trial35(65:128) = data_trial35(64:-1:1);

data_trial36 = data45(:, 2);
data_trial36(65:128) = data_trial36(64:-1:1);

data_trial37 = data46(:, 2);
data_trial37(65:128) = data_trial37(64:-1:1);

data_trial38 = data47(:, 2);
data_trial38(65:128) = data_trial38(64:-1:1);

data_trial39 = data48(:, 2);
data_trial39(65:128) = data_trial39(64:-1:1);

data_trial40 = data49(:, 2);
data_trial40(65:128) = data_trial40(64:-1:1);

data_trial41 = data50(:, 2);
data_trial41(65:128) = data_trial41(64:-1:1);

data_trial42 = data51(:, 2);
data_trial42(65:128) = data_trial42(64:-1:1);

data_trial43 = data52(:, 2);
data_trial43(65:128) = data_trial43(64:-1:1);

data_trial44 = data53(:, 2);
data_trial44(65:128) = data_trial44(64:-1:1);

data_trial45 = data54(:, 2);
data_trial45(65:128) = data_trial45(64:-1:1);

data_trial46 = data55(:, 2);
data_trial46(65:128) = data_trial46(64:-1:1);

data_trial47 = data56(:, 2);
data_trial47(65:128) = data_trial47(64:-1:1);

data_trial48 = data57(:, 2);
data_trial48(65:128) = data_trial48(64:-1:1);

data_trial49 = data58(:, 2);
data_trial49(65:128) = data_trial49(64:-1:1);

data_trial50 = data59(:, 2);
data_trial50(65:128) = data_trial50(64:-1:1);

data_trial51 = data60(:, 2);
data_trial51(65:128) = data_trial51(64:-1:1);

data_trial52 = data61(:, 2);
data_trial52(65:128) = data_trial52(64:-1:1);

data_trial53 = data62(:, 2);
data_trial53(65:128) = data_trial53(64:-1:1);

data_trial54 = data63(:, 2);
data_trial54(65:128) = data_trial54(64:-1:1);

%% Sorting the symmetrized data
% Starts from grid point 10
data_sym(:, 1) = data_trial1;
data_sym(:, 2) = data_trial2;
data_sym(:, 3) = data_trial3;
data_sym(:, 4) = data_trial4;
data_sym(:, 5) = data_trial5;
data_sym(:, 6) = data_trial6;
data_sym(:, 7) = data_trial7;
data_sym(:, 8) = data_trial8;
data_sym(:, 9) = data_trial9;
data_sym(:, 10) = data_trial10;
data_sym(:, 11) = data_trial11;
data_sym(:, 12) = data_trial12;
data_sym(:, 13) = data_trial13;
data_sym(:, 14) = data_trial14;
data_sym(:, 15) = data_trial15;
data_sym(:, 16) = data_trial16;
data_sym(:, 17) = data_trial17;
data_sym(:, 18) = data_trial18;
data_sym(:, 19) = data_trial19;
data_sym(:, 20) = data_trial20;
data_sym(:, 21) = data_trial21;
data_sym(:, 22) = data_trial22;
data_sym(:, 23) = data_trial23;
data_sym(:, 24) = data_trial24;
data_sym(:, 25) = data_trial25;
data_sym(:, 26) = data_trial26;
data_sym(:, 27) = data_trial27;
data_sym(:, 28) = data_trial28;
data_sym(:, 29) = data_trial29;
data_sym(:, 30) = data_trial30;
data_sym(:, 31) = data_trial31;
data_sym(:, 32) = data_trial32;
data_sym(:, 33) = data_trial33;
data_sym(:, 34) = data_trial34;
data_sym(:, 35) = data_trial35;
data_sym(:, 36) = data_trial36;
data_sym(:, 37) = data_trial37;
data_sym(:, 38) = data_trial38;
data_sym(:, 39) = data_trial39;
data_sym(:, 40) = data_trial40;
data_sym(:, 41) = data_trial41;
data_sym(:, 42) = data_trial42;
data_sym(:, 43) = data_trial43;
data_sym(:, 44) = data_trial44;
data_sym(:, 45) = data_trial45;
data_sym(:, 46) = data_trial46;
data_sym(:, 47) = data_trial47;
data_sym(:, 48) = data_trial48;
data_sym(:, 49) = data_trial49;
data_sym(:, 50) = data_trial50;
data_sym(:, 51) = data_trial51;
data_sym(:, 52) = data_trial52;
data_sym(:, 53) = data_trial53;
data_sym(:, 54) = data_trial54;

data_orig = data_sym; %% Storing the original data for future use

% Rescaling data_orig()
for j = 1:54
        data_orig(:, j) = (data_orig(:, j) - min(data_orig(:, j)))./(max(data_orig(:, j)) - min(data_orig(:, j)));
end

%% Visualization
% Introducing a constant vertical shift between the plots
m = 54;
for i = 1:54
    data_sym(:, i) = normalize(data_sym(:, i), 'range') + m*k1;
    m = m - 1;
end

% Plotting
figure(2);
for i = 1:54
    plot(volt_val, smooth(data_sym(:, i)), 'b');
    xlabel('Voltage [mV]', 'FontSize', 18);
    ylabel('$\sigma \equiv dI/dV$', 'Interpreter', 'latex', 'FontSize', 18);
    title('STM tunneling spectra at various grid points', 'FontSize', 20);
    hold on
end

%% Defining parameters
kB = (1.38/1.6)*1.0e-4; % Boltzmann's constant (in units of eV)
e = 1; % Electronic charge (in units of eV)
Ef = 0; % We take the Fermi energy level (reference point) as zero
T = 1.54; % Temperature [K]

%% Generating the BCS Dynes density of states fit
% Note: For this, we make use of the BCS Dynes density of states expression

options = optimoptions('lsqcurvefit', 'Algorithm', 'trust-region-reflective', 'OptimalityTolerance', 1e-6, 'FunctionTolerance', 1e-6); % Setting optimization options
x0 = [1e-4 5.5e-4]; % Initial guesses for the parametric fit
BCS_sol = zeros(54, 2); % Vector containing the fitting parameter values (order: [\Delta \Gamma])
gamma = zeros(54, 1); % Vector containing the thermal broadening parameter values
delta = zeros(54, 1); % Vector containing the superconducting energy gap values

BCS_fit = @(gamma, delta, tt) arrayfun(@(V) integral(@(E) real((E - 1i*gamma)./sqrt((E - 1i*gamma).^2 - delta.^2)).*(cosh((E + e*V)./2*kB*T)).^-2, Ef, Ef + e*V), tt);

for j = 1:length(gamma)
    BCS_sol(j, :) = lsqcurvefit(@(para_val, V) BCS_fit(para_val(1), para_val(2), V), x0, volt_val', data_orig(:, j)', [], [], options);
    delta(j) = BCS_sol(j, 1);
    gamma(j) = BCS_sol(j, 2);
end

BCS_data = zeros(128, 54);

for i = 1:54
    BCS_data(:, i) = BCS_fit(gamma(i), delta(i), data_orig(:, i));
end

delta_reserve = delta; %% Temporary storage of the superconducting energy gap data

%% Initializing the range of interest
d = 5.34e-11; % Spatial distance between two adjacent gridpoints
r = zeros(1, 54);
r(1) = 0.534e-9;
a = r(1);
for i = 2:length(r)
    r(i) = a + (i - 1)*d;
end

len = r; %% Temporary storage of the distance data

%% Fitting the \Delta values to an exponential model
% Rescaling delta and len
delta_reserve = (delta_reserve - min(delta_reserve))./(max(delta_reserve) - min(delta_reserve));
len = (len - min(len))./(max(len) - min(len));

% Find \Delta = 0, thereafter eliminate
delta_0 = find(delta_reserve == 0);
delta_reserve(delta_0) = [];
len(delta_0) = [];

% Locate max(delta) and initialize
idx_delta_data = find(delta_reserve == max(delta_reserve));
delta_new = delta_reserve(idx_delta_data:end);
len_new = len(idx_delta_data:end);

% Reshaping r
r(delta_0) = [];
r_new = r(idx_delta_data:end);

% Denormalizing r (for extracting the superconductinh coherence length value)
mu = mean(r_new);
standard_dev = std(r_new);

% Defining the exponential model
beta0 = [50; 0.3; 0.5]; % Initial guesses for the parametric fit
f = @(b, len_new) b(1).*exp(-len_new./b(2)) + b(3); % Exponential model

% Parametric fitting using fminsearch()
para_val = fminsearch(@(b) norm(delta_new' - f(b, len_new)), beta0); % Estimating fitting parameters

% Extract the value of \xi from the normalized fitted parameters
xi = para_val(2)*standard_dev + mu; 

%% Visualization of the fits (for normalized values of \Delta and r)
figure(3);
plot(len_new, delta_new', 'pg');
hold on
plot(len_new, f(para_val, len_new), '-r');
hold off
grid
xlabel('Distance [m] (normalized)', 'FontSize', 17);
ylabel('$\Delta$ (normalized)', 'Interpreter', 'latex', 'FontSize', 18);
title('Superconducting energy gap fitting', 'FontSize', 16);

%% BCS Dynes fit to the STM tunneling spectral data (for demonstration purposes -- does not capture all the gapped STM data)
m1 = 0;
m2 = 0;
m3 = 0;
shift = 0.4;

% First set (figure 4)
for i = 9:14
    data_orig(:, i) = data_orig(:, i) + m1*shift;
    BCS_data(:, i) = BCS_data(:, i) + m1*shift;
    m1 = m1 + 1;
end
    
% Second set (figure 5)
for i = 20:25
    data_orig(:, i) = data_orig(:, i) + m2*shift;
    BCS_data(:, i) = BCS_data(:, i) + m2*shift;
    m2 = m2 + 1;
end

% Third set (figure 6)
for i = 45:50
    data_orig(:, i) = data_orig(:, i) + m3*shift;
    BCS_data(:, i) = BCS_data(:, i) + m3*shift;
    m3 = m3 + 1;
end

% Visualization
figure(4);
for i = 9:14
    plot(volt_val, smooth(data_orig(:, i)), 'bo');
    plot(volt_val, smooth(BCS_data(:, i)) , 'r');
    xlabel('Bias voltage [mV]', 'FontSize', 18);
    ylabel('$\sigma \equiv dI/dV$', 'Interpreter', 'latex', 'FontSize', 18);
    title('BCS Dynes fit to the fully gapped STM data', 'FontSize', 15);
    hold on
end

figure(5);
for j = 20:25
    plot(volt_val, smooth(data_orig(:, j)), 'bo');
    plot(volt_val, smooth(BCS_data(:, j)) , 'r');
    xlabel('Bias voltage [mV]', 'FontSize', 18);
    ylabel('$\sigma \equiv dI/dV$', 'Interpreter', 'latex', 'FontSize', 18);
    title('BCS Dynes fit to the fully gapped STM data', 'FontSize', 15);
    hold on
end

figure(6);
for j = 45:50
    plot(volt_val, smooth(data_orig(:, j)), 'bo');
    plot(volt_val, smooth(BCS_data(:, j)) , 'r');
    xlabel('Bias voltage [mV]', 'FontSize', 18);
    ylabel('$\sigma \equiv dI/dV$', 'Interpreter', 'latex', 'FontSize', 18);
    title('BCS Dynes fit to the fully gapped STM data', 'FontSize', 15);
    hold on
end

%% Visualization of the variation of \Delta with distance

% For the actual (estimated) data
figure(7);
hold on
plot(r', delta, 'r');
title('Variation of the SC gap magnitude with distance', 'FontSize', 16);
xlabel('Distance [m]', 'FontSize', 18);
ylabel('$\Delta$ [eV]', 'Interpreter', 'latex', 'FontSize', 18);
hold off

% For the normalized values
figure(8);
hold on
plot(len_new', delta_reserve, 'r');
title('Variation of the SC gap (normalized with distance', 'FontSize', 16);
xlabel('Distance (normalized)', 'FontSize', 18);
ylabel('$\Delta$ (normalized)', 'Interpreter', 'latex', 'FontSize', 18);
hold off




















    









    