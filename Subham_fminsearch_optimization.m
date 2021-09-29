function [ds1] = ds1(C)
%{
the variable in the square brakets "ds" will be what the variable the
function outputs. The function is called using the file name "ds" and
inputting a variable that has at least 4 elements.
%}


%{
Loads data.mat so that we can use the x values in the definition of f
%}
load('2.75.mat');


%%%====================================== 
T = 2; % Temperature at which the experimental dataset was obtained (in K)
kb = 1.3807*1.0e-4; % Boltzmann's constant (in units of eV)
e = 1; 
deltamV= C(1);
taomV= C(2);            
valmV= C(3);
%sqrt2=sqrt(2);

delta=deltamV*0.001;
tao=taomV*0.001;
val=valmV*0.001;
%Vmod=VmodmV*0.001;


cond = abs(real(complex(x,tao)./sqrt(complex(x,tao).^2-delta^2))).*(1+cosh((e*x)./2*kb*T)).^-1;

tmp=x-val;
[idx idx]= min(abs(tmp));                       %Normalization
norm = y(idx);
cond=cond/cond(idx)*norm;


%f = @(V) BCS_fit(C(1), C(2), V);

%f = @(gamma, delta, tt) arrayfun(@(V) integral(@(E) real(abs(E - 1i*gamma)./sqrt((E - 1i*gamma).^2 - delta.^2)).*(cosh((E + e.*V)./2*k*T)).^-2, Ef, Ef + e*V), tt);
%{
The general form of a sine function. Be sure that all parentheses are
where they need to be. Alos, notice that the 1,2,3, and 4th constants are
referred to as C(1), C(2), C(3), and C(4).
%}

ds1 = (cond-y).^2;
%{
Defines an array where each element is equal to the difference squared
between corresponding elements of f and y. Notice the use of ".^2" instead
of "^2". The period tells Matlab to perform the calculation element by
element instead of through matrix multiplication.
%}

ds1 = sum(ds1);
%{
Adds all elements of ds together. This last line will be the variable that
the function "ds" returns, and is the function to be minimized by
fminsearch.
%}