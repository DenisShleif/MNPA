function [] = pa()
clc
clear
close all
warning off

comp.R1 = 1;
comp.R2 = 2;
comp.R3 = 10;
comp.R4 = 0.1;
comp.R5 = 1000;
comp.a = 100;
comp.C = 0.25;
comp.L = 0.2;

circ = initializeCircuit(comp);
DCsweep(circ);
ACSweepStable(circ);
ACSweepUnStable(comp)
end

function [o] = initializeCircuit(c)
N = 10;
o.G = zeros(N,N);
o.C = zeros(N,N);
o.F = zeros(N,1);

o.G(1,1) = 1;
o.G(2,[1,2,6]) = [-1,1,-c.R1];
o.G(3,10) = -1;
o.G(4,[6,7,8,10]) = [-1,1,1,-1];
o.G(5,[2,3]) = [-1,1];
o.G(6,[3,8]) = [1,c.R3];
o.G(7,[4,8]) = [1,c.a];
o.G(8,[4,5,9]) = [-1,1,-c.R4];
o.G(9,[5,9]) = [1,c.R5];
o.G(10,[2,7]) = [1,c.R2];

o.C(3,[1,2]) = [-c.C,c.C];
o.C(5,8) = -c.L;

o.F(1) = 1;
end

function DCsweep(circ)
N = 100;
V = linspace(-10,10,100);
Vo = zeros(1,N);
V3 = zeros(1,N);
for n = 1:N
    tmp = solve(circ,V(n));   
    Vo(n) = tmp(5);
    V3(n) = tmp(3);
end

figure;
plot(V,Vo,'.-');
hold on
plot(V,V3,'.-');
xlabel('V_{in} (V)');
ylabel('Voltage (V)');
legend({'V_o','V_3','V_oT','V_3T'});
title(sprintf('DC Voltage Sweep From %i V to %i V',min(V),max(V)));
grid on
end

function ACSweepStable(circ)
N = 100;
f = logspace(-2,1,N-1);
f = [0,f];
Vin = 5;
w = 2*pi*f;
Vo = zeros(1,N);

for n = 1:N
    tmp = solve(circ,Vin,w(n));
    Vo(n) = tmp(5);
end

figure
plot(f,abs(Vo),'.-');
xlabel('Frequency (Hz)');
ylabel('V_o (V)');
title('AC Sweep Output Voltage');
grid on

bode(f,Vo,Vin);
end

function [] = bode(f,Vo,Vin)
figure
subplot(2,1,1)
plot(f,20*log10(abs(Vo/Vin)),'.-');
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
title('Bode Plot');
grid on

subplot(2,1,2)
plot(f,rad2deg(angle(Vo/Vin)),'.-');
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
grid on
end

function ACSweepUnStable(comp)
Vin = 1;
trials = 10000;
G = zeros(1,trials);

for n = 1:trials
    comp.C = nrmrnd(0.25,0.05);
    circ = initializeCircuit(comp);
    tmp = solve(circ,Vin,pi);
    G(n) = 20*log10(abs(tmp(5)/Vin)); 
end

figure
hist(G);
xlabel('Gain (dB)');
ylabel('Count');
title('Distribution of Gain with Small Pertubations in C (\omega = \pi)');
grid on
end

function nrmmatrix = nrmrnd(mu, sigma)
nrmmatrix = mu+sigma*randn(1,1);
end


function [V] = solve(circ,Vin,w)
if nargin == 3
    V = mldivide(circ.G + w*1j*circ.C,circ.F*Vin);
elseif nargin == 2
    V = mldivide(circ.G,circ.F*Vin);
else
    V = mldivide(circ.G,circ.F);
end
end