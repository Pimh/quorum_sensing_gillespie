function delta = tauLeap()

N = 10;                     %number of cell
x = [0,0,125,125,0];       %molecular numbers, [AHLc,LacRs,LacRt,lambCl,AHLe]
x = [x;zeros(N-1,5)];
tau = .00001;                    %[sec]
t=0; 

while t<= 100
    %Update time and number of molecules 
    delta = calDelta(N,x,tau);
    x = x + delta
    t = t+tau;
    pause
end

function delta = calDelta(N,x,tau)
%Initiation
Vext = 10;               %Culture volume, [uL]
Vratio = 3.5*10^6;          %Vext/Vc
Vc = Vext/Vratio;           %cell volume, [uL]
b = [0.45,0.2,0.2,0.2,0]*500;   %synthesis rate, [uM/min]
d = [1,1,1,0,0.01];         %degradation rat, [1/min]
mu = [2, 2/Vratio];         %membrane permeabilit per volume [1/min]
k = [500, 500*Vratio];      %factor for calculating AHL diffusing in and out of a cell
beta = [0.97,4,4]*500;          %[uM,uM/min,uM/min]
e = [0.5,1];                %copy numbers of the sensor & toggle switch plasmid
K = [0.11,1,1]*500;             %equilibrium constant[uM]
delta = zeros(N,5);

%Calculate propensity

%change in AHL level within a cell
Psyn1 = poissrnd(b(1)*tau,N,1);
PdiffIn1 = poissrnd(mu(1)*k(1)*x(5)*tau/Vext,N,1);
PdiffOut1 = poissrnd(mu(1)*k(1)*x(1)*tau/Vc,N,1);
Pdeg1 = poissrnd(d(1)*x(1)*tau,N,1);
delta(:,1) = Psyn1+PdiffIn1-PdiffOut1-Pdeg1;
%change in the level of LacR expressed from the sensor plasmid
Pactivator2= poissrnd(e(1)*beta(1)*x(1)^2*tau/(K(1)^2+x(1)^2),N,1);
Psyn2 = poissrnd(e(1)*b(2)*tau,N,1);
Pdeg2 = poissrnd(x(2)*d(2)*tau,N,1);
delta(:,2) = Pactivator2+Psyn2-Pdeg2;
%change in the level of LacR expressed from the toggle switch plasmid
Prepressor3 = poissrnd(e(2)*beta(2)*K(2)^3*tau/(K(2)^3+x(4)^3),N,1);
Psyn3 = poissrnd(e(2)*b(3)*tau,N,1);
Pdeg3 = poissrnd(d(3)*x(3)*tau,N,1);
delta(:,3) = Prepressor3+Psyn3-Pdeg3;
%change in lambdaCI level
Prepressor4 = poissrnd(e(2)*beta(3)*K(3)^3*tau/(K(3)^3+(x(2)+x(3))^3),N,1);
Psyn4 = poissrnd(e(2)*b(4)*tau,N,1);
Pdeg4 = poissrnd(x(4)*d(3)*tau,N,1);
delta(:,4) = Prepressor4+Psyn4-Pdeg4;
%change in external AHL level 
PdiffIn5 = sum(poissrnd(mu(2)*k(1)*x(1)*tau/Vc,N,1));
PdiffOut5 = sum(poissrnd(mu(2)*k(2)*x(5)*tau/Vext,N,1));
Pdeg5 = poissrnd(x(5)*d(5)*tau);
delta(1,5) = PdiffIn5-PdiffOut5-Pdeg5; 

