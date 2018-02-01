function [time,conc]= ODE_QS()


x = [0,0,5*10^-6,6*10^-6,0];       %molecular numbers, [AHLc,LacRs,LacRt,lambCl,AHLe]
Vext = 10;                  %Culture volume, [uL]
Vratio = 3.5*10^6;          %Vext/Vc
Vc = Vext/Vratio;           %cell volume, [uL]
c = [x(1:4)/Vc x(5)/Vext];

Tspan = 0:0.01:10;
IC = c; 
[time,conc] = ode23(@dynamics,Tspan,IC);
plot(time,conc(:,1)+conc(:,2),time,conc(:,3))

function dx = dynamics(t,x)

N = 10000;                     %number of cell
Vratio = 3.5*10^6;          %Vext/Vc
b = [0.45,0.2,0.2,0.2,0];   %synthesis rate, [uM/min]
d = [1,1,1,0,0.01];         %degradation rat, [1/min]
mu = [2, 2/Vratio];         %membrane permeabilit per volume [1/min]
%k = [500, 500*Vratio];      %factor for calculating AHL diffusing in and out of a cell
beta = [0.97,4,4];          %[uM,uM/min,uM/min]
e = [0.5,1];                %copy numbers of the sensor & toggle switch plasmid
K = [0.11,1,1];             %equilibrium constant[uM]


dx = zeros(5,1);
dx(1) = b(1)+ mu(1)*x(5) - mu(1)*x(1) - d(1)*x(1);
dx(2) = e(1)*(b(2)+beta(1)*x(1)^2/(K(1)^2+x(1)^2)) - d(2)*x(2);
dx(3) = e(2)*(b(3)+beta(2)*K(2)^3/(K(2)^3+x(5)^3)) - d(2)*x(3);
dx(4) = e(2)*(b(4)+beta(3)*K(3)^3/(K(3)^3+(x(2)+x(3))^3)) - d(3)*x(4);
dx(5) = mu(2)*N*(x(1)-x(5))-d(5)*x(5);
