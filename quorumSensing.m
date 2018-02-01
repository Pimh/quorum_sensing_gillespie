function [time,conc,Nswitch]=quorumSensing()


x = [0,0,100,1000,0];       %molecular numbers, [AHLc,LacRs,LacRt,lambCl,AHLe]


b = [0.2,.2,.54,.5,.2]*500;
beta = [.97,4,4]*500;
K = [0.11,1,1]*500;  
Tspan = 0:0.01:1000;
e = [0.5,1,5];                %copy numbers of the sensor & toggle switch plasmid

%b1 = [.679,.776,.873,.97,1.067,1.164,1.261]*500;
b1=[.5,1,1.5,2,2.5,3,3.5];
colr=['.r','.g','.b'];
IC = x;
for i=1:7   
    e(2)= b1(i)*e(1);
    [time,conc] = ode23(@(t,x)dynamics(t,x,b,beta,K,e),Tspan,IC);
    conc=conc/500;
    tswitch = length(time)-sum(conc(:,2)+conc(:,3)>conc(:,4));
    N = 3*10^6*2.^(.0024*(Tspan));
    if tswitch >0
        Nswitch(i) = N(tswitch)/10
    else
        Nswitch(i) = 0
    end
    %figure(1)
%     plot(time,conc(:,3)+conc(:,2),'.b',time,conc(:,4),'.m',time,conc(:,1),'.r',time,conc(:,5),'.y')
%     figure(2)
    %plot(b1(i),Nswitch(i),'*m')
    %hold on
    %pause
end


function dx = dynamics(t,x,b,beta,K,e)

N0 = 3*10^6;                     %number of cell
d = [1,1,1,0,1];         %degradation rat, [1/min]

Vext = 10;                  %Culture volume, [uL]
Vratio = 3.5*10^6;          %Vext/Vc
Vc = Vext/Vratio;           %cell volume, [uL]
mu=2;AHLc=100;
N = N0*2^(.0024*(t));
AHLe= N*AHLc*Vc/Vext;

dx = zeros(5,1);
%AHL
dx(1) = b(4)- d(1)*x(1)+ mu*AHLe; %+ (beta(1)*x(1)^2/(K(1)^2+x(1)^2))) - d(1)*x(1);
%LacR sensor
dx(2) = e(1)*(b(1)+beta(1)*x(1)^2/(K(1)^2+x(1)^2)) - d(2)*x(2);
%LacR toggle
dx(3) = e(2)*(b(2)+beta(2)*K(2)^3/(K(2)^3+x(4)^3)) - d(2)*x(3);
%lambda CI
dx(4) = e(2)*(b(3)+beta(3)*K(3)^3/(K(3)^3+(x(2)+x(3))^3)) - d(3)*x(4);
dx(5) = e(3)*(b(5) + (beta(2)*K(2)^3/(K(2)^3+x(4)^3)))-d(5)*x(5);
