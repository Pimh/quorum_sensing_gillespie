function [x,tvec]=gillespieQS2()

%%%%%%%%%%%%%%%%%
%initiation

%molecules involved
LacR_s=0;
LacR_t=20;
CI=1000;
AHL=0;
GFP=0;

%constants
b = [1,1,1,0.1,1];   %basal synthesis rate, [uM/min]
d = [1,1,1,1,1];         %degradation rat, [1/min]
beta = [2,20,20];          %[uM,uM/min,uM/min]
e = [5,50,100];                %copy numbers of the sensor & toggle switch plasmid
K = [.5,5,5];             %equilibrium constant[uM]
N=1;tau=0;
%time
t=0; tend=100;
tvec=[];x=[];

while t <= tend
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate propensity
    
    %rxn1: Synthesis of LacR_s
    a(1)= e(1)*(b(1) + (beta(1)*AHL^2/(K(1)^2+AHL^2)));
    %rxn2: Degradation of LacR_s
    a(2)= LacR_s*d(1);
    %rxn3: Synthesis of LacR_t
    a(3)= e(2)*(b(2) + (beta(2)*K(2)^3/(K(2)^3+CI^3)));
    %rxn4: Degradation of LacR_t
    a(4)= LacR_t*d(2);
    %rxn5: Synthesis of CI
    a(5)= e(2)*(b(3) + (beta(3)*K(3)^3/(K(3)^3+(LacR_s+LacR_t)^3)));
    %rxn6: Degradation of CI
    a(6)= CI*d(3);
    %rxn7: Synthesis of AHL
    a(7)= e(1)*(b(4) + (beta(1)*AHL^2/(K(1)^2+AHL^2)));
    %rxn8: Degradation of AHL
    a(8)= AHL*d(4);
    %rxn9: Synthesis of GFP
    a(9)= e(3)*(b(5) + (beta(2)*K(2)^3/(K(2)^3+CI^3)));
    %rxn10: Degradation of GFP
    a(10)= GFP*d(5);
    
    atot=sum(a);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %Choose a reaction and update molecule number
    r1 = rand(1);
    
    if r1<= a(1)/atot
        LacR_s= LacR_s+1;
    elseif r1<= sum(a(1:2))/atot
        LacR_s= LacR_s-1;
    elseif r1<= sum(a(1:3))/atot
        LacR_t= LacR_t+1;
    elseif r1<= sum(a(1:4))/atot
        LacR_t= LacR_t-1;
    elseif r1<= sum(a(1:5))/atot
        CI = CI+1;
    elseif r1<= sum(a(1:6))/atot
        CI = CI-1;
    elseif r1<= sum(a(1:7))/atot
        AHL = AHL+1;
    elseif r1<= sum(a(1:8))/atot
        AHL = AHL-1;
    elseif r1<= sum(a(1:9))/atot
        GFP = GFP+1;
    elseif r1<= sum(a(1:10))/atot
        GFP = GFP-1;
    end
    
    AHL = AHL+ N*exp(.05*tau);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate tau and update time
    tau = -log(rand)/atot;
    t=t+tau;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Record molecular numbers at a given timepoint 
    tvec = [tvec t];
    x= [x;LacR_s,LacR_t,CI,AHL,GFP];
end