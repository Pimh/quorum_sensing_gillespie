function cell=gillespieQS()
N=50;
for i=1:N
    [x,tvec]= quorumSensing();
    cell(i,:)= x(length(tvec),:); 
end
