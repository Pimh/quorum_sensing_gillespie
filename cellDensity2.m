function y = cellDensity2()

for i=1:100
    [x,tvec,rac]=gillespieQS3();
    y(i) = x(length(x(:,5)),5);
end