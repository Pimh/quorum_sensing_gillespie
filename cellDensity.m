function y = cellDensity(x,t)
time=10:1:200;
j=1;
for i=1:length(t)
    if i>1 && j<= length(time)&& t(i) > time(j)
        y(j) = x(i-1);
        j=j+1;
    end
end
