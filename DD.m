function D = DD(n)
diaga=ones(n,1);diagb=ones(n-1,1);
D=diag(diaga)+diag(-diagb,1);
D(end,1)=-1;
end

