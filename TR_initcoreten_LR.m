% this is for simulation data, randomly generate tensor cores Z 
function Z=TR_initcoreten_LR(S,r)
N=numel(S);
Z=cell(N,1);
z1 = cell(N,1);
for i=1:N-1
    Z{i}=randn(r(i),S(i),r(i+1));
    z1{i} = Gunfold(Z{i},2);
    [u,s,v]=svd(z1{i});   
    S1=max(s-s(15,15),0);
    Z{i} = Gfold(u*S1*v',[r(i),S(i),r(i+1)],2);
end
    Z{N}=randn(r(N),S(N),r(1));
    % below is from Qibin
%   for i=1:N
%             Z{i} = Z{i}./max(abs(Z{i}(:)));
%   end
    z1{N} = Gunfold(Z{N},2);
    [u,s,v]=svd(z1{N});
    S1=max(s-s(5,5),0);
    Z{N} =Gfold( u*S1*v',[r(N),S(N),r(1)],2);


end