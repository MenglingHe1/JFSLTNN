%%
% load coil100_128x128.mat
 load cifar_10_all.mat
X=double(testdata);
X=reshape(X,[32,32,3,10000]);
X_data=X(:,:,:,1:500);
% X=uint8(X(:,:,:,1:200));
% imdisp(X)
X=reshape(X_data,[32,32,3,500]);
% X=reshape(X,[4,4,2,4,4,2,3,500]);
% X=permute(X,[1,4,2,5,3,6,7,8]);
% X=reshape(X,[16,16,4,3,500]);
X_std=std(X(:));
X=X./X_std;
tr = tensor_ring(X,'Tol',1e-1,'Alg','SGD','Rank',[5     3    50    62]);
node=tr.node;
Np=0;
for i=1:numel(node)
    Np=Np+numel(node{i});
end
cr=Np/numel(X);
mean_r=mean(tr.r);
S=size(X);
r=tr.r';
fprintf('----evaluation-------\n');
disp(S); disp(r);
X_hat=reshape(full(tr),S);
err=X_hat(:)-X(:);
RSE=sqrt(sum(err.^2)/sum(X(:).^2));
fprintf('Np=%d, mean=%g, compress rate=%f, RSE=%f\n',Np,mean_r,cr,RSE);
