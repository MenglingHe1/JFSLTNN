%% 9_D reshape 4,4,4,4,4,4,4,4,3
image=imread('2_lena.bmp');
image=imresize(image,[256 256]);
X=double(image);
% using permute and reshape to make the 1st order of the tensor contains a
% unit information, then the following orders are the expandion
index_arrangement=[2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3];
index_permute=[1 9 2 10 3 11 4 12 5 13 6 14 7 15 8 16 17];
index_repermute=[1 3 5 7 9 11 13 15 2 4 6 8 10 12 14 16 17];
index_reshape=[4 4 4 4 4 4 4 4 3];
X=permute(reshape(X,index_arrangement),index_permute);
X=reshape(X,index_reshape)/255; % load data and reshape
S=size(X);
%%  5_D reshape 
image=imread('facade.png');
image=imresize(image,[256 256]);
X=double(image);
index_arrangement=[4,4,4,4,4,4,4,4,3];
index_permute=[1 5 2 6 3 7 4 8 9];
index_repermute=[1 3 5 7 2 4 6 8 9];
index_reshape=[16 16 16 16 3];
X=permute(reshape(X,index_arrangement),index_permute);
X=reshape(X,index_reshape)/255; % load data and reshape
S=size(X);
%% 3_D original 
image=imread('facade.png');
image=imresize(image,[256 256]);
X=double(image)/255;
S=size(X);
%% Hyperspectral Image 960x1200x33
load HSI.mat
X=reshape(X,[32,30,40,30,33])/std(X_(:));
S=size(X);
%% Hyperspectral Image 256x256x33
load HSI.mat
image=X__;
index_arrangement=[4,4,4,4,4,4,4,4,33];
index_permute=[1 5 2 6 3 7 4 8 9];
index_repermute=[1 3 5 7 2 4 6 8 9];
index_reshape=[16 16 16 16 33];
X=permute(reshape(image,index_arrangement),index_permute);
X=reshape(X,index_reshape)/std(X(:)); % load data and reshape
S=size(X);
%% Hyperspectral Image 64x64x33
load HSI.mat
image=X___;
index_arrangement=[2 2 2 2 2 2 2 2 2 2 2 2 33];
index_permute=[1 7 2 8 3 9 4 10 5 11 6 12 13];
index_repermute=[1 3 5 7 9 11 2 4 6 8 10 12 13];
index_reshape=[4 4 4 4 4 4 33];
X=permute(reshape(image,index_arrangement),index_permute);
X=reshape(X,index_reshape)/std(X(:)); % load data and reshape
S=size(X);
%%
tic
%tr = tensor_ring(X,'Tol',1e-2,'Alg','ALS','Rank',r)
%tr = tensor_ring(X,'Tol',1e-4,'Alg','ALSAR')
r=36*ones(1,9);
%  S=size(X);
%  Z=TR_initcoreten_abs(S,r);

    % ,'Initcore',Z
%tr = tensor_ring(X,'Tol',1e-2,'Alg','SGD','Rank',r)
tr = tensor_ring(X,'Tol',1e-10,'Alg','SGD','Rank',r)
% r=36*ones(1,9);
% tr =tensor_ring(X,'Tol',1e-1,'Alg','SGD','Rank',r)
t=toc;
X_hat=reshape(full(tr),S);
err = X_hat(:) - X(:);
RSE = sqrt(sum(err.^2)/sum(X(:).^2))
disp(['Optimization time is  ',num2str(t)]);
%% Evaluation
X_hat=reshape(full(tr),S);
X_hat_=permute(reshape(X_hat,index_arrangement),index_repermute);
image_hat=uint8(reshape(X_hat_,[256,256,3])*255);
figure('position', [250, 400, 900, 300]);
subplot(1,2,1);  imshow(image);title('Original Data');set(gca,'Fontsize',20);
subplot(1,2,2); imshow(image_hat);title('Decomposition Result');set(gca,'Fontsize',20);
%
err = X_hat(:) - X(:);
error_X=X-mean(X(:)).*ones(S);
RSE = sqrt(sum(err.^2)/sum(X(:).^2));
RRSE = sqrt(sum(err.^2)/sum((error_X(:)).^2));
RMSE=sqrt((sum(err).^2)/numel(X));
PSNR = PSNR_RGB(double(image_hat),double(image));
SSIM = ssim_index(rgb2gray(uint8(image_hat)),rgb2gray(uint8(image)));
fprintf('Tensor size');disp(S);
fprintf('TR rank');disp(r);
fprintf('Data evaluation: RSE = %g, RRSE = %g, RMSE = %g\n',RSE,RRSE,RMSE);
fprintf('Image evaluation: PSNR = %g, SSIM = %g\n',PSNR,SSIM);
disp(['Optimization time is  ',num2str(t)]);
%% add noise 
% core
N=numel(S);
Z=tr.node;
Z_=tr.node;
for core_num=1:N-1
Z_{core_num} = Z{core_num} + 0.1*randn(size(Z_{core_num}));
X_noise=permute(reshape(coreten2tr(Z_),index_arrangement),index_repermute);
image_noise=uint8(reshape(X_noise,[256,256,3])*255);
%imshow(image_noise,'InitialMagnification', 800)
subtightplot(1,N-1,core_num);  imshow(image_noise,'InitialMagnification', 1800);
end
%% Hyperspectral Image Evaluation
X_hat=reshape(full(tr),S);
err = X_hat(:) - X(:);
error_X=X-mean(X(:)).*ones(S);
RSE = sqrt(sum(err.^2)/sum(X(:).^2));
RRSE = sqrt(sum(err.^2)/sum((error_X(:)).^2));
RMSE=sqrt((sum(err).^2)/numel(X));
fprintf('Tensor size');disp(S);
%fprintf('TR rank');disp(r);
fprintf('Data evaluation: RSE = %g, RRSE = %g, RMSE = %g\n',RSE,RRSE,RMSE);
%disp(['Optimization time is  ',num2str(t)]);
X_hat_=permute(reshape(X_hat,index_arrangement),index_repermute)*std(X(:));
X_=reshape(X_hat_,[256,256,33]);
dispRGB(X_,600)
%% Display Sampling results of HSI
load HSI_8image_data.mat
subtightplot(1,5,1);  dispRGB(it4e5,100);title('4e5, RSE=1.09');set(gca,'Fontsize',10); xlabel('20%') 
subtightplot(1,5,2);  dispRGB(it8e5,100);title('8e5, RSE=0.41');set(gca,'Fontsize',10); xlabel('40%') 
subtightplot(1,5,3);  dispRGB(it12e5,100);title('1.2e6, RSE=0.27');set(gca,'Fontsize',10); xlabel('60%') 
subtightplot(1,5,4);  dispRGB(it16e5,100);title('1.6e6, RSE= 0.20');set(gca,'Fontsize',10); xlabel('80%') 
% subtightplot(1,8,5);  dispRGB(it2e6,100);title('2e6, RSE=0.23');set(gca,'Fontsize',10);
% subtightplot(1,8,6);  dispRGB(it4e6,100);title('4e6, RSE=0.20');set(gca,'Fontsize',10);
%subtightplot(1,5,5);  dispRGB(it10e6,100);title('1e7, RSE=0.18');set(gca,'Fontsize',10); xlabel('500%') 
%%
 load coil100_128x128.mat
 load coil100_32x32.mat
X=reshape(imgs,[32,32,3,72,100]);
X=X(:,:,:,1:12,1:100);
% X=uint8(X(:,:,:,1:200));
% imdisp(X)
X=reshape(X,[32,32,3,12,4,25]);
X_std=std(X(:));
X=X/X_std;
r=24*ones(1,6);
S=size(X);
% imshow(uint8(imgs(:,:,:,73)),'InitialMagnification', 400)
tic
%tr = tensor_ring(X,'Tol',1e-2,'Alg','ALS','Rank',r)
%tr = tensor_ring(X,'Tol',1e-2,'Alg','ALSAR')
tr = tensor_ring(X,'Tol',1e-2,'Alg','SGD','Rank',r)
t=toc;
X_hat=reshape(full_tr,S);
err = X_hat(:) - X(:);
RSE = sqrt(sum(err.^2)/sum(X(:).^2))
X_=reshape(X_hat,[32,32,3,288]).*X_std;
imdisp(uint8(X_))
disp(['Optimization time is  ',num2str(t)]);
%% KTH
 load KTH-TVT-2013.mat
X=double(reshape(X1(:,:,:,1:180),[20,20,32,15,12]));
X_std=std(X(:));
X=X/X_std;
r=30*ones(1,10);
S=size(X);
tic
%tr = tensor_ring(X,'Tol',1e-2,'Alg','ALS','Rank',r)
%tr = tensor_ring(X,'Tol',1e-2,'Alg','ALSAR')
tr = tensor_ring(X,'Tol',1e-2,'Alg','SGD','Rank',r)
t=toc;
X_hat=reshape(full(tr),S);
err = X_hat(:) - X(:);
RSE = sqrt(sum(err.^2)/sum(X(:).^2))
X_=reshape(X_hat,[20,20,32,15,12]).*X_std;
disp(['Optimization time is  ',num2str(t)]);