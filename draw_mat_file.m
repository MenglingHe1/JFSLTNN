
%% disp mat file lena
index_arrangement=[2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3];
index_permute=[1 9 2 10 3 11 4 12 5 13 6 14 7 15 8 16 17];
index_repermute=[1 3 5 7 9 11 13 15 2 4 6 8 10 12 14 16 17];
index_reshape=[4 4 4 4 4 4 4 4 3];
image=imread('2_lena.bmp');
X_hat=reshape(full_tr,index_reshape);
X_hat_=permute(reshape(X_hat,index_arrangement),index_repermute);
image_hat=uint8(reshape(X_hat_,[256,256,3])*255);
figure('position', [250, 400, 900, 300]);
subplot(1,2,1);  imshow(image);title('Original Data');set(gca,'Fontsize',20);
subplot(1,2,2); imshow(image_hat);title('Decomposition Result');set(gca,'Fontsize',20);

%% disp mat file KTH
load KTH-TVT-2013.mat
X1=double(reshape(X1(:,:,:,1:180),[20,20,32,15,12]));
X_std=std(X1(:));
X=uint8(reshape(full_tr,[20,20,32,180]));
X=X.*X_std;
imshow(X(:,:,20,1))