
clc
clear  

addpath( 'data','quality_assess', 'TR_function', 'Toolbox')
addpath('D:\程序\JFSLTNN\TR_function\support_tools')
addpath('D:\程序\JFSLTNN\TR_function\corefun')
addpath('D:\程序\JFSLTNN\Toolbox\tensorlab')
addpath('D:\程序\JFSLTNN\Toolbox\tensor_toolbox')
addpath('D:\程序\JFSLTNN\Toolbox\poblano_toolbox_1.1')
load('lemon slices.mat')
 S=HSI; 
S=double(S);
S=S/max(S(:));
sf = 8;
s0= sf/2;

  [M, N, L]=size(S); 
  S1=hyperConvert2D(S);
F = [2  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  2  6 11 17 21 22 21 20 20 19 19 18 18 17 17;...
     1  1  1  1  1  1  2  4  6  8 11 16 19 21 20 18 16 14 11  7  5  3  2  2  1  1  2  2  2  2  2;...
     7 10 15 19 25 29 30 29 27 22 16  9  2  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1];
  for band = 1:3
        div = sum(F(band,:));
        for i = 1:31
            F(band,i) = F(band,i)/div;
        end
  end

 psf=fspecial('gaussian',[7,1],2);
 % psf=ones(8,1)/8;
  
 BW1=psf2otf(psf,[M 1]);
 S_w=ifft(fft(S).*repmat(BW1,1,N,L)); 
 BH1=psf2otf(psf,[N 1]);
aa=fft(permute(S_w,[2 1 3]));
  S_h=(aa.*repmat(BH1,1,M,L));
 S_h= permute(ifft(S_h),[2 1 3]);  
  Y_h=S_h(s0:sf:end,s0:sf:end,:);
  Y_h_bar=hyperConvert2D(Y_h);
  
  SNRh=  30 ;
sigmam = sqrt(sum(Y_h_bar(:).^2)/(10^(SNRh/10))/numel(Y_h_bar));
rng(10,'twister')
   Y_h_bar = Y_h_bar+ sigmam*randn(size(Y_h_bar));
HSI=hyperConvert3D(Y_h_bar,M/sf, N/sf );

 %%  simulate HR-MSI
 rng(10,'twister')
Y = F*S1;
SNRm=  35  ;
sigmam = sqrt(sum(Y(:).^2)/(10^(SNRm/10))/numel(Y));
Y = Y+ sigmam*randn(size(Y));
MSI=hyperConvert3D(Y,M,N);

%% FSTRD_TV_PAM
%%%%%%%%%%%%%%%%%%%%%%   参数设置     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r    =   5 ;
par.TRrank = [r,300,r];
%%
par.lambda =  1.1; 
par.mu     =   0.009;
par.eta   =  1.5; 
par.rho    =   6; 

%%%%%%%%%%%%%%%%%%%%%%   模型算法   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Z0,runtime] = WSLTNN( HSI,MSI,F,BW1,BH1,S,par ,s0,sf);
[psnr4,rmse4,ergas4,sam4,uiqi4,ssim4,DD4,CC4] = quality_assessment(Z0*255, S*255, 0, 1.0/sf);
fprintf('================== Result =====================\n');
fprintf(' %8.8s  %5.4s  %5.4s  %6.5s  %4.3s  %5.4s  \n','method','PSNR', 'SSIM', 'ERGAS', 'SAM', 'UIQI' );
fprintf(' %8.8s  %2.3f  %2.3f  %2.3f  %2.3f  %2.3f  \n','Ours', psnr4, ssim4,ergas4,sam4,uiqi4);


% 
% % band_psnr
% e1=[]; 
%  a=double(im2uint8(S));
%  b=double(im2uint8(Z0)); 
%  for i=1:31
%          e=a-b; 
%           e=e(1:256,1:256,i); 
%           mse = mean(mean(e.^2));
%          PSNR= 10*log10(255^2/mse); 
%           e1(i)=PSNR;
%  end
%    save('LOGLRTR_band_psnr_paviau','e1');
%   x=1:1:31;        %x轴上的数据
%  load('LOGLRTR_band_psnr_paviau.mat');%读取mat文件
%  plot(x,e1,'-r');         %线性，颜色，标记
%  axis([0,35,30,50]);       %确定x轴和y轴框图的大小
%  set(gca,'XTick',[0:5:35]); 
%  set(gca,'YTick',[30:5:50]);
%  xlabel('Band');
%  ylabel('PSNR');
% 
% % %%band_uiqi
% e2=[];
% a=double(im2uint8(S));
% b=double(im2uint8(Z0));
% q_band = zeros(1,93);
% sum=0;
%  for i=1:31
%     q_band(i)=img_qi(a(:,:,i),b(:,:,i), 32);
%     c=img_qi(a(:,:,i),b(:,:,i), 32);
%  e2(i)=c;
%  end
%  x=1:1:31;        %x轴上的数据
%  plot(x,e2,'-r');         %线性，颜色，标记
%  axis([0,35,0.56,0.96]);        %确定x轴和y轴框图的大小
%  set(gca,'XTick',[0:5:35]);
%  set(gca,'YTick',[0.56:0.05:0.96]);
%  xlabel('Band');
%  ylabel('UIQI');
%   save('LOGLRTR_band_uiqi_paviau','e2');


