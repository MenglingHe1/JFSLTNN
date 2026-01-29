
clc
clear  

addpath( 'data','quality_assess', 'TR_function', 'Toolbox')
addpath('D:\程序\JFSLTNN\TR_function\support_tools')
addpath('D:\程序\JFSLTNN\TR_function\corefun')
addpath('D:\程序\JFSLTNN\Toolbox\tensorlab')
addpath('D:\程序\JFSLTNN\Toolbox\tensor_toolbox')
addpath('D:\程序\JFSLTNN\Toolbox\poblano_toolbox_1.1')
 load('Pavia.mat')
 %S=HSI; 
S=double(S);
S=S/max(S(:));
sf = 4;
s0= sf/2;

  [M, N, L]=size(S); 
  S1=hyperConvert2D(S);
% F = [2  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  2  6 11 17 21 22 21 20 20 19 19 18 18 17 17;...
%      1  1  1  1  1  1  2  4  6  8 11 16 19 21 20 18 16 14 11  7  5  3  2  2  1  1  2  2  2  2  2;...
%      7 10 15 19 25 29 30 29 27 22 16  9  2  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1];
%   for band = 1:3
%         div = sum(F(band,:));
%         for i = 1:31
%             F(band,i) = F(band,i)/div;
%         end
%   end

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
r    =  6  ;
par.TRrank = [r,300,r];
%%
par.lambda =   0.006   ; 
par.mu     =   0.002        ;
par.eta   =   2       ; 
par.rho    =  31      ; 

%%%%%%%%%%%%%%%%%%%%%%   模型算法   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Z0,runtime] = WSLTNN( HSI,MSI,F,BW1,BH1,S,par ,s0,sf);
[psnr4,rmse4,ergas4,sam4,uiqi4,ssim4,DD4,CC4] = quality_assessment(Z0*255, S*255, 0, 1.0/sf);
fprintf('================== Result =====================\n');
fprintf(' %8.8s  %5.4s  %5.4s  %6.5s  %4.3s  %5.4s  \n','method','PSNR', 'SSIM', 'ERGAS', 'SAM', 'UIQI' );
fprintf(' %8.8s  %2.3f  %2.3f  %2.3f  %2.3f  %2.3f  \n','Ours', psnr4, ssim4,ergas4,sam4,uiqi4);


% 
% % % band_psnr
% % e1=[]; 
% %  a=double(im2uint8(S));
% %  b=double(im2uint8(Z0)); 
% %  for i=1:93
% %          e=a-b; 
% %           e=e(1:256,1:256,i); 
% %           mse = mean(mean(e.^2));
% %          PSNR= 10*log10(255^2/mse); 
% %           e1(i)=PSNR;
% %  end
% %    save('LOGLRTR_band_psnr_paviau','e1');
% %   x=1:1:93;        %x轴上的数据
% %  load('LOGLRTR_band_psnr_paviau.mat');%读取mat文件
% %  plot(x,e1,'-r');         %线性，颜色，标记
% %  axis([0,35,30,50]);       %确定x轴和y轴框图的大小
% %  set(gca,'XTick',[0:5:35]); 
% %  set(gca,'YTick',[30:5:50]);
% %  xlabel('Band');
% %  ylabel('PSNR');
% % 
% % % %%band_uiqi
% % e2=[];
% % a=double(im2uint8(S));
% % b=double(im2uint8(Z0));
% % q_band = zeros(1,93);
% % sum=0;
% %  for i=1:93
% %     q_band(i)=img_qi(a(:,:,i),b(:,:,i), 32);
% %     c=img_qi(a(:,:,i),b(:,:,i), 32);
% %  e2(i)=c;
% %  end
% %  x=1:1:93;        %x轴上的数据
% %  plot(x,e2,'-r');         %线性，颜色，标记
% %  axis([0,35,0.56,0.96]);        %确定x轴和y轴框图的大小
% %  set(gca,'XTick',[0:5:35]);
% %  set(gca,'YTick',[0.56:0.05:0.96]);
% %  xlabel('Band');
% %  ylabel('UIQI');
% %   save('LOGLRTR_band_uiqi_paviau','e2');
% 
% R = Z0(:, :, 3);
% G = Z0(:, :, 31);
% B = Z0(:, :, 59);
% 
% % 进行归一化处理，将像素值调整到[0, 1]区间
% R_normalized = (R - min(R(:))) / (max(R(:)) - min(R(:)));
% G_normalized = (G - min(G(:))) / (max(G(:)) - min(G(:)));
% B_normalized = (B - min(B(:))) / (max(B(:)) - min(B(:)));
% % 合成伪彩色图像
% pseudo_color_image = cat(3, R_normalized, G_normalized, B_normalized);
% % 显示伪彩色图像
% % figure;
% % 指定要放大的区域
% x_start = 190;  % 区域的左上角x坐标
% y_start = 55;  % 区域的左上角y坐标
% width = 35;    % 放大区域的宽度
% height =35;   % 放大区域的高度
% % 提取该区域
% zoomed_region = pseudo_color_image(y_start:(y_start+height-1), x_start:(x_start+width-1), :);
% % 创建一个新的figure
% figure;
% imshow(pseudo_color_image,'Border','tight');
% % imshow(pseudo_color_image);
% % axis tight; % 去除轴周围的空白区域
% % 在伪彩色图像上画出放大区域的矩形框
% hold on;
% rectangle('Position', [x_start, y_start, width, height], 'EdgeColor', 'r', 'LineWidth', 1);
% % 创建一个嵌套的axes用于显示放大后的区域
% zoom_factor = 4; % 放大倍数
% zoomed_region_resized = imresize(zoomed_region, zoom_factor, 'nearest');
% % 放大区域的位置和大小
% zoom_width = size(zoomed_region_resized, 1);
% zoom_height = size(zoomed_region_resized, 1);
% x_zoom = x_start + width + 222;  % 放大区域的位置, 相对于原图的偏移
% y_zoom = y_start+60;  % 调整这个值来设置放大图的位置
% % 确保放大区域不会超出图像边界
% if x_zoom + zoom_width > size(pseudo_color_image, 2)
%     x_zoom = size(pseudo_color_image, 2) - zoom_width;
% end
% if y_zoom + zoom_height > size(pseudo_color_image, 1)
%     y_zoom = size(pseudo_color_image, 1) - zoom_height;
% end
% % 显示放大后的区域
% axes('Position', [(x_zoom/size(pseudo_color_image, 2)), (1-y_zoom/size(pseudo_color_image, 1)-zoom_height/size(pseudo_color_image, 1)), (zoom_width/size(pseudo_color_image, 2)), (zoom_height/size(pseudo_color_image, 1))]);
% imshow(zoomed_region_resized);
% hold on;
% line([0, zoom_width], [1, 1], 'Color', 'w', 'LineWidth', 1); % 上沿
% line([1, 1], [0, zoom_height], 'Color', 'w', 'LineWidth', 1); % 左沿
% hold off;
% frame = getframe(gcf);
% imwrite(frame.cdata, 'LOGLRTR_color_PaviaU.png');
% 
% 
% a = S;
% b = Z0;
% c = a - b;
% figure
% imagesc(c(:, :, 59));  
% axis image;
% caxis([0 0.05]);
% colormap(jet);
% set(gca, 'XTickLabel', [], 'YTickLabel', []);
% axis off; 
% 
% % figure;
% % caxis([0 0.1]);
% % colormap(jet);
% % axis off;
% % 
% % % 创建色柱并获取句柄
% % c = colorbar('southoutside');
% % 
% % % 调整色柱的大小和位置
% % % 'Position' 的值是 [left bottom width height]，根据需要调整
% % set(c, 'Position', [0.1 0.4 0.8 0.04]);
% %  saveas(gcf, 'colorbar-Yellow peppers.png');
