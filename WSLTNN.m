function [imgBest,runtime] = WSLTNN(HSI,MSI,F,BW,BH,S,par ,s0,sf)
Rank = par.TRrank;psnr_max=0;
%%  simulate LR-HSI
HSI1 = tenmat_sb(HSI,1);HSI2 = tenmat_sb(HSI,2);HSI3 = tenmat_sb(HSI,3);
%%  simulate HR-MSI
MSI1 = tenmat_sb(MSI,1);MSI2 = tenmat_sb(MSI,2);MSI3 = tenmat_sb(MSI,3);
%% inilization D1 D2 D3
    core = tensor_ring1(MSI,'Tol',1e-3,'Alg','ALS','Rank',Rank); %1e-2  
G = core.node;

D1 = G{1}; D2 = G{2};
D_1 = ifft(fft(Gunfold(D1,2)).*repmat(BW,[1 Rank(1)*Rank(2)]));
D_1 = Gfold(D_1(s0:sf:end,:),[Rank(1),size(HSI,1),Rank(2)],2);

D_2 = ifft(fft(Gunfold(D2,2)).*repmat(BH,[1 Rank(2)*Rank(3)]));
D_2 = Gfold(D_2(s0:sf:end,:),[Rank(2),size(HSI,2),Rank(3)],2);

D3 = vca(HSI3,Rank(3)*Rank(1));
D_3 = Gfold(F*D3,[Rank(3),size(MSI,3),Rank(1)],2);
D3 =  Gfold(D3,[Rank(3),size(HSI,3),Rank(1)],2); 

D11{1} = D_1;D11{2} = D_2;D11{3} = D3;
D22{1} = D1;D22{2} = D2;D22{3} = D_3;
%% iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = zeros(size(S));
imax = 30;i = 0;tic 
while i<imax
    i=i+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  update D1
P1 = tenmat_sb(Z_neq(D11',1),2);  P1=P1';
P2 = tenmat_sb(Z_neq(D22',1),2);  P2=P2';
DG=DD(size(D1,2));


[ D1_M ] = WST1(  D1,HSI1, MSI1,P1,P2,BW,s0,sf,par,DG);

D1 = Gunfold(D1_M,2);
D_1 = ifft(fft(D1).*repmat(BW,[1 Rank(1)*Rank(2)]));
D1 =  Gfold(D1,[Rank(1),size(MSI,1),Rank(2)],2); 
D_1 = Gfold(D_1(s0:sf:end,:),[Rank(1),size(HSI,1),Rank(2)],2);

D11{1} = D_1;D22{1} = D1;
%%  update D2
P1 = tenmat_sb(Z_neq(D11',2),2);  P1 = P1';
P2 = tenmat_sb(Z_neq(D22',2),2);  P2 = P2';
DG=DD(size(D2,2));

[ D2_M ] = WST1(D2,HSI2, MSI2,P1,P2,BH,s0,sf,par,DG);

D2 = Gunfold(D2_M,2);
D_2 = ifft(fft(D2).*repmat(BH,[1 Rank(2)*Rank(3)]));
D2 =  Gfold(D2,[Rank(2),size(MSI,2),Rank(3)],2); 
D_2 = Gfold(D_2(s0:sf:end,:),[Rank(2),size(HSI,2),Rank(3)],2);

D11{2} = D_2;D22{2} = D2;
%%  update D3
P1 = tenmat_sb(Z_neq(D11',3),2);  P1=P1';
P2 = tenmat_sb(Z_neq(D22',3),2);  P2=P2';
DG=DD(size(D3,2));

[ D3_M ] =WST3(D3, HSI3,MSI3,P1,P2,F,par,DG );

D3 = Gunfold(D3_M,2);
D_3 = Gfold(F*D3,[Rank(3),size(MSI,3),Rank(1)],2);
D3 =  Gfold(D3,[Rank(3),size(HSI,3),Rank(1)],2); 

D11{3} = D3;D22{3} = D_3;

 %%  reconstruct
GG = []; GG{1} = D1; GG{2} = D2; GG{3} = D3;
HR_HSI = coreten2tr(GG);
HR_HSI(HR_HSI<0)=0;
HR_HSI(HR_HSI>1)=1;
HR_HSI = double(HR_HSI);
% rmse1=getrmse(double(S*255),double(HR_HSI*255));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err_x=norm(HR_HSI(:)-U(:))/norm(HR_HSI(:));
PSNR=csnr(S*255,HR_HSI*255,0,0);
     if i==1 || mod(i, 1) == 0
             fprintf('I=%d,psnr=%f,diff=%f\n',i,PSNR,err_x);
     end
    if psnr_max < PSNR
        psnr_max = PSNR;
        imgBest = HR_HSI;
%     else
%         break;
    end
    if (err_x<1e-4)     
        fprintf('Stop at %.0f\n',i); break
    end
    U = HR_HSI;
    
%     
    Relative_error(i) = err_x;
if i==30;
 figure;
plot(1:i, Relative_error(1:i)); 
 axis([0,30,0,0.08]);%确定x轴和y轴框图的大小 
 set(gca,'XTick',[0:5:30]); 
 set(gca,'YTick',[0:0.02:0.08]);
 xlabel('The Number of Iterations');
 ylabel('Reative error');
 grid on;
%  
%   i2= 1:i; 
%   LOGLRTR_error_iter_paviaC = Relative_error(1:i);
%  filename = 'LOGLRTR_error_iter_paviaC';
%  save(filename, 'i2', 'LOGLRTR_error_iter_paviaC');
% % 
% % %saveas(gcf, 'peppers(eror).png');   % 将图形保存为PNG文件
  end
% % %%%%%%%%%%%%%%%%%%%%%%%%%  PSNR vs iter   %%%%%%%%%%%%%%%%
psnr_values(i) = PSNR;
if i==30;
figure;
plot(1:i, psnr_values(1:i));

 axis([0,30,25,50]);%确定x轴和y轴框图的大小 
 set(gca,'XTick',[0:5:30]); 
 set(gca,'YTick',[25:5:50]);
 xlabel('The Number of Iterations');
 ylabel('PSNR');
 grid on;
%  
%  i2= 1:i;
% LOGLRTR_psnr_iter_paviaC= psnr_values(1:i);
% filename = 'LOGLRTR_psnr_iter_paviaC.mat';
% save(filename, 'i2', 'LOGLRTR_psnr_iter_paviaC');
% %saveas(gcf, 'psnr_iter_Yellow peppers.png');   % 将图形保存为PNG文件
 end 
% 
%     
    
    
    
end
runtime = toc;

end
