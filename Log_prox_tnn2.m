
function [X,tnn,trank] = Log_prox_tnn2(Y,rho)

% The proximal operator of the tensor nuclear norm of a 3 way tensor
%
% min_X rho*||X||_*+0.5*||X-Y||_F^2
%
% Y     -    n1*n2*n3 tensor
%
% X     -    n1*n2*n3 tensor
% tnn   -    tensor nuclear norm of X
% trank -    tensor tubal rank of X
%
% version 2.1 - 14/06/2018
%
% Written by Canyi Lu (canyilu@gmail.com)
%
%
% References: 
% Canyi Lu, Tensor-Tensor Product Toolbox. Carnegie Mellon University. 
% June, 2018. https://github.com/canyilu/tproduct.
%
% Canyi Lu, Jiashi Feng, Yudong Chen, Wei Liu, Zhouchen Lin and Shuicheng
% Yan, Tensor Robust Principal Component Analysis with A New Tensor Nuclear
% Norm, arXiv preprint arXiv:1804.03728, 2018
%

A = fft(Y,[],3);
[n1,n2,n3] = size(Y);    
X = zeros(n1,n2,n3);      
Y = fft(Y,[],2);        
tnn = 0;
trank = 0;
        
% first frontal slice


X(:,1,:) = WNNM( squeeze(Y(:,1,:)), rho, 1);
% i=2,...,halfn3
halfn2 = round(n2/2);
for i = 2 : halfn2
   X(:,i,:) = WNNM( squeeze(Y(:,i,:)), rho, 1);
    X(:,n2+2-i,:) = conj(X(:,i,:));
end

% if n3 is even
if mod(n2,2) == 0
    i = halfn2+1;
    X(:,i,:) = WNNM( squeeze(Y(:,i,:)), rho, 1);
end
tnn = tnn/n2;
X = ifft(X,[],2);

