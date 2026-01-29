function [X] = Log_prox_tnn(Y,rho)

[n1,n2,n3] = size(Y); 
X = zeros(n1,n2,n3);   

if (n1 <= n2 && n1 < n3)
   Y = fft(Y,[],1);  
        X(1,:,:) = WNNM(squeeze(Y(1,:,:)), rho, 1);
        halfn1 = round(n1/2);
        for i = 2 : halfn1
           X(i,:,:) = WNNM(squeeze(Y(i,:,:)), rho, 1);
            X(n1+2-i,:,:) = conj(X(i,:,:));
        end

        % if n3 is even
        if mod(n1,2) == 0
            i = halfn1+1;
            X(i,:,:) = WNNM(squeeze(Y(i,:,:)), rho, 1);
        end
        X = ifft(X,[],1);
end
   
if (n2 <= n1 && n2 <= n3)
   Y = fft(Y,[],2);  
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
        X = ifft(X,[],2);
end
if (n3 <= n2 && n3 <= n1)
   Y = fft(Y,[],3);  
        X(:,:,1) = WNNM( Y(:,:,1), rho, 1);
        % i=2,...,halfn3
        halfn3 = round(n3/2);
        for i = 2 : halfn3
           X(:,:,i) = WNNM( Y(:,:,i), rho, 1);
            X(:,:,n3+2-i) = conj(X(:,:,i));
        end

        % if n3 is even
        if mod(n3,2) == 0
            i = halfn3+1;
            X(:,:,i) = WNNM( Y(:,:,i), rho, 1);
        end
        X = ifft(X,[],3);
end