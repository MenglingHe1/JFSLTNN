function [X] = Log_prox_tnn3(Y,rho)

[n1,n2,n3] = size(Y); 
X = zeros(n1,n2,n3);   

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


