function [G3] = WST3(G3,Y3,Z3,Q1,Q2,U3,par,D)

mu    =  par.mu;
eta      =  par.eta ;
U = zeros(size(G3));
M3 = U;N3 = U;
T3 = double(ttm(tensor(G3),D,2));

iter = 0;
while iter < 10
        iter = iter + 1; 
     
        [ G3] = WST33(G3,Y3,Z3,Q1,Q2,U3,T3,M3,par,D); 
         
        
        %J3 = double(ttm(tensor(G3),D,2)) - (M3./mu); 
        %W3 = 1./(abs(J3)+eps);
        W3 = ones(size(T3));
        
       [ H3] =  Log_prox_tnn2( W3.*T3 - (N3/mu) , eta/mu  );

          R3 = real(fft(W3, [], 2)) .* real(fft(H3 + N3 / mu, [], 2));
       T3 = real(ifft((real(fft(double(ttm(tensor(G3), D, 2)) - M3 / mu, [], 2)) + R3) ./ (1 +real( fft((W3).^2, [], 2))), [], 2)); 
      
        
    
        
       
        M3 = M3 + mu*(T3 - double(ttm(tensor(G3),D,2)));
        N3 = N3 + mu * (H3 - W3.*T3);
     

end

