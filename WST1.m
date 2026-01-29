
function [G1] = WST1(G1,Y1,Z1,P1,P2,U1,s0,sf,par,D)
mu    =  par.mu;
eta      =  par.eta ;
U = zeros(size(G1));
M1 = U;N1 = U;
T1 = double(ttm(tensor(G1),D,2));

iter = 0;
while iter < 10
        iter = iter + 1; 
     
        [ G1] = WST11(G1,Y1,Z1,P1,P2,U1,T1,M1,par,s0,sf,D); 
        
        
       % J1 = double(ttm(tensor(G1),D,2)) - (M1./mu);     
        %W1 = 1./(abs(J1)+eps);
        W1 = ones(size(T1));
      [ H1] =  Log_prox_tnn2( W1.*T1 - (N1/mu) , eta/mu  );
        
      
      R1 = real(fft(W1, [], 2) .* real(fft(H1 + N1 / mu, [], 2)));
       T1 =real(ifft((real(fft(double(ttm(tensor(G1), D, 2)) - M1 / mu, [], 2) )+ R1) ./ (1 + real(fft((W1).^2, [], 2))), [], 2));  
         
        
      M1 = M1 + mu*(T1 - double(ttm(tensor(G1),D,2)));
        N1 = N1 + mu * (H1 - W1.*T1);
        
 
end


    