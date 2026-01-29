
function [G0] = WST11(G1,Y1,Z1,P1,P2,U1,T1,M1,par,s0,sf,D)
lambda  =  par.lambda;
rho     =  par.rho;
mu    =  par.mu;
G0 = Gunfold(G1,2);


A1 = P1*P1';  DG=D'*D;   
A2 = lambda*(P2*P2')+rho*eye(size(P2,1));
H3 = PTx( Y1,U1,s0,sf )*P1' + lambda*Z1*P2' + rho*G0 + D'*Gunfold(mu*T1+M1,2);
r0 = H3 - PTPx(G0,U1,s0,sf)*A1 - G0*A2 - mu*DG*G0;

p0=r0;
maxIter =30;
for i=1:maxIter  
    pp = PTPx(p0,U1,s0,sf)*A1 + p0*A2 + mu*DG*p0 ;  
    pp1=p0(:)'*pp(:);
    a = (r0(:)')*r0(:)/pp1;
    G = G0+a*p0;
    r1 = r0-a*pp;
      b1=(r1(:)'*r1(:))/(r0(:)'*r0(:));
       p1=r1+b1*p0;
 p0=p1;
 r0=r1;
 G0=G;            
end
G0 = double(Gfold(G0,size(G1),2));
end

  
