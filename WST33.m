function [G] = WST33(G3,Y3,Z3,P1,P2,U3,T3,M3,par,D)
lambda  =  par.lambda;
rho     =  par.rho;
mu    =  par.mu;
G4 = Gunfold(G3,2);


A1 = P1*P1'+rho*eye(size(P1,1)); 
A2 = P2*P2'; B1=U3'*U3;   DG=D'*D;         
H3 = Y3*P1' + lambda*U3'*Z3*P2' + rho*G4  + D'*Gunfold(mu*T3+M3,2);
r0 = H3 - lambda*B1*G4*A2 - G4*A1 - mu*DG*G4;
p0=r0;
maxIter =30;
for i=1:maxIter
    pp = (lambda*B1*p0*A2 +  p0*A1 + mu*DG*p0);    
    pp1=p0(:)'*pp(:);
    a = (r0(:)')*r0(:)/pp1;
    G = G4+a*p0;
    r1 = r0-a*pp;
      b1=(r1(:)'*r1(:))/(r0(:)'*r0(:));
       p1=r1+b1*p0;
 p0=p1;
 r0=r1;
 G4=G;
end
G = double(Gfold(G,size(G3),2));
end
