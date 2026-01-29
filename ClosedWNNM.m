function [SigmaX,svp]=ClosedWNNM(SigmaY,C,oureps)
temp=(SigmaY-oureps).^2-4*(C-oureps*SigmaY);  %C2
ind=find (temp>0);
svp=length(ind);
SigmaX= sign(SigmaY(ind)).*max(SigmaY(ind)-oureps+sqrt(temp(ind)),0)/2;

end
