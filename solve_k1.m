% function [K] = solve_k1(x1,X,H_sol)
I = eye(3,3);

detH = det(H_sol); 
H_sol = H_sol/nthroot(detH,3); % normalize the homography matrix 
% so that the determinate is 1 

[egvec,egval] = eig(H_sol); 
angle = asin(imag(egval(2,2))); 
s = sin(angle); 
c = cos(angle); 

A1 = H_sol - I; 
A2 = [H_sol-c*I, -s*I;s*I, H_sol-c*I]; 

[~,~,V] = svd(A1); 
h1 = V(:,3); 
[~,~,V] = svd(A2); 
h23 = V(:,6); 
h2 = h23(1:3); 
h3 = h23(4:6); 

HH = [h1,h2,h3]; 

[UMatrix,Mupper] = qr(inv(HH)); 
KK = inv(Mupper); 
RR = inv(UMatrix); 
r11 = RR(1,1); 
r21 = RR(2,1); 
k11 = KK(1,1); 
k12 = KK(1,2); 
nom = k11*r11*r21 - k12*r11*r11; 
denom = k11*r11*r21 - k12*(r11*r11-1); 
alpha = sqrt(abs(nom/denom)); 
diaga = eye(3,3); 
diaga(1,1) = alpha; 
Hnew = HH*diaga; 

[UMatrix,Mupper] = qr(inv(Hnew)); 
KK = inv(Mupper); 
camM = KK/KK(3,3); 
camM(1,1)= abs(camM(1,1)); 
camM(2,2)= abs(camM(2,2)); 
% end