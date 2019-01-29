clear all;
addpath(genpath('./images'));

[ip1,b1] = detectCheckerboardPoints('img1.png');
[ip2,b2] = detectCheckerboardPoints('img3.png');
[ip3,b3] = detectCheckerboardPoints('img5.png');
obj   = [2.4*2 2.4;6 * 2.4 2.4;2.4 8*2.4;6*2.4 8*2.4];

img_1 = [ip1(2,1) ip1(2,2);ip1(6,1) ip1(6 ,2);ip1(43,1) ip1(43,2);ip1(48,1) ip1(48,2)];
img_2 = [ip2(2,1) ip2(2,2);ip2(6,1) ip2(6 ,2);ip1(43,1) ip2(43,2);ip2(48,1) ip2(48,2)];
img_3 = [ip3(2,1) ip3(2,2);ip3(6,1) ip3(6 ,2);ip1(43,1) ip3(43,2);ip3(48,1) ip3(48,2)];
ax_1  = horzcat(-1*obj(1,:),-1,0,0,0,img_1(1,1)*horzcat(obj(1,:),1));
ax_2  = horzcat(-1*obj(2,:),-1,0,0,0,img_1(2,1)*horzcat(obj(2,:),1));
ax_3  = horzcat(-1*obj(3,:),-1,0,0,0,img_1(3,1)*horzcat(obj(3,:),1));
ax_4  = horzcat(-1*obj(4,:),-1,0,0,0,img_1(4,1)*horzcat(obj(4,:),1));
ay_1  = horzcat(0,0,0,-1*obj(1,:),-1,img_1(1,2)*horzcat(obj(1,:),1));
ay_2  = horzcat(0,0,0,-1*obj(2,:),-1,img_1(2,2)*horzcat(obj(2,:),1));
ay_3  = horzcat(0,0,0,-1*obj(3,:),-1,img_1(3,2)*horzcat(obj(3,:),1));
ay_4  = horzcat(0,0,0,-1*obj(4,:),-1,img_1(4,2)*horzcat(obj(4,:),1));
M = [ax_1;ay_1;ax_2;ay_2;ax_3;ay_3;ax_4;ay_4];
[U,S,V] = svd(M);
H1 = V(:,end);
h1 = [H1(1) H1(2) H1(3) ;H1(4) H1(5) H1(6) ;H1(7) H1(8) H1(9)];
ax_1  = horzcat(-1*obj(1,:),-1,0,0,0,img_2(1,1)*horzcat(obj(1,:),1));
ax_2  = horzcat(-1*obj(2,:),-1,0,0,0,img_2(2,1)*horzcat(obj(2,:),1));
ax_3  = horzcat(-1*obj(3,:),-1,0,0,0,img_2(3,1)*horzcat(obj(3,:),1));
ax_4  = horzcat(-1*obj(4,:),-1,0,0,0,img_2(4,1)*horzcat(obj(4,:),1));
ay_1  = horzcat(0,0,0,-1*obj(1,:),-1,img_2(1,2)*horzcat(obj(1,:),1));
ay_2  = horzcat(0,0,0,-1*obj(2,:),-1,img_2(2,2)*horzcat(obj(2,:),1));
ay_3  = horzcat(0,0,0,-1*obj(3,:),-1,img_2(3,2)*horzcat(obj(3,:),1));
ay_4  = horzcat(0,0,0,-1*obj(4,:),-1,img_2(4,2)*horzcat(obj(4,:),1));
M = [ax_1;ay_1;ax_2;ay_2;ax_3;ay_3;ax_4;ay_4];
[U,S,V] = svd(M);
H2 = V(:,end);
h2 = [H2(1) H2(2) H2(3) ;H2(4) H2(5) H2(6) ;H2(7) H2(8) H2(9) ];
ax_1  = horzcat(-1*obj(1,:),-1,0,0,0,img_3(1,1)*horzcat(obj(1,:),1));
ax_2  = horzcat(-1*obj(2,:),-1,0,0,0,img_3(2,1)*horzcat(obj(2,:),1));
ax_3  = horzcat(-1*obj(3,:),-1,0,0,0,img_3(3,1)*horzcat(obj(3,:),1));
ax_4  = horzcat(-1*obj(4,:),-1,0,0,0,img_3(4,1)*horzcat(obj(4,:),1));
ay_1  = horzcat(0,0,0,-1*obj(1,:),-1,img_3(1,2)*horzcat(obj(1,:),1));
ay_2  = horzcat(0,0,0,-1*obj(2,:),-1,img_3(2,2)*horzcat(obj(2,:),1));
ay_3  = horzcat(0,0,0,-1*obj(3,:),-1,img_3(3,2)*horzcat(obj(3,:),1));
ay_4  = horzcat(0,0,0,-1*obj(4,:),-1,img_3(4,2)*horzcat(obj(4,:),1));
M = [ax_1;ay_1;ax_2;ay_2;ax_3;ay_3;ax_4;ay_4];
[U,S,V] = svd(M);
H3 = V(:,end);
h3 = [H3(1) H3(2) H3(3) ;H3(4) H3(5) H3(6) ;H3(7) H3(8) H3(9) ];

% creating linear equation of V.b = 0
h = h1;
v11 = [h(1,1)*h(1,1) , h(1,1)*h(2,1) + h(2,1)*h(1,1) , h(3,1)*h(1,1)+h(1,1)*h(3,1) , h(2,1)*h(2,1) , h(3,1)*h(2,1) + h(2,1)*h(3,1) , h(3,1)*h(3,1)];
   
v22 = [h(1,2)*h(1,2) , h(1,2)*h(2,2) + h(2,2)*h(1,2) , h(3,2)*h(1,2)+h(1,2)*h(3,2) , h(2,2)*h(2,2) , h(3,2)*h(2,2) + h(2,2)*h(3,2) , h(3,2)*h(3,2)];
   
v12 = [h(1,1)*h(1,2) , h(1,1)*h(2,2) + h(2,1)*h(1,2) , h(3,1)*h(1,2)+h(1,1)*h(3,2) , h(2,1)*h(2,2) , h(3,1)*h(2,2) + h(2,1)*h(3,2) , h(3,1)*h(3,2)];

V_1 = [v12;v11 - v22 ];
% h2
h = h2;
v11 = [h(1,1)*h(1,1) , h(1,1)*h(2,1) + h(2,1)*h(1,1) , h(3,1)*h(1,1)+h(1,1)*h(3,1) , h(2,1)*h(2,1) , h(3,1)*h(2,1) + h(2,1)*h(3,1) , h(3,1)*h(3,1)];
   
v22 = [h(1,2)*h(1,2) , h(1,2)*h(2,2) + h(2,2)*h(1,2) , h(3,2)*h(1,2)+h(1,2)*h(3,2) , h(2,2)*h(2,2) , h(3,2)*h(2,2) + h(2,2)*h(3,2) , h(3,2)*h(3,2)];
   
v12 = [h(1,1)*h(1,2) , h(1,1)*h(2,2) + h(2,1)*h(1,2) , h(3,1)*h(1,2)+h(1,1)*h(3,2) , h(2,1)*h(2,2) , h(3,1)*h(2,2) + h(2,1)*h(3,2) , h(3,1)*h(3,2)];

V_2 = [v12;v11 - v22 ];

% h3
h = h3;
v11 = [h(1,1)*h(1,1) , h(1,1)*h(2,1) + h(2,1)*h(1,1) , h(3,1)*h(1,1)+h(1,1)*h(3,1) , h(2,1)*h(2,1) , h(3,1)*h(2,1) + h(2,1)*h(3,1) , h(3,1)*h(3,1)];
   
v22 = [h(1,2)*h(1,2) , h(1,2)*h(2,2) + h(2,2)*h(1,2) , h(3,2)*h(1,2)+h(1,2)*h(3,2) , h(2,2)*h(2,2) , h(3,2)*h(2,2) + h(2,2)*h(3,2) , h(3,2)*h(3,2)];
   
v12 = [h(1,1)*h(1,2) , h(1,1)*h(2,2) + h(2,1)*h(1,2) , h(3,1)*h(1,2)+h(1,1)*h(3,2) , h(2,1)*h(2,2) , h(3,1)*h(2,2) + h(2,1)*h(3,2) , h(3,1)*h(3,2)];

V_3 = [v12;v11 - v22 ];

% svd(V)
V_m = [V_1;V_2;V_3];
[U,S,V] = svd(V_m);
b = V(:,end);
B = [b(1) b(2) b(3) ; b(2) b(4) b(5); b(3) b(5) b(6)];
disp(B);
T1 = (B + B')/2;
% tm = eig(T1);
% T=T1;
% mn = min(tm);
% while true
%     if min(eig(T))>0
%         break;
%     end
%     if mn ==0
%         T = T + 10e-10 * eye(3);
%         break;
%     else
%         T = T + mn .* mn .*eye(3);
%         T = (T + T')/2;
%         tm = eig(T);
%         mn = min(tm);
%     end
% end
k= chol(T1);
k = inv(k);
k(:,:) = k(:,:) /k(3,3);
disp(k);