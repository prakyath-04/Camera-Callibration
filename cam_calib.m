clear all;close all;
%% extracting the pixel co-ordinates of 4 same corners of 4 different checkerboard images
a = detectCheckerboardPoints('./images/img1.png');
b = detectCheckerboardPoints('./images/img3.png');
c = detectCheckerboardPoints('./images/img5.png');
d = detectCheckerboardPoints('./images/img4.png');
x1 = zeros(16,2);
x1(1:4,1:2) = [a(3,1) a(3,2);a(6,1) a(6,2);a(43,1) a(43,2);a(48,1) a(48,2)];
x1(5:8,1:2) = [b(3,1) b(3,2);b(6,1) b(6,2);b(43,1) b(43,2);b(48,1) b(48,2)];
x1(9:12,1:2) = [c(3,1) c(3,2);c(6,1) c(6,2);c(43,1) c(43,2);c(48,1) c(48,2)];
x1(13:16,1:2) = [d(3,1) d(3,2);d(6,1) d(6,2);d(43,1) d(43,2);d(48,1) d(48,2)];
pixels = x1';
%considering the world center at the bottom left corner of the checkeboard
real_pts = 0.01*[2.4*3 2.4;6 * 2.4 2.4;2.4 8*2.4;6*2.4 8*2.4]; % size of each block on the checkerboard = 2.4cm
X = real_pts;
%% finding Homography matrix
h = {};
for i=1:4
    h{i} = solve_h(x1((i-1)*4+1:4*i,:),X);
end
for i=1:4
    h{i} = h{i} ./ h{i}(3,3); %% normalize the Homography matrix
end
h1 = h{1};h2 = h{2};h3 = h{3};h4 = h{4};
%% finding 'B' matrix B = (inv(K)'*inv(K))
v = cell(3,2,2);
for k =1:3
    for i=1:2
        for j=1:2
            v{k,i,j} = [h{k}(1,i)*h{k}(1,j); h{k}(1,i)*h{k}(2,j) + h{k}(2,i)*h{k}(1,j);...
                h{k}(3,i)*h{k}(1,j) + h{k}(1,i)*h{k}(3,j); h{k}(2,i)*h{k}(2,j);...
                h{k}(3,i)*h{k}(2,j) + h{k}(2,i)*h{k}(3,j); h{k}(3,i)*h{k}(3,j)]';
        end
    end
end
V = [v{1,1,2};v{1,1,1}-v{1,2,2};v{2,1,2};v{2,1,1}-v{2,2,2};v{3,1,2};v{3,1,1}-v{3,2,2}];
aaa = v;
[~,~,v] = svd(V);
b = zeros(1,6); b(1,:) = v(:,6); B = zeros(3,3);
B(1,:) = b(1:3); B(2,1) = b(1,2); B(2,2) = b(1,4);B(2,3) = b(1,5);
B(3,1) = B(1,3); B(3,2) = B(2,3); B(3,3) = b(1,6);
%% finding matrix K
T1 = (B + B')/2;
tm = eig(T1);
T=T1;
mn = min(tm);
i=100000;
while i>0
    if min(eig(T))>0
        break;
    end
    T = T - (mn*mn + 10e-12) .* eye(3);
    T = (T + T')/2;
    tm = eig(T);
    mn = min(tm);
    i=i-1;
end
if(mn < 0)
    disp("loop for making the 'B' matrix positive definite takes a long time to converge");
end
T = (T + T')/2;
k= chol(T);
k = inv(k);
k = k./k(3,3);
disp(k);