function h = solve_h(x1,X)
ax = zeros(9,4);
ay = zeros(9,4);
t1 = [-1 0 0 0];
t2 = [0 0 0];
ax(3:6,:) = [t1;t1;t1;t1]';
ay(1:3,:) = [t2;t2;t2;t2]';
ay(6,:) = [-1 -1 -1 -1];
for i=1:4
    ax(1:2,i) = -[X(i,1) X(i,2)];ax(7:9,i) = [x1(i,1)*X(i,1) x1(i,1)*X(i,2) x1(i,1)];
    ay(4:5,i) = -[X(i,1) X(i,2)];ay(7:9,i) = [x1(i,2)*X(i,1) x1(i,2)*X(i,2) x1(i,2)];
end
a = [ax ay];
a = a*a'; % a square matrix 
[~,~,v] = svd(a);
h = (reshape(v(:,9),3,3)).';
end