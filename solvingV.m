function [x,obj,lambda] = solvingV(A,b)
%% min_x 1/2*x'*A*x - b'*x, s.t. x'*x =1.

A = (A+A')/2;
num = size(A,1);
% v0 = ones(num,1)*(1/sqrt(num));
% flag =1;
% iter = 0;
% while flag
%     iter = iter +1;
%     u = A*v0 + b;
%     v = u/norm(u);
%     obj(iter) = 0.5*(v'*A*v) + b'*v;
%     if (iter>2 && abs(obj(iter-1)-obj(iter))/obj(iter)<1e-4)
%         flag = 0;
%     end
%     v0 = v;
% end
G = [A,-eye(num); -(b*b'),A];
% G = (G+G')/2;
[D,Lambda] = eig(G);
diagL = diag(Lambda);
indx = find(abs(imag(diagL))<eps);
lambda = min(diagL(indx));
x = (A-lambda*eye(num))\b;
obj = (1/2)*x'*A*x - b'*x;
