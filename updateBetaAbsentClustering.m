function beta = updateBetaAbsentClustering(HP,WP,Hstar,qnorm)
%% function beta = updateBetaAbsentClustering(HP,WP,Hstar,HP00,lambda)
numker = size(WP,3);
HHPWP = zeros(numker,1);
for  p=1:numker
    HHPWP(p) = trace(Hstar'*(HP(:,:,p)*WP(:,:,p)));
end
% beta = HHPWP./norm(HHPWP);

beta = HHPWP.^(1/(qnorm-1))/sum(HHPWP.^(qnorm/(qnorm-1)))^(1/qnorm);

% numker = size(HP,3);
% M = zeros(numker);
% for p=1:numker
%     for q = p:numker
%         M(p,q) = trace(WP(:,:,p)'*HP(:,:,p)'*HP(:,:,q)*WP(:,:,q));
%     end
% end
% M = (M+M')-diag(diag(M))+1e-8*eye(numker);
% f = zeros(numker,1);
% for p=1:numker
%     f(p) = - trace(HP(:,:,p)'*(Hstar*WP(:,:,p)'+lambda*HP00(:,:,p)));
% end
% 
% A = -eye(numker);
% b = zeros(numker,1);
% Aeq = ones(numker,1)';
% beq = 1;
% LB = zeros(numker,1);
% UB = ones(numker,1);
% 
% beta = quadprog(M,f,A,b,Aeq,beq,LB,UB);
% beta((beta<eps))=0;
% beta = beta/sum(beta);