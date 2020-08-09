function [beta,obj] = updateBetaAbsentClusteringBeta(HP,WP,Hstar,lambda)

numker = size(HP,3);
M = zeros(numker);
for p=1:numker
    for q = p:numker
%         M(p,q) = trace(WP(:,:,p)'*HP(:,:,p)'*HP(:,:,q)*WP(:,:,q))/sqrt(trace((WP(:,:,p)'*HP(:,:,p)'*HP(:,:,p)*WP(:,:,p)))...
%             *trace(WP(:,:,q)'*HP(:,:,q)'*HP(:,:,q)*WP(:,:,q)) );
        M(p,q) = trace(WP(:,:,p)'*HP(:,:,p)'*HP(:,:,q)*WP(:,:,q));
    end
end
M = (M+M')-diag(diag(M))+1e-8*eye(numker);
f = zeros(numker,1);
for p=1:numker
    f(p) = -trace(HP(:,:,p)'*(Hstar*WP(:,:,p)'))/lambda;
end

% [beta,obj] = solvingV(M,f);

A = -eye(numker);
b = zeros(numker,1);
Aeq = ones(numker,1)';
beq = 1;
LB = zeros(numker,1);
UB = ones(numker,1);

[beta,obj]= quadprog(M,f,A,b,Aeq,beq,LB,UB);
beta((beta<eps))=0;
beta = beta/sum(beta);