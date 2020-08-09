function [H_normalized,WP,HP,beta,obj] = incompleteLateFusionMKCOrthHp_lambda(KH,S,k,qnorm,lambda)

num = size(KH, 2); %the number of samples
numker = size(KH, 3); %m represents the number of kernels
maxIter = 100; %the number of iterations
[HP,WP] = myInitializationHp(KH,S,k);
HP00 = HP;
beta = ones(numker,1)*(1/numker)^(1/qnorm);
%%%---H0
KA = feval('algorithm2',KH,S);
%% combining the base kernels
KC  = mycombFun(KA,beta);
H0 = mykernelkmeans(KC,k);

%%%%%%%%%%
flag = 1;
iter = 0;
RpHpwp = zeros(num,k); % k - clusters, N - samples
for p=1:numker
    RpHpwp = RpHpwp +  beta(p)*(HP(:,:,p)*WP(:,:,p));
end
RpHpwp_lambda = RpHpwp +lambda*H0;  
while flag
    iter = iter +1;
    %---the first step-- optimize H_star with given (HP, WP and beta)
    [Uh,Sh,Vh] = svd(RpHpwp_lambda,'econ');
    Hstar = Uh*Vh';
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     res90(iter,:) = myNMIACC(Hstar,Y,k);
    
    %---the second step-- optimize WP with (HP, H_star and beta)
    WP = updateWPabsentClusteringV1(HP,Hstar);
    
    %---the third step-- optimize HP with (WP, H_star and beta)
    HP = updateHPabsentClusteringOrthHp(WP,Hstar,S,HP00);
    
    %---the fourth step-- optimize beta with (WP, H_star and HP)
    beta = updateBetaAbsentClustering(HP,WP,Hstar,qnorm);
    
    %---Calculate Obj--
    RpHpwp = zeros(num,k);
    for p = 1:numker
        RpHpwp = RpHpwp + beta(p)*HP(:,:,p)*WP(:,:,p);
    end
    RpHpwp_lambda = RpHpwp +lambda*H0;
    obj(iter) = trace(Hstar'*RpHpwp_lambda);
    if (iter>2) && (abs((obj(iter)-obj(iter-1))/(obj(iter)))<1e-4 || iter>maxIter)
        flag =0;
    end
%     if (iter>2) && (iter>maxIter)
%         flag =0;
%     end
end
H_normalized = Hstar./ repmat(sqrt(sum(Hstar.^2, 2)), 1,k);