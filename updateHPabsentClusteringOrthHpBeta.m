function HP = updateHPabsentClusteringOrthHpBeta(WP,Hstar,beta,S,HP00,HP_0,lambda)

num = size(Hstar,1);
k = size(Hstar,2);
numker = size(WP,3);
HP0 = zeros(num,k);
for p = 1:numker
    HP0 = HP0 + beta(p)*HP_0(:,:,p)*WP(:,:,p);
end

HP = zeros(num,k,numker);
for p = 1:numker
    Ap = Hstar - lambda*(HP0 - beta(p)*HP_0(:,:,p)*WP(:,:,p));
    
    mis_indx = S{p}.indx';
    obs_indx = setdiff(1:num,mis_indx);
    %%%%
    Vp = Ap(mis_indx,:)*WP(:,:,p)';
    [Up,Sp,Vp] = svd(Vp,'econ');
    HP(mis_indx,:,p) = Up*Vp';
    %%%%
    HP(obs_indx,:,p) = HP00(obs_indx,:,p);
end