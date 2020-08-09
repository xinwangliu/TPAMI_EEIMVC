function WP = updateWPabsentClusteringV1Beta(HP,Hstar,beta,lambda)

k = size(HP,2);
num = size(HP,1);
numker = size(HP,3);
WP = zeros(k,k,numker);

HW = zeros(num,k);
for p = 1 : numker
  HW = HW + beta(p)*HP(:,:,p)*WP(:,:,p);
end

for p = 1 : numker
    Tp = HP(:,:,p)'*(Hstar - lambda*(HW - beta(p)*HP(:,:,p)*WP(:,:,p)));
    [Up,Sp,Vp] = svd(Tp,'econ');
    WP(:,:,p) = Up*Vp';
end