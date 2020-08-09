function KH2 = algorithm2(KH,S)

num = size(KH,1);
numker = size(KH,3);
KH2 = zeros(num,num,numker);
for p =1:numker
    %% missing index: S{p}.indx
    obs_indx = setdiff(1:num,S{p}.indx');
    KAp = KH(obs_indx,obs_indx,p);
    KH2(obs_indx,obs_indx,p) = (KAp+KAp')/2;
end
clear KH