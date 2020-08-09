clear
clc
warning off;

path = './';
addpath(genpath(path));

dataName = 'flower17'; %%% flower17; flower102; CCV; caltech101_numofbasekernel_10
%% %% washington; wisconsin; texas; cornell
load([dataName,'_Kmatrix'],'KH','Y');
% load([path,'datasets\',dataName,'_Kmatrix'],'KH','Y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numclass = length(unique(Y));
numker = size(KH,3);
num = size(KH,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KH = kcenter(KH);
KH = knorm(KH);

qnorm = 2;

epsionset = 0.1;
for ie =1:length(epsionset)
    for iter = 1:10
        load(['./generateAbsentMatrix/',dataName,'_missingRatio_',num2str(epsionset(ie)),...
            '_missingIndex_iter_',num2str(iter),'.mat'],'S');
        lambda = 2.^(-15:3:15);
        for lam = 1:length(lambda)
            H_normalized = incompleteLateFusionMKCOrthHp_lambda(KH,S,numclass,qnorm,lambda(lam));
            res(lam,:) = myNMIACC(H_normalized, Y, numclass);
        end
    end
end