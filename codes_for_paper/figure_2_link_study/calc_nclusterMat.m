nClustMat_init = nan(nVar,nMC);
nClustMat_final = nan(nVar,nMC);

wbar = waitbar(0, 'Starting');

for iVar=1:1:nVar
    waitbar(iVar/nVar, wbar, sprintf('Progress: %d %%', floor(iVar/nVar*100)));
    
    parfor MC_ctr=1:nMC
        % Initial pos
        tmp1 = posArr(iVar,1,:,1,MC_ctr);
        tmp2 = posArr(iVar,2,:,1,MC_ctr);
        ps = [reshape(tmp1,[1,length(tmp1)]);reshape(tmp2,[1,length(tmp2)])];

        linkThresh = linkThreshArr(iVar);
        clustIDs = numCluster_rad(ps(1,:)',ps(2,:)',linkThresh);

        nClust = max(clustIDs);
        nClustMat_init(iVar,MC_ctr) = nClust;


        % final pos
        tmp1 = posArr(iVar,1,:,end,MC_ctr);
        tmp2 = posArr(iVar,2,:,end,MC_ctr);
        ps = [reshape(tmp1,[1,length(tmp1)]);reshape(tmp2,[1,length(tmp2)])];

        linkThresh = linkThreshArr(iVar);
        clustIDs = numCluster_rad(ps(1,:)',ps(2,:)',linkThresh);

        nClust = max(clustIDs);
        nClustMat_final(iVar,MC_ctr) = nClust;
    end
end
close(wbar)
save(mainStrng+".mat",'-v7.3');