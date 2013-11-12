
function downSampMask = downSamp_full(maskROI,voxNum)
% subsample voxels from each of the regions of interest for simex
%% takes input voxNum and estiamtes what fraction of total voxels that is. then randomly samples from maskROI at that percentage rate
totvox = sum(sum(maskROI));
perc = voxNum/totvox;

RN = rand(size(maskROI)); RN2 = (RN<perc);

 downSampMask=(maskROI==1 & RN2==1);

%% downSampMask, each column is one ROI
test=sum(downSampMask,1);  t=(test<50 & sum(maskROI,1)>50); 
while sum(t)>0
    remaining_perc=(50-test)./sum(maskROI,1);  remaining_perc(remaining_perc<0)=0;%number voxels needed/number voxelsin region
    perc_max=perc+remaining_perc; %new maximum percentage, unique for each ROI, if ROI<50, then perc>1 and all will be chosen
    for k=1:25
        RN2(:,k) = (RN(:,k)<perc_max(k));
    end
    downSampMask=(maskROI==1 & RN2==1);
    test=sum(downSampMask,1); t=(test<50 & sum(maskROI,1)>50);
    perc=perc_max;

end

t=(sum(maskROI,1)<51);  %any regions with 50 or less voxels
t=repmat(t,size(maskROI,1),1); t=(maskROI==1 & t==1); downSampMask=downSampMask+t;
downSampMask=(downSampMask>0);

