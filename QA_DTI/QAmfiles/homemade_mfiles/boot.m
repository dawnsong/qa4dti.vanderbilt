function FAboot=boot(ModelData,Errors,b,g,bootNum,FAobs)


b=b(1:end-1);
grads=g';
gtable=[grads.^2  2*grads(:,1).*grads(:,2) 2*grads(:,1).*grads(:,3) 2*grads(:,3).*grads(:,2)];
vox=size(ModelData,2);
% % % % % % % % % % % % % % % % % 
% Resample Errors, populate FA 
% % % % % % % % % % % % % % % % % 
FA = zeros(bootNum+1,size(FAobs,2)); FA(1,:)=FAobs;
BootData=zeros(size(ModelData));
checkA=zeros(size(FA)); checkA(1,:)=ones(1,size(FA,2));

try
    matlabpool open
    parfor boot=2:bootNum+1
           neg = (randn(size(ModelData))<0); %50% chance positive, 50% chance negative
           randsign=ones(size(ModelData));
           randsign(neg)=-1;
        save(sprintf('~/tmp/boot.bf.%d.mat', boot))
        for m=1:vox
            BootData(:,m) = ModelData(:,m)+(randsample(Errors(:,m),size(ModelData,1),'true').*randsign(:,m)); 
        end
        save(sprintf('~/tmp/boot.af.%d.mat', boot))

      
       %FA(boot,:)=DTIfit_A(vox,BootData,b,gtable);
       FA(boot,:)=dawn_DTIfit_A(vox,BootData,b,gtable);
       
    end

except err,
end
matlabpool close


FAboot=FA;
    
% % %     D = gtable\(-log(ModelData)/1000);   % [ 6 by #voxels]
% % %   
% % %     D = reshape(D,6,1,[]); % [6 by 1 by #voxels]
% % %     Dtensor = [D(1,1,:) D(4,1,:) D(5,1,:); D(4,1,:) D(2,1,:) D(6,1,:); D(5,1,:) D(6,1,:) D(3,1,:)]; % [3 by 3 by #voxels];   
% % %     for m=1:size(ModelData,2);
% % %         if sum(sum(isnan(Dtensor(:,:,m)))) == 0
% % %             checkA(boot,m)=1;
% % %             [Vec Lambda]=eig(Dtensor(:,:,m)); 
% % %             elle=diag(Lambda)'; index=[1 2 3]; [selle sIndex]= sortBbyA(elle,index);
% % %             sIndex=sIndex(end:-1:1);
% % %             L = [Lambda(sIndex(1),sIndex(1)) Lambda(sIndex(2),sIndex(2)) Lambda(sIndex(3),sIndex(3))];% lambda for this voxel for this simulation number
% % %             MD = (sum(L)/3); % MD for this voxel this simulation number
% % %             FA(boot,m) = sqrt(3/2)* (sqrt(sum((L-MD).^2)) /sqrt(sum(L.^2))); % FA for this voxel this simulation number
% % %         end
% % %         
% % %     end
% % % end
