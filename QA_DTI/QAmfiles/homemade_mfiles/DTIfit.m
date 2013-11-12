function [FAsim ]=DTIfit(vox,nosim,dwi,bo,sigma_hat,b,g,n_bo)

%dwibreak = dwi; bobreak = bo; g = gtable;     
     
     FAsim = zeros(vox,1);
  
  
    for k=1:nosim

        Si = ricernd(dwi,sigma_hat); %...... [#grad_dir by #voxels] 
        So = ricernd(bo(1,:), sigma_hat/sqrt(n_bo)); So=repmat(So,size(dwi,1),1);  % ..... [#grad_dir by #voxels]  should be 1 by #voxels, but make bigger for easy division
    
        b2=repmat(b',1,vox);
        D = g\(-log(Si./So)./b2);  % [ 6 by #voxels]
        D = reshape(D,6,1,[]); % [6 by 1 by #voxels]
        Dtensor = [D(1,1,:) D(4,1,:) D(5,1,:); D(4,1,:) D(2,1,:) D(6,1,:); D(5,1,:) D(6,1,:) D(3,1,:)]; % [3 by 3 by #voxels];   
      
         for m=1:vox;
             tstartm=tic;
             if sum(sum(isnan(Dtensor(:,:,m)))) ==0
                [Vec Lambda]=eig(Dtensor(:,:,m)); 
                elle=diag(Lambda)'; index=[1 2 3]; [selle sIndex]= sortBbyA(elle,index);
                sIndex=sIndex(end:-1:1);
                L = [Lambda(sIndex(1),sIndex(1)) Lambda(sIndex(2),sIndex(2)) Lambda(sIndex(3),sIndex(3))];% lambda for this voxel for this simulation number
                MD = (sum(L)/3); % MD for this voxel this simulation number
                FA = sqrt(3/2)* (sqrt(sum((L-MD).^2)) /sqrt(sum(L.^2))); % FA for this voxel this simulation number
                FAsim(m) = FAsim(m) +FA/nosim; % FA for this voxel averaged over all simulations
             end
         end
    end
  

  