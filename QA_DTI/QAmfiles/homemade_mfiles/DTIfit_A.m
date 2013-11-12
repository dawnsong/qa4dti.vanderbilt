function [FA check]=DTIfit_A(vox,data,b,g)
         FA=zeros(vox,1);FA=FA./0;
        b2=repmat(b',1,vox);
        D = g\(-log(data)./b2);  % [ 6 by #voxels]
        D = reshape(D,6,1,[]); % [6 by 1 by #voxels]
        Dtensor = [D(1,1,:) D(4,1,:) D(5,1,:); D(4,1,:) D(2,1,:) D(6,1,:); D(5,1,:) D(6,1,:) D(3,1,:)]; % [3 by 3 by #voxels];   
        check=zeros(vox,1);
         for m=1:vox;
             if sum(sum(isnan(Dtensor(:,:,m)))) ==0
                 check(m)=1;
                [Vec Lambda]=eig(Dtensor(:,:,m)); 
                elle=diag(Lambda)'; index=[1 2 3]; [selle sIndex]= sortBbyA(elle,index);
                sIndex=sIndex(end:-1:1);
                L = [Lambda(sIndex(1),sIndex(1)) Lambda(sIndex(2),sIndex(2)) Lambda(sIndex(3),sIndex(3))];% lambda for this voxel for this simulation number
                MD = (sum(L)/3); % MD for this voxel this simulation number
                FA(m) = sqrt(3/2)* (sqrt(sum((L-MD).^2)) /sqrt(sum(L.^2))); % FA for this voxel this simulation number
                
             end
           
         end
 

  