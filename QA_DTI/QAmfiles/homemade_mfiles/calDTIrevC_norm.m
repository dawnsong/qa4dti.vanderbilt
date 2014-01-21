function [chi_sq_p Se Errors]=calDTIrevC_norm(T,grad_file,bvals,S_measured,mask)
%normalized according to 2003 "A measure of curve fitting error for noise
%filtering diffusion tensor MRI data"
Nx=size(S_measured,1); Ny=size(S_measured,2); Nz=size(S_measured,3);
bo=reshape(S_measured(:,:,:,end),1,[]);
 
 D_vec=[];
  
    D_vec(1,:)=reshape(T(:,:,:,1),1,[]);  %Dxx
    D_vec(2,:)=reshape(T(:,:,:,3),1,[]); %Dyy
    D_vec(3,:)=reshape(T(:,:,:,6),1,[]); %Dzz
    D_vec(4,:)=reshape(T(:,:,:,2),1,[]);%Dxy
    D_vec(5,:)=reshape(T(:,:,:,4),1,[]);%Dxz
    D_vec(6,:)=reshape(T(:,:,:,5),1,[]); %Dyz

%create gradient vector
grads=grad_file';
gs=[grads.^2  2*grads(:,1).*grads(:,2) 2*grads(:,1).*grads(:,3) 2*grads(:,3).*grads(:,2)];

S_model=[];
for vol=1:size(gs,1)-1
    gb(vol,1:6)=bvals(vol)*gs(vol,:);
end
Sm=S_measured(:,:,:,1:end-1)./repmat(S_measured(:,:,:,end),[1 1 1 size(S_measured,4)-1]);
Sm=reshape(Sm,Nx*Ny*Nz,size(gs,1)-1); Sm=permute(Sm,[2 1]);
Se=exp(-gb*D_vec);
Errors=Sm-Se;
    ee=repmat(bo,size(gs,1)-1,1).*Se;

S_model=reshape(permute(ee,[2 1]),Nx,Ny,Nz,size(gs,1)-1);
delta_S=(S_model-S_measured(:,:,:,1:end-1)).^2;
S_fi=S_measured(:,:,:,1:end-1);
normS=sum(S_fi.^2,4);
chi_sq_p=zeros(Nx*Ny*Nz,size(gs,1)-1);
mask=mask(:); 
mask=(mask==1);
delta_S=reshape(delta_S,Nx*Ny*Nz,size(gs,1)-1);
normS=reshape(normS,Nx*Ny*Nz,1);
for vol=1:size(gs,1)-1
    chi_sq_p(mask,vol)=delta_S(mask,vol)./normS(mask);
end



% sd(Sd>=1)=NaN;
% sd(Sd==Inf)=NaN;
% warning=zeros(size(sd));
% warning(Sd>=1)=1;
% warning(Sd==Inf)=1;



% % % 
% % % 
