%get median bias and standard deviation values to represent entire region

Y=load_nii(segs{1}); labelBo=Y.img; labelBo=labelBo(:); 
Bm=zeros(size(labelBo)); Sm=Bm;
clear segs

for k=1:25
    Bm(labelBo==k)=Broi_b(1,k);
    Sm(labelBo==k)=Fboot(1,k);
end

Bm=reshape(Bm,Nx,Ny,Nz); Sm=reshape(Sm,Nx,Ny,Nz);


