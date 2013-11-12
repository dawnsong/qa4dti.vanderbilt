
Reg_Im=Reg_ImA;
%h=waitbar(0,'Registering Images');
affMat=zeros(4,4,Ng);
trackMove=zeros(1,Ng);
rotation=zeros(Ng,3); translation=rotation;
for vol=1:Ng
    dwi=make_nii(Y_data(:,:,:,vol),resolution,[0 0 0],16); dwi.hdr.dime.xyzt_units=2; name_dwi = sprintf('%s%sdwi.nii',tmp,filesep);
    save_nii(dwi,name_dwi)  %convert to nii file for FLIRT
    name_temp = sprintf('%s%stemp',tmp,filesep); name_aff = sprintf('%s%saff',tmp,filesep);
    
    cmmd=sprintf('!flirt -in %s -ref %s -out %s -searchrx -10 10 -searchry -10 10 -searchrz -10 10 -omat %s', name_dwi, name_bo_ref,name_temp, name_aff);
    eval(cmmd)
    
    cmmd=sprintf('avscale --allparams %s %s',name_aff,name_dwi);
    [ah avs]=system(cmmd);
    [translation(vol,:) rotation(vol,:)]=readAVS(avs);
    
    cmmd=sprintf('!gunzip %s.nii.gz',name_temp);
    eval(cmmd)
    
    load(name_aff); affMat(:,:,vol)=aff; %convert affine matrix FLIRT output to matlab variable affMat
    name_temp_nii=sprintf('%s.nii',name_temp);
    T1=load_nii(name_temp_nii); Reg_Im(:,:,:,vol)=T1.img; %convert image output to matlab variable Reg_im
    
    cmmd=sprintf('!rm %s %s', name_temp_nii,name_aff);
   eval(cmmd)
    %waitbar(vol/Ng,h);
end
%close(h)
%-coarsesearch 15 -finesearch 6 