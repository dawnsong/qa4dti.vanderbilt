function [scan_info_text grad_file bval_vec resolution name_Y_data]= getNIIinfo(tmp,file_name)

Y =load_nii(file_name);
Ydata=Y.img; 

%read grad and bvalues file
[P,Fname,E]=fileparts(file_name);
if isempty(P)==1
    grad_name=sprintf('%s_bvecs',Fname);
    b_name=sprintf('%s_bvals',Fname);
else
    grad_name=sprintf('%s%s%s_bvecs',P,filesep,Fname);
    b_name=sprintf('%s%s%s_bvals',P,filesep,Fname);
end
grad_file=load(grad_name,'-ascii'); 
if size(grad_file,1)==3 % make grad file three columns
    grad_file=grad_file'; 
end
bval_vec=load(b_name,'-ascii');
if size(bval_vec,1)>1       % bvector one row
    bval_vec=bval_vec';
end

%% identify all images that are dwi and identify Bo (e.g. remove averaged image)
norm_gf=zeros(1,size(grad_file,1));
for k=1:size(grad_file,1)
    norm_gf(k) = norm(grad_file(k,:));
end
dwi_loc =(norm_gf >.98 & norm_gf < 1.02);
bo_loc = (bval_vec==0);

%% load dwi, place Bo as last volume
Y_data = Ydata(:,:,:,dwi_loc);
Y_data(:,:,:,end+1)=Ydata(:,:,:,bo_loc);

%% fix gradient table and bo to match identified dwi and Bo as last volume
gf=grad_file(dwi_loc,:); gf(end+1,:)=grad_file(bo_loc,:);
bv=bval_vec(dwi_loc); bv(end+1)=bval_vec(bo_loc);
bval_vec=bv;

%% fix gradient table to match changes made to image upon loading into matlab
rot=Y.hdr.hist.rot_orient;
flip=Y.hdr.hist.flip_orient;
[t s]=sort(rot);
%dawnsong
if isempty(s),
    grad_file=gf;
else
    grad_file=gf(:,s);
end
k=(flip~=0);
grad_file(:,k)=-grad_file(:,k);
   
%% save data as RAS-NII file
name_Y_data=sprintf('%s%sY_data.nii',tmp,filesep);
resolution=Y.hdr.dime.pixdim(2:4); 
nii=make_nii(Y_data,resolution,[0 0 0],16); nii.hdr.dime.xyzt_units=10; save_nii(nii,name_Y_data);

%% Parse header
scan_info_text = {};
        scan_info_text{end+1} = sprintf('Exam File Name:  %s', Fname );
        scan_info_text{end+1} = sprintf('Exam Header Name:  %s', Y.hdr.hk.db_name);
        scan_info_text{end+1} = sprintf('Header Descriptiton:  %s', Y.hdr.hist.descrip);
        scan_info_text{end+1} = sprintf('Slices: %d', Y.hdr.dime.dim(4));
        scan_info_text{end+1} = sprintf('Scan Resolution: %d, %d', Y.hdr.dime.dim(2), Y.hdr.dime.dim(3));
       
