function [Ydata H grad_table bval slice_order ]=readRECv3(filename, output_folder, Fname, P)
T=loadPARREC(filename);
H=T.hdr;
novols = H.NumberOfVolumes;
Ydata(:,:,:,1)=T.scans{1};
for k=2:novols
    Ytemp=T.scans{k};
    check=(size(Ytemp)==size(Ydata(:,:,:,1)));
    if(sum(check)~=3)  %grad file should be 3 columns
        disp('Not a dti');
        mkdir(output_folder);
        write_file = [output_folder filesep Fname '_' strrep(P,filesep,'-') ];
        fp = fopen([write_file '.txt'],'wt');
        fprintf(fp,'This file does not contain consistent volume sizes. Punting.');
        fclose(fp);
        delete(handles.figure1)
        clear handles f
        return;
    else
        Ydata(:,:,:,k)=Ytemp;
    end
end
        
   
slice_order=H.img(3).orient.slice_orientation;  %chose 3 incase Bo or other is in 1 or 2. max would be 7
if slice_order == 3
    Ydata=permute(Ydata,[1 3 2 4]);
    Ydata=flipdim(Ydata,3);Ydata=flipdim(Ydata,2); Ydata=flipdim(Ydata,1);
end
if slice_order == 1
    Ydata=Ydata(end:-1:1,end:-1:1,:,:);
end
if slice_order ==2
    Ydata=permute(Ydata,[3 1 2]);
end
for i=1:length(H.img), 
    GRAD(H.img(i).info.dynamic_scan_num,:) = H.img(i).special.diffusion_ap_fh_lr; 
    bval(H.img(i).info.dynamic_scan_num)=H.img(i).special.diffusion_b_factor;
end


table=GRAD;     



% gradients are always stored: ap, fh, rl   need to make  rl,ap,hf to match
% RAS of saved nii data
table(:,1)=GRAD(:,3); table(:,2)= GRAD(:,1); table(:,3)=-GRAD(:,2);
grad_table=table;
