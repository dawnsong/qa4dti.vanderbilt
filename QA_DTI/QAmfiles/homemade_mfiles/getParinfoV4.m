
[P,Fname,E]=fileparts(file_name);
[Ydata H grad_f bval slice_order ]=readRECv3(file_name, output_folder, Fname,P);
 %%for correct RL AP FH directions !!! 

% and Bo, and check for appropiate slice and gradient order
if(size(grad_f,1)<7)
    disp('Not a dti');
    mkdir(output_folder);
    write_file = [output_folder filesep Fname '_' strrep(P,filesep,'-') ];
    fp = fopen([write_file '.txt'],'wt');
    fprintf(fp,'This files does not appear to be a DTI file. Punting.');
    fclose(fp);
        delete(handles.figure1)
    clear handles f
    return;
end

        %%%%%%%%%%%%%%Scan information
        scan_info_text = {};
        scan_info_text{end+1} = sprintf('Exam Name:  %s', Fname );
        scan_info_text{end+1} = sprintf('Protocol: %s', H.scn.protocol_name );
        scan_info_text{end+1} = sprintf('Date/Time: %s', H.info.exam_datetime);
        scan_info_text{end+1} = sprintf('Maximum Slices: %d', H.max.num_slices);
        scan_info_text{end+1} = sprintf('Maximum Dynamics: %d', H.max.num_dynamics);
        scan_info_text{end+1} = sprintf('Prepare Direction: %s', H.orient.prep_dir);
        scan_info_text{end+1} = sprintf('Technique: %s', H.scn.technique);
        scan_info_text{end+1} = sprintf('Scan Duration: %.1f', H.scn.scan_dur );
        scan_info_text{end+1} = sprintf('TR [msec]: %.1f', H.scn.rep_time );
        scan_info_text{end+1} = sprintf('Scan Resolution: %d, %d', H.scn.scan_res);
        scan_info_text{end+1} = sprintf('FOV (AP, FH, RL) [mm]: %.1f, %.1f, %.1f', H.scn.fov);
        scan_info_text{end+1} = sprintf('Water Fat Shift: %.3f', H.scn.water_fat_shift);
        scan_info_text{end+1} = sprintf('EPI Factor: %.1f', H.special.epi_factor);
        scan_info_text{end+1} = sprintf('Par/Rec Version: %.2f', H.version);
        
        resolution=H.scn.fov(3)/size(Ydata,1);  %%RL divide by xdimension
        resolution(2)=H.scn.fov(1)/size(Ydata,2); %%AP divide by ydimension
        resolution(3)=H.scn.fov(2)/size(Ydata,3); %%FH divide by number of slices, slices MUST be FH!!
   
% % load only dwi associated slices and find Bo slice
Ng = size(grad_f,1);
grad_file = [];
bval_vec = [];
Y_dwi=[];
count = 1;
for vol = 1:Ng
    if norm(grad_f(vol,:))<1.001 && norm(grad_f(vol,:))>.998
        Y_dwi(:,:,:,count) = Ydata(:,:,:,vol);
        grad_file(count,:) = grad_f(vol,:);
        bval_vec(count)=bval(vol);
        count = count+1;
    end
    if norm(grad_f(vol,:)) == 0 && bval(vol)==0
        Y_Bo = Ydata(:,:,:,vol);
    end
end
if exist('Y_Bo')~=1
    mkdir(output_folder);
    write_file = [output_folder filesep Fname '_' strrep(P,filesep,'-') ];
    fp = fopen([write_file '.txt'],'wt');
    fprintf(fp,'This files does not contain a Bo image. Punting.');
    fclose(fp);
        delete(handles.figure1)
    clear handles f
    return;
end
Y_data = Y_dwi;
Y_data(:,:,:,end+1) = Y_Bo;
bval_vec(end+1)=0;
grad_file(end+1,:) = [0 0 0];
clear Y_Bo bval grad_f%..................................................................................................................clear line
    
    name_Y_data=sprintf('%s%sY_data.nii',tmp,filesep);
    nii=make_nii(Y_data, resolution,[0 0 0],16); nii.hdr.dime.xyzt_units=2; save_nii(nii,name_Y_data)

