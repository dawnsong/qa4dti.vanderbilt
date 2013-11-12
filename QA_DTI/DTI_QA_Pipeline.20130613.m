function [reportFile, state] = DTI_QA_Pipeline(file_name,output_folder,path_to_QA_DTI,n_bo,target)
%
% this is version 1.2
% last updated: June 13, 2013
% Past History
% 3- version 1.# release: June 13, 2013 - values saved in a text file for redcap
% 2- version 1.2 release: June 3, 2013 - reduced memory footprint
% 1- verson 1.1 released on nitrc: August 10, 2012
% 0- verson 1.0 released on nitrc: August 10, 2012
% Updates in version 1.1
% 1- added center of mass calculation
% 2- use center of mass to find mid coronal,axial,and sagital slice for
% page 3 and 4
% 3- page-3 bottom row of images now uses  Bo instead of the vector colormap
% 4- use robust BET options for FSL brain extraction (-R)
% 5- corrected bug. now normalizing chi-squared using the correct mask
% 6- added comments to clarify when using the statistcal mask, and when using the dti-mask
% 7- parse header to inlcude voxel mm resolution on page-1 of PDF
% 8- move n_bo (number of Bo averages) to user input option.
% 9- now store registration in output folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% TROUBLESHOOTING VALUES
% flag=1;
% bootnum=5;
% numsim = [1 3 3 3 3];

disp(['EBRL Version'])
% %%% NORMAL RUN VALUES
flag=18;
bootnum=1000;
numsim = [1 2000 4000 6000 8000];

if exist('n_bo','var')==0
    n_bo=1; % number of averaged Bo in data
end

%dawnsong added from runme
QA_pathname='/home/dawnsong/test/qa.dti/vanderbilt/Demo_DTI_QA_120807/QA_DTI';  % path to QA_DTI, end with QA_DTI, not QA_DTI/ 
%% ADD path to java and matlab 
java_path=sprintf('%s%smulti-atlas%smasi-fusion%sbin%s',QA_pathname,filesep,filesep,filesep,filesep);
addpath(genpath(QA_pathname))


QAmfiles_loc=sprintf('%s/QAmfiles/',QA_pathname);

%%flag = 0 for short run, else flag = time in hours for BOOT+SIMEX (should be 6+)

limchi=.2;
trble=sprintf('%s%sextra',output_folder,filesep);
tmp=sprintf('%s%stemp_folder',output_folder,filesep);
mkdir(trble)
mkdir(tmp)
testps=sprintf('%s%sQA_DTI.ps',output_folder,filesep);
%%%%%%%%%%%%%%%%%%%%%%%
%%%EXTRACT FILE PARTS
%%%%%%%%%%%%%%%%%%%%%
[P,Fname,E]=fileparts(file_name);
whos %% WHOS
name_fig1=sprintf('%s/homemade_mfiles/DTI_calp1.fig',QAmfiles_loc);
fig = openfig(name_fig1);
handles = guihandles(fig);
name_text={};
name_text{end+1}=sprintf('filename: %s',Fname);
name_text{end+1} = sprintf('folder: %s',P);
set(handles.text18,'String',name_text,'fontsize',9);
verStr = '1.0';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Prepatory Processing Depedning on if PARREC OR NII
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
niiname=file_name;
if isempty(E)==1
    cmmd=sprintf('%s.PAR',file_name);
    if exist(cmmd)~=0
        E='PAR';
    end
    cmmd=sprintf('%s.par',file_name);
    if exist(cmmd)~=0
        E='par';
    end
    cmmd=sprintf('%s.nii',file_name);
    
    if exist(cmmd)~=0
        E='nii';
        niiname=sprintf('%s.nii',file_name);
        if exist('target','var')~=0
            target=sprintf('%s.nii',target);
        end
    end
end

if E(1)=='P' | E(2)=='P' |E(1)=='p' | E(2)=='p'
    getParinfoV4
    if exist('target','var')
        target=loadPARREC(target);
        HT=target.hdr; targ=target.scans;
        if HT.img.orient.slice_orientation==1
            targ=targ(end:-1:1,end:-1:1,:);
        end
        if HT.img.orient.slice_orientation==2
            targ=permute(targ,[1 3 2]);
            targ=flipdim(targ,3); targ=flipdim(targ,2); targ=flipdim(targ,3);
        end
        
        resolutionT=HT.scn.fov(3)/size(Ydata,1);  %%RL divide by xdimension
        resolutionT(2)=HT.scn.fov(1)/size(Ydata,2); %%AP divide by ydimension
        resolutionT(3)=HT.scn.fov(2)/size(Ydata,3); %%FH divide by number of slices, slices MUST be FH!!
        nii=make_nii(targ, resolutionT,[0 0 0],16); save_nii(nii,'targ.nii')
    end
else
    
    [scan_info_text grad_file bval_vec resolution name_Y_data]= getNIIinfo(tmp,niiname);
    
    if exist('target','var')
        HT=load_nii(target); targ=HT.img;
        resolutionT=HT.hdr.dime.pixdim(2:4);
        nii=make_nii(targ, resolutionT,[0 0 0],16); save_nii(nii,'targ.nii')
    end
end
whos %% WHOS
set(handles.text2,'String',scan_info_text,'fontsize',7);

%grad)file is 3 columns, bval_vec is one row%  .................................................................................................................clear line

% load data from nii
Y=load_nii(name_Y_data);
Y_data=single(Y.img);
clear Y %_________________________________________________________________________clear line
pack
Nx = size(Y_data,1);
Ny = size(Y_data,2);
Nz = size(Y_data,3);
Nt = size(Y_data,4);
Ng = Nt-1; %number dwi weighted images
Reg_Im=zeros(size(Y_data),'single'); % last entry in Reg_im will be Bo
Reg_Im(:,:,:,end)=single(Y_data(:,:,:,end));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             DTI PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% REGISTRATION Too Target If Exists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bo_ref= make_nii(Y_data(:,:,:,end),resolution,[0 0 0], 16); bo_ref.hdr.dime.xyzt_units=2; name_bo_ref=sprintf('%s%sbo_ref.nii',tmp,filesep);
save_nii(bo_ref,name_bo_ref);
tranTARG=zeros(1,3,'single'); rotTARG=tranTARG;

if exist('target','var')
    !flirt -in bo_ref -ref targ -out temp -omat Taff -cost mutualinfo
    [ah avs]=system('avscale --allparams aff targ');
    [tranTARG rotTARG]=readAVS(avs);
    !gunzip temp.nii.gz
    load Taff;
    T1=load_nii('temp.nii'); Reg_Im(:,:,:,end)=single(T1.img);
    !rm temp*
end
whos %% WHOS
name_reg=sprintf('%s%s%s_registered.mat',tmp,filesep,Fname);
name_reg2=sprintf('%s%sRegistration_motion.mat',trble,filesep);
if exist(name_reg)~=0
    load(name_reg)
else
    register_c
    save(name_reg,'Reg_Im', 'affMat', 'rotation', 'translation')
    save(name_reg2,'affMat','rotation','translation');
end

clear Y_data dwi%...........................................................................................................................clear line
pack
nii=make_nii(Reg_Im,resolution,[0 0 0],16); nii.hdr.dime.xyzt_units=2; name_RegIm=sprintf('%s%sRegIm.nii',tmp,filesep); save_nii(nii,name_RegIm);
cmmd=sprintf('!cp %s %s',name_RegIm,trble);
eval(cmmd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PREPARE GRAD TABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grad_file=grad_file';
for vol=1:Ng
    F=inv(affMat(1:3,1:3,vol));
    FT=transpose(F);
    R=((F*FT)^-.5)*F;
    grad_file(1:3,vol)=R*grad_file(1:3,vol);
end
clear  R FT F%.........
pack

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MASK, strict one for stats loose one for DTI calc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_Mask2=sprintf('%s%sMask2',tmp,filesep); name_Mask=sprintf('%s%sMask',tmp,filesep);

cmmd=sprintf('!bet %s %s -f 0.2 -g 0 -m -R ',name_bo_ref, name_Mask); %for camino
eval(cmmd);
cmmd=sprintf('!bet %s %s -f 0.4 -g 0 -m -R ',name_bo_ref,name_Mask2);%for statistical analysis
eval(cmmd)

%clear any existing mask files so code doesnt need user input
cmmd=sprintf('!rm %s%sMas*.nii',tmp,filesep);
eval(cmmd)

%unzip mask files
cmmd=sprintf('!gunzip %s%s*.gz',tmp,filesep);
eval(cmmd);
whos %% WHOS
%find center of mass
tem_name=sprintf('%s_mask.nii',name_Mask);
temp=load_nii(tem_name); temp=temp.img;
%[centerX centerY centerZ]=com(temp); 
%dawnsong revise
[centerX centerY centerZ]=temp.hdr.dime(2:4)/2; 
centerX=round(centerX); centerY=round(centerY); centerZ=round(centerZ); %find better center of image
clear temp tem_name
pack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MASK ROIs and estimate sigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andrew Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file locations
atlases_file = '/scratch/mcr/Bo-atlas-final/atlases.txt';
labels_file = '/scratch/mcr/Bo-atlas-final/labels.txt';
targets = {name_bo_ref};
out_dir = sprintf('%s%smulti-atlas',output_folder,filesep);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% registration options
regtype = 'art';
regloc = '/scratch/mcr/art';
regopts = struct;

% fusion options
fusetype = 'NonLocalSpatialSTAPLE';
fuseloc = '/scratch/mcr/mipav/';
fuseopts = struct;
fuseopts.nc = [3 3 3 0];
%fuseopts.use_xvfb = true;

% general settings
runtype = 'single';

% read the atlas txt files
atlases = textread(atlases_file, '%s');
labels = textread(labels_file, '%s');

temp=tempname;
STATE_Atlas_Temp = [temp '.mat'];
save(STATE_Atlas_Temp);
clearvars -except atlases labels targets out_dir regtype runtype regloc regopts fuseloc  fusetype fuseopts STATE_Atlas_Temp
pack 

% register atlas and labels to target, save in out_dir
[out_atlas_dir out_labels_dir out_warp_dir] = ...
    run_registrations(atlases, labels, targets, out_dir, ...
    regtype, runtype, regloc, regopts);

% fuse atlases, segs is file location of final label
seg_dir = run_fusions(atlases, labels, targets, out_dir, out_labels_dir, ...
    out_atlas_dir, fuseloc, fusetype, runtype, fuseopts);
load(STATE_Atlas_Temp)
% get the location of the output file
segs{1} = sprintf('%s%smulti-atlas/fusion-%s/bo_ref_est.nii', ...
    output_folder, filesep, fusetype);
whos %% WHOS
cmmd=sprintf('!cp %s %s/multi_atlas_labels.nii', segs{1}, trble);
eval(cmmd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End Andrew Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y=load_untouch_nii(segs{1}); Y.hdr.dime.dim(6)=1; save_untouch_nii(Y,segs{1});
clear Y
pack
Y=load_nii(segs{1}); labelBo=Y.img; labelBo=labelBo(:);
maskROI=zeros(length(labelBo),25,'single'); ROI_sig=zeros(25,1,'single');
for roi=1:25
    maskROI(:,roi)=(labelBo==roi);
end
maskROI=(maskROI==1);
for gr=1:Ng
    RR=Reg_Im(:,:,:,gr); RR=RR(:);
    for roi=1:25
        ROI_sig(roi,gr)=std(RR(maskROI(:,roi)));
    end
end

m=median(ROI_sig,2); s=sort(m); sigmaEst=s(2);
cmmd=sprintf('save %s%ssigmaEst sigmaEst',trble,filesep);
eval(cmmd)

clear Y ROI_sig s RR labelBo%__________________________________________________________clearline
pack
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%CAMINO DTI FIT
% %%%%%%%%%%%%%%%%%%%%%%%%%%

%make scheme table
scheme=[grad_file' bval_vec'];
name_temp_text=sprintf('%s%stemp.txt',tmp,filesep);
save(name_temp_text,'scheme','-ascii');

name_temp2_text=sprintf('%s%stemp2.txt',tmp,filesep);
mim=fopen(name_temp2_text,'w+');
fprintf(mim,'VERSION:BVECTOR\n')
fclose(mim)
name_scheme_text=sprintf('%s%sscheme.txt',tmp,filesep);
cmmd=sprintf('!cat %s %s > %s',name_temp2_text,name_temp_text,name_scheme_text);
eval(cmmd)


cmmd=sprintf('!cp %s %s',name_scheme_text,trble);
eval(cmmd)
whos %% WHOS
%convert DWIdata to camino float
name_dwi_Bfloat=sprintf('%s%sdwi.Bfloat',tmp,filesep);

% cmmd=sprintf('!image2voxel -4dimage %s -outputfile %s',name_RegIm,name_dwi_Bfloat);
% eval(cmmd)

permRegIM = permute(Reg_Im,[4 1 2 3]);
fp = fopen(name_dwi_Bfloat,'wb','ieee-be');
fwrite(fp,permRegIM(:),'single');
fclose(fp);
clear permRegIM;

temp=tempname;
STATE_CAMINO_Temp = [temp '.mat'];
save(STATE_CAMINO_Temp)
clearvars -except name* sigmaEst trble tmp Ng STATE_CAMINO_Temp

%fit tensor with ROBUST fit
name_outliermap=sprintf('%s%soutliermap',tmp,filesep);
name_Mask_mask_nii=sprintf('%s%sMask_mask.nii',tmp,filesep);
name_dt_Bdouble=sprintf('%s%sdt.Bdouble',tmp,filesep);
cmmd=sprintf('!restore %s %s %d -bgmask %s -outliermap %s > %s',name_dwi_Bfloat, name_scheme_text, sigmaEst,name_Mask_mask_nii,name_outliermap,name_dt_Bdouble);
eval(cmmd)

%Convert tensor outputs to nii
name_tensor=sprintf('%s%stensor_',tmp,filesep);
cmmd=sprintf('!dt2nii -inputfile %s -inputdatatype double -header %s -outputroot %s',name_dt_Bdouble, name_RegIm,name_tensor);
eval(cmmd)

%cacluate FA and MD, output fa.nii and md.nii
name_fa=sprintf('%s%sfa',tmp,filesep); name_md=sprintf('%s%smd',tmp,filesep);
cmmd=sprintf('!cat %s | fa | voxel2image -outputroot %s -header %s',name_dt_Bdouble, name_fa, name_RegIm);
eval(cmmd)
cmmd=sprintf('!cat %s | md | voxel2image -outputroot %s -header %s',name_dt_Bdouble, name_md, name_RegIm);
eval(cmmd)

%calculate eigenvalues
name_eig_Bdouble=sprintf('%s%seig.Bdouble',tmp,filesep);
cmmd=sprintf('!dteig < %s > %s',name_dt_Bdouble,name_eig_Bdouble);
eval(cmmd)

name_V1_Bdouble=sprintf('%s%sV1.Bdouble',tmp,filesep);
cmmd=sprintf('!shredder 8 24 72 < %s > %s', name_dt_Bdouble, name_V1_Bdouble);
eval(cmmd)

name_dteig=sprintf('%s%sdteig',tmp,filesep);
cmmd=sprintf('!cat %s | dteig | voxel2image -components 12 -outputroot %s -header %s',name_dt_Bdouble, name_dteig, name_RegIm);
eval(cmmd)

%%%report outlier number
name_test=sprintf('%s%stest_',tmp,filesep);
cmmd=sprintf('!voxel2image -inputfile %s -outputroot %s -header %s -inputdatatype byte -components %d',name_outliermap,name_test,name_RegIm,Ng);
eval(cmmd)

load(STATE_CAMINO_Temp)

OutlierMask=ones(Nx,Ny,Nz,Ng,'int8');
for g=1:Ng
    name=sprintf('%s%04d.nii',name_test,g);
    T=load_nii(name); 
    OutlierMask(:,:,:,g)=int8(T.img);
end
nii=make_nii(OutlierMask,resolution,[0 0 0],2); nii.hdr.dime.xyzt_units=2;
nm=sprintf('%s%sOutlierMask.nii',trble,filesep);
whos %% WHOS
save_nii(nii,nm);
clear nii H Y_dwi Ydata art_home atlases atlases_file g gr scan_info_text scheme out_warp_dir out_labels_dir out_dir Onm
clear T %_______________________________________________________________________________clear line
pack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              STATISTICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PAGE 1%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %PLOT PATIENT MOTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rotation=rotation-repmat(rotTARG,Ng,1); translation=translation-repmat(tranTARG,Ng,1);
xt=translation(:,1); yt=translation(:,2); zt=translation(:,3);
xr=abs(rotation(:,1)); yr=abs(rotation(:,2)); zr=abs(rotation(:,3));


axes(handles.axes22);
plot(1:Ng,xt,'b-',1:Ng,yt,'c-',1:Ng,zt,'m-','linewidth',2);
hold
ylabel('Translation (mm)','fontweight','bold');
set(gca,'Yaxislocation','right','box','off','xlim',[1 Ng],'linewidth',2,'tickDir','out','fontweight','demi');
s = ['x';'y';'z'];
legend(s, 'Location','NorthOutside','orientation','horizontal')
axis tight
model=get(gca,'Position');
axes(handles.axes23);
plot(1:Ng,xr,'b-',1:Ng,yr,'c-',1:Ng,zr,'m-','linewidth',2);
hold
ylabel('Rotation (deg)','fontweight','bold');
set(gca,'Yaxislocation','right','box','off','xlim',[1 Ng],'linewidth',2,'TickDir','out','fontweight','demi')
xmarks=get(gca,'XTick');
fax=get(gca,'Position'); fax(1)=model(1); fax(3)=model(3); fax(4)=model(4);
set(gca,'Position',fax);
axis tight
clear affMat Taff translation rotation xt yt zt xr yr zr bo_ref%...................................................................................................clear line
pack
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%MASK IMAGES brain=stat_mask, dti_mask=mask for camino
%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_Mask2_mask_nii=sprintf('%s%sMask2_mask.nii',tmp,filesep);
stat_mask=load_nii(name_Mask2_mask_nii);stat_mask=stat_mask.img;stat_mask=stat_mask(:);
dti_mask = load_nii(name_Mask_mask_nii); dti_mask=dti_mask.img; dti_mask=dti_mask(:); dti_mask=(dti_mask==1);
brain=(dti_mask==1 & stat_mask==1);%%should be equal to stat_mask, but just in case.

whos %% WHOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULATE CHI_2 - uses dti-mask so all data is saved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_tensor_dt_nii=sprintf('%s%stensor_dt.nii',tmp,filesep);
T=load_nii(name_tensor_dt_nii); T=T.img;
name_md_nii=sprintf('%s%smd.nii',tmp,filesep); name_fa_nii=sprintf('%s%sfa.nii',tmp,filesep);
ADC=load_nii(name_md_nii); ADC=reshape(ADC.img,[],1);
FA=load_nii(name_fa_nii); FA=reshape(FA.img,[],1);

mADC=ADC(brain); mFA=FA(brain);
stdFA=std(FA(brain)); stdADC=std(ADC(brain));

[sd2 ModelData Errors]=calDTIrevC_norm(single(T),single(grad_file),single(bval_vec),Reg_Im,dti_mask);


clear T ADC FA mADC mFA stdFA stdADC %______________________________________________________-clear line
pack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %   MAKE MAP OF CHI-2 (use stat_mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
map=zeros(Ng,Nz,'single');

sd2(brain==0)=0; %use more restricted mask
brain_vol=reshape(brain,Nx,Ny,Nz);
sd2=reshape(sd2,Nx,Ny,Nz,Ng);

for volume=1:Ng
    for slice=1:Nz
        mask_slice=reshape(brain_vol(:,:,slice),1,[]);
        sss=reshape(sd2(:,:,slice,volume),1,[]);
        map(volume,slice)=sum(sss(mask_slice))./sum(mask_slice); %normalize by number of voxels in each image
        map(volume,slice)=map(volume,slice)*Ng;  % should be in units of chi_sq_p
    end
end
clear sss sdf k ks brain_vol %..................................................................................................................clear line
pack
savename=sprintf('%s/voxel_wise_chi_q',trble);
save(savename,'sd2','map')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %  PLOT Projections of MAP (map uses stat_mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % %
hite=trimmean(map,3); %%remove from ev
mite=trimmean(map',3)';

axes(handles.axes36) %%gradient number
plot(1:Ng,mite,'k','linewidth',2);
%ylabel('mean Chi Sq','fontweight','bold')
set(gca,'Yaxislocation','right','box','off','linewidth',2,'fontweight','demi')
fax=[0.4405    0.4990    0.4550    0.0769];
set(gca,'Position',fax,'TickDir','in','XTickLabel','');
axis tight
glim=get(gca,'ylim'); set(gca,'ylim',[0 glim(2)],'linewidth',2)
whos %% WHOS
axes(handles.axes29)  %%axial slice
plot(hite, 1:Nz,'k','linewidth',2);
set(gca,'XDir','reverse','Yaxislocation','right','Xaxislocation','bottom','TickDir','in','fontweight','demi','YTickLabel','')
rotateXLabels(gca,90);
set(gca,'box','off')
sliceAx=get(gca,'Position');
axis tight
set(gca,'xlim',[0 glim(2)],'linewidth',2,'XTick',[0 floor(10000*glim(2))/10000],'XTickLabel',{'0',floor(10000*glim(2))/10000});
pos2=[0.3360    0.0697    0.0969    0.4259];
set(gca,'Position',pos2);

clear hite mite %___________________________________________________________________________-clear line
pack
%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT RESTORE OUTLIERS--use stat_mask
%%%%%%%%%%%%%%%%%%%%%%%%
for g=1:Ng
    
    T=OutlierMask(:,:,:,g);
    T=T(:); T=T(brain);
    outs(g)=(sum(T)/sum(brain))*100;
end


axes(handles.axes40) %%gradient number


plot(1:Ng,outs,'k','linewidth',2);
set(gca,'Yaxislocation','right','box','off','fontweight','demi','TickDir','out','linewidth',2)
fax=get(gca,'Position'); fax(1)=model(1); fax(3)=model(3); fax(4)=model(4);
set(gca,'Position',fax);
xlabel('gradient #')
ylabel('% outliers','fontweight','demi')
axis tight
save(sprintf('%s/Outliers', trble), 'outs')
clear outs %____________________________________________________________________________________-clear line
pack
%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE HIST OF CHI_2--use stat_mask
%%%%%%%%%%%%%%%%%%%%%%%%

ChiSqHist


clear MM chisq chi_sq_p xx h%..........................................%%%..............................................................clear line
pack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PLOT MAP-EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num=5;
div=1/num;
st=1;

mx=max(max(max(max(Reg_Im)))); RI=Reg_Im/mx;
whos %% WHOS
% determine caxis from several values
figure(100)
delta=round(Nz/6); deltag=round(Ng/5);
RR=[];
for zc=1:5
    for gc=1:4
        RR=[RR RI(:,:,delta*zc,deltag*gc)];
    end
end
imagesc(RR); ct=get(gca,'clim'); ct(2)=ct(2)/2;
close(100)
clear RR M mm a d
pack
for b=1:num
    t=1:num; t=flipdim(t,2);
    name=['handles.axes' num2str(t(b)*2-1)];
    en=round(div*b*Nz); %determine end value for sampling M
    M=map(:,st:en);
    mm=min(reshape(M,[],1));
    [a d]=find(M==mm);
    axes(eval(name))
    RR = squeeze(RI(:,:,st+d(1)-1,a(1)));
    RR = permute(RR, [2 1]); RR=flipdim(RR,1); RR=repmat(RR,[1 1 3]);
    RRn=imadjust(RR, ct, [0 1]);
    imagesc(RRn);
    
    if b==num
        text(Nx*.02,Ny*.1,'(grd,slice)','fontsize',9,'color',[1 0 0],'fontweight','bold')
        text(Nx*.66,Ny*.1,'Chi-Sq','fontsize',9,'color',[1 0 0])
    end
    set(gca,'visible','off','clim',ct);
    mm=round(mm*1000)/1000;
    if mm>1
        value='> 1.0';
    else
        value=num2str(mm);
    end
    text(Nx*.7,Ny*.94,value,'fontsize',9,'color',[1 0 0])
    location=['(' num2str(a(1)) ',' num2str(st+d(1)-1) ')'];
    text(Nx*.03,Ny*.94,location,'fontsize',9,'color',[1 0 0])
    mm=max(reshape(M,[],1));
    [f g]=find(M==mm);
    name2=['handles.axes' num2str( t(b)*2 )];
    axes(eval(name2))
    RR2 = RI(:,:,st+g(1)-1,f(1));
    RR2 = permute(RR2, [2 1]); RR2=flipdim(RR2,1);RR2=repmat(RR2,[1 1 3]);
    RR2n=imadjust(RR2,ct,[0 1]);
    imagesc(RR2n); %colormap(gray);
    set(gca,'clim',ct);
    set(gca,'visible','off')
    location=['(' num2str(f(1)) ',' num2str(st+g(1)-1) ')'];
    mm=round(mm*1000)/1000;
    if mm>1
        value='> 1.0';
    else
        value=num2str(mm);
    end
    text(Nx*.03,Ny*.94,location,'fontsize',9,'color',[1 0 0])
    text(Nx*.7,Ny*.94,value,'fontsize',9,'color',[1 0 0])
    st=en+1;
end
whos %% WHOS
clear RR RI RR2 RRn mm value st en %____________________________________________________________________clear line
pack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PLOT MAP  -- uses stat_mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.axesM)
imagesc(map');
colorbar('location','WestOutside');
load QAmfiles/homemade_mfiles/colmap
set(gcf,'colormap',colmap)
caxis([0 limchi]);
cbar=findobj('tag','Colorbar');
fax=get(gca,'Position'); fax(1)=model(1); fax(3)=model(3); fax(2)=sliceAx(2); fax(4)=sliceAx(4);
cpos=[0.30400    0.0682    0.0265    0.4000];
set(cbar,'Position',cpos);
set(gca,'Position',fax);
set(cbar,'YtickMode','Manual','YTick',[0 .05 .1 .15 .2],'linewidth',2)
set(gca,'YDir','normal')
set(gca,'Xaxislocation','bottom','Yaxislocation','right','linewidth',2,'TickDir','out','box','off','fontweight','bold')
xlabel('Gradient #','fontsize',14,'fontweight','bold')
ylabel('Axial Slice #','fontsize', 14,'fontweight','bold')
clear map Reg_Im M%...............................................................................................................clear line
pack
% F=getframe(fig);
% close(fig)
% figure(1)
% imagesc(F.cdata);
% colormap(F.colormap);
set(gcf,'PaperPosition',[.25 .15 8.25 10.7]);
print('-dpsc2',testps,gcf);
close(gcf)
% clear F


clear stat_mask brain  yh ym xxx slice slice_order slice_Ax mask_slice m labels RR2n labels_file hhh h ho csqp colmap
pack


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PAGE 2    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
whos %% WHOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FA and MD boxplots --uses mask from multi-atlas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FA=load_nii(name_fa_nii); MD=load_nii(name_md_nii);
FA=FA.img; MD=MD.img;
name_fig2=sprintf('%s/homemade_mfiles/DTI_calp2.fig',QAmfiles_loc);
fig2 = openfig(name_fig2);
handles2=guihandles(fig2);
page2FAMD_KKI_array %_____________________________________________clear lines within code


clear grp
pack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SHRINK DATA FOR SIMEX AND BOOT --
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thours=flag;
tsec=thours*60*60;
voxNum=tsec/(2.2); % ~ 2.2 seconds per voxel second with Boot (at 5000reps)
NewMask = downSamp_full(maskROI,voxNum);


cmmd=sprintf('save(''%s%sNewMask.mat'',''NewMask'')',trble,filesep);
eval(cmmd)
NewMask=sum(NewMask,2);
NewMask=(NewMask==1);
percentage=sum(NewMask)/sum(sum(maskROI)); percentage=percentage*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SIMEX not with RESTORE (yet?)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=load_nii(name_RegIm);R=R.img;
R=reshape(R,Nx*Ny*Nz,Ng+1);
sampleVox=R(NewMask,:);
sampleVox=permute(sampleVox,[2 1]);
clear R %_____________________________________________________________________clear line
[FAsmx Bias]=simex(sampleVox,bval_vec,grad_file(:,1:end-1),sigmaEst,n_bo, numsim);
clear FAsmx %____________________________________________________________________-clear line
clear sampleVox
pack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BOOT (not with RESTORE yet?)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ModelData=reshape(ModelData,Ng,Nx*Ny*Nz);
sampleVox2=ModelData(:,NewMask); Errors=Errors(:,NewMask);
FAsample=FA(NewMask);
FAboot=boot(sampleVox2, Errors, bval_vec,grad_file(:,1:end-1),bootnum,FAsample');
ffstd=nanstd(FAboot);
nm=sprintf('%s/ModelData',trble);
save(nm,'ModelData')
whos %% WHOS
clear FAboot ModelData Errors;
pack
page2SimAndBoot_KKI_array %______________________________________________________________clear lines within code


set(gcf,'PaperPosition',[.25 .15 8.25 10.7]);
print('-dpsc2','-append',testps,gcf);
close(gcf)
% clear F

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PAGE 3   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_fig3 = sprintf('%s/homemade_mfiles/DTI_calp3.fig',QAmfiles_loc);
fig3 = openfig(name_fig3);
handles3 = guihandles(fig3);
page3pow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Show Bo DIFFUSION DIRCTION MAP--uses dti_mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bzero=load_nii(name_bo_ref); Bzero=Bzero.img;
Bzero=Bzero(:); Bzero(dti_mask==0)=0;Bzero=reshape(Bzero,Nx,Ny,Nz);
axial=Bzero(:,:,centerZ);
coronal=squeeze(Bzero(:,centerY,:));
sag=squeeze(Bzero(centerX,:,:));
findclim=[coronal sag axial]; figure(100); imagesc(findclim); ca=get(gca,'clim'); close(100); ca(2)=ca(2)/1.5;
axes(handles3.axes4(3))
imagesc(coronal')
hold on
set(gca,'XTick',[],'YTick',[],'YDir','normal','clim',ca)
xlabel('INFERIOR','fontweight','demi','fontsize',12)
ylabel('LEFT','fontweight','demi','fontsize',12)
labelTitle = sprintf('Coronal Slice, #%02d',round(Nx/2));
title(labelTitle,'fontweight','demi','fontsize',12)
pos=get(gca,'Position');


%Plot axial slice
axes(handles3.axes4(2))
imagesc(axial')
hold on
set(gca,'XTick',[],'YTick',[],'YDir','normal','clim',ca)
xlabel('POSTERIOR','fontweight','demi','fontsize',12)
%ylabel('LEFT','fontweight','demi','fontsize',12)
labelTitle = sprintf('Axial Slice, #%02d',round(Nz/2));
title(labelTitle,'fontweight','demi','fontsize',12)
p2=get(gca,'Position'); p2(2:end)=pos(2:end); set(gca,'Position',p2);
whos %% WHOS
% Sagittal Slice
axes(handles3.axes4(1))
imagesc(sag')
hold on

xlabel('INFERIOR','fontweight','demi','fontsize',12)
ylabel('ANTERIOR','fontweight','demi','fontsize',12);
labelTitle = sprintf('Sagittal Slice, #%02d',round(Ny/2));
title(labelTitle,'fontweight','demi','fontsize',12)
p2=get(gca,'Position'); p2(2:end)=pos(2:end); set(gca,'Position',p2);
set(gca,'XTick',[],'YTick',[],'Yaxislocation','right','YDir','normal','clim',ca)


set(gcf,'PaperPosition',[.2 .15 8.2 10.7]);
print('-dpsc2','-append',testps,gcf);
close(gcf)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PAGE 4  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name_fig4 = sprintf('%s/homemade_mfiles/DTI_calp4.fig',QAmfiles_loc);
fig4 = openfig(name_fig4);
handles4 = guihandles(fig4);


V1=ones(Nx,Ny,Nz,3);
for k=2:4
    dteig=sprintf('%s000%d.nii',name_dteig,k);
    temp=load_nii(dteig);
    V1(:,:,:,k-1)=temp.img;
end
name_v1_nii=sprintf('%s%sV1.nii',tmp,filesep);
nii=make_nii(V1,resolution,[0 0 0],16); save_nii(nii,name_v1_nii); nii.hdr.dime.xyzt_units=2;
clear V1
V1=load_nii(name_v1_nii); av=V1.img; fa=load_nii(name_fa_nii);
f=fa.img; f(f>1)=1; f(f<0)=0;
faV=repmat(f,[1 1 1 3]);
dtmap=faV.*abs(av);
dtmap(dtmap==0)=1;
Nx=size(av,1); Ny=size(av,2); Nz=size(av,3);

name_Mask_mask_nii=sprintf('%s%sMask_mask.nii',tmp,filesep);
mask=load_nii(name_Mask_mask_nii);mask=mask.img;

makeVectorMap(centerZ, centerY,dtmap,av,handles4.axes53,handles4.axes56,mask,tmp,1)
makeVectorMap(centerZ, centerY,dtmap,av,handles4.axes52,handles4.axes55,mask,tmp,2)
makeVectorMap(centerZ, centerY,dtmap,av,handles4.axes49,handles4.axes54,mask,tmp,3)
whos %% WHOS
clear av dtmap
pack
set(gcf,'PaperPosition',[.2 .15 8.2 10.7]);
print('-dpsc2','-append',testps,gcf);
close(gcf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  PUT OUTPUTS IN FOLDER %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name=[output_folder filesep 'QA_maps'];
mkdir(name)

%save good stuff
cmmd=sprintf('!cp %s %s %s %s %s %s %s%sQA_maps', name_Mask_mask_nii, name_Mask2_mask_nii, name_dt_Bdouble, name_md_nii, name_fa_nii, name_v1_nii,output_folder,filesep);
eval(cmmd)

%save QA report
% write_file = [output_folder filesep Fname '_' strrep(P,filesep,'-') '.pdf' ];
% ps2pdf('psfile',testps,'pdffile',write_file);

cmmd=sprintf('!cp %s%spowerData.mat %s',tmp,filesep,trble);
eval(cmmd)

close all
%
cmmd=sprintf('!rm -r %s',tmp);
eval(cmmd)

whos %% WHOS

