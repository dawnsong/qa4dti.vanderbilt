function [reportFile, state] = DTI_QA_Pipeline(file_name,output_folder,path_to_QA_DTI)
%note for author: function from DTI_QA_Report_V17, 120730_QA_DTI.tar
%this is release version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TROUBLESHOOTING VALUES
% flag=1;
% bootnum=5;
% numsim = [1 3 3 3 3];


% %%% NORMAL RUN VALUES
flag=18; % parameter controlling voxel sub-sampling for Bootstrap and SIMEX
bootnum=1000; % bootstrap monte-carlo simulation numbers
numsim = [1 2000 4000 6000 8000]; %SIMEX monte-carlo simulation numbers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%n_bo=5;  %% number of averaged Bo images in Bo.
n_bo=1;  %% dawnsong, 20140105
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% BEGIN MAIN PROGRAM %%%%%%% 

progress='Loading data. Converting to RAS nifti'
QAmfiles_loc=sprintf('%s%sQAmfiles%s',path_to_QA_DTI,filesep,filesep);

%%flag = 0 for short run, else flag = time in hours for BOOT+SIMEX (should be 6+)

limchi=.2;
trble=sprintf('%s%sextra',output_folder,filesep);
tmp=sprintf('%s%stemp_folder',output_folder,filesep);
mkdir(trble)
mkdir(tmp)
testps=sprintf('%s%stemp.ps',tmp,filesep);
%%%%%%%%%%%%%%%%%%%%%%%
%%%EXTRACT FILE PARTS
%%%%%%%%%%%%%%%%%%%%%
[P,Fname,E]=fileparts(file_name);
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

    set(handles.text2,'String',scan_info_text,'fontsize',9);

   %grad)file is 3 columns, bval_vec is one row%  .................................................................................................................clear line

% load data from nii
Y=load_nii(name_Y_data);
Y_data=Y.img;
Nx = size(Y_data,1);
Ny = size(Y_data,2);
Nz = size(Y_data,3);
Nt = size(Y_data,4);
Ng = Nt-1; %number dwi weighted images
Reg_ImA=zeros(size(Y_data)); % last entry in Reg_im will be Bo
Reg_ImA(:,:,:,end)=Y_data(:,:,:,end);
clear Y %_________________________________________________________________________clear line

progress='Finished loading data. Registering DWI to Bo'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             DTI PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% REGISTRATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bo_ref= make_nii(Y_data(:,:,:,end),resolution,[0 0 0], 16); bo_ref.hdr.dime.xyzt_units=2; name_bo_ref=sprintf('%s%sbo_ref.nii',tmp,filesep);
save_nii(bo_ref,name_bo_ref);

name_reg=sprintf('%s%s%s_registered.mat',tmp,filesep,Fname);
if exist(name_reg)~=0
        load(name_reg)
    else
        register_c
        save(name_reg,'Reg_Im', 'affMat', 'rotation', 'translation')
end


clear Y_data dwi Reg_ImA%...........................................................................................................................clear line
nii=make_nii(Reg_Im,resolution,[0 0 0],16); nii.hdr.dime.xyzt_units=2; name_RegIm=sprintf('%s%sRegIm.nii',tmp,filesep); save_nii(nii,name_RegIm);
cmmd=sprintf('!cp %s %s',name_RegIm,trble);
eval(cmmd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PREPARE GRAD TABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progress = 'Finished registration. Rotating Gradient table'
grad_file=grad_file';
for vol=1:Ng
    F=inv(affMat(1:3,1:3,vol));
    FT=transpose(F);
    R=((F*FT)^-.5)*F;
    grad_file(1:3,vol)=R*grad_file(1:3,vol);
end
    clear  R FT F%.........

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MASK, strict one for stats loose one for DTI calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progress='Finished gradient table rotation. Masking bo'
name_Mask=sprintf('%s%sMask',tmp,filesep); name_Mask2=sprintf('%s%sMask2',tmp,filesep);

cmmd=sprintf('!bet %s %s -f 0.2 -g 0 -m',name_bo_ref,name_Mask); %for use with camino
eval(cmmd)
cmmd=sprintf('!bet %s %s -f 0.4 -g 0 -m',name_bo_ref,name_Mask2);%for statistical analysis 
eval(cmmd) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MASK ROIs and estimate sigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progress='Finished masking. Performing multi-atlas segmentation and noise estimation'
% ANDREW CODEsettings
art_home = sprintf( '%s/../multi-atlas/art', QAmfiles_loc );
nlsloc = sprintf('%s/../multi-atlas/nls', QAmfiles_loc);   
atlases_file = sprintf('%s/../multi-atlas/atlas_dir/atlases.txt', QAmfiles_loc);
labels_file = sprintf('%s/../multi-atlas/atlas_dir/labels.txt', QAmfiles_loc    );
%art_home = '/scratch/mcr/bin/art';  %%%%%%%%%% andrew codeSERVER
%nlsloc = '/scratch/mcr/bin/nls';    %%%%%%%%%%SERVER
%atlases_file = '/scratch/mcr/atlas_dir/atlases.txt'; %%% location of atlases
%labels_file = '/scratch/mcr/atlas_dir/labels.txt'; %%% location of labels
targets = {name_bo_ref};
out_dir = sprintf('%s%smulti-atlas',output_folder,filesep);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

atlases = textread(atlases_file, '%s');
labels = textread(labels_file, '%s');
atlases  =regexprep(atlases,  '^(.*)',sprintf('%s/../multi-atlas/atlas_dir/$1', QAmfiles_loc));
labels =regexprep(labels, '^(.*)',sprintf('%s/../multi-atlas/atlas_dir/$1',     QAmfiles_loc));

[out_atlas_dir out_labels_dir out_warp_dir] = run_art_registrations_single(atlases, labels, targets, out_dir, art_home); % register atlas and labels to target, save in out_dir
segs = run_nls_fusions_single(atlases, labels, targets, out_dir, out_atlas_dir, out_labels_dir, nlsloc); % fuses atlases, segs is file location of final label
cmmd=sprintf('!cp %s%smulti-atlas/nls-fusion/bo_ref_est.nii %s/multi_atlas_labels.nii',output_folder,filesep,trble);
eval(cmmd)

% % % %   END andrew code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y=load_untouch_nii(segs{1}); Y.hdr.dime.dim(6)=1; save_untouch_nii(Y,segs{1});
clear Y
Y=load_nii(segs{1}); labelBo=Y.img; labelBo=labelBo(:); 
maskROI=zeros(length(labelBo),25); ROI_sig=zeros(25,1);
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%CAMINO DTI FIT
% %%%%%%%%%%%%%%%%%%%%%%%%%%
progress='Finished segmentation. Running CAMINO-RESTORE diffusion tensor fit'
%disable since CAMINO NEEDED!!!
%clear any existing mask files so code doesnt need user input
    cmmd=sprintf('!rm %s%sMas*.nii',tmp,filesep);
    %eval(cmmd)%dawnsong
    
%unzip mask files
cmmd=sprintf('!gunzip -f %s%s*.gz',tmp,filesep);
eval(cmmd);

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
   
    %convert DWIdata to camino float
    name_dwi_Bfloat=sprintf('%s%sdwi.Bfloat',tmp,filesep);
    cmmd=sprintf('!image2voxel -4dimage %s -outputfile %s',name_RegIm,name_dwi_Bfloat);
    eval(cmmd)

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

    %calculate FA and MD, output fa.nii and md.nii
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

dti_mask = load_nii(name_Mask_mask_nii); dti_mask=dti_mask.img; dti_mask=dti_mask(:);tot=sum(dti_mask);
OutlierMask=ones(Nx,Ny,Nz,Ng);
for g=1:Ng
    name=sprintf('%s%04d.nii',name_test,g);
    T=load_nii(name); T=T.img;
    OutlierMask(:,:,:,Ng)=T;
    outs(g)=(sum(T(:))/tot)*100;
end
nii=make_nii(OutlierMask,resolution,[0 0 0],2); nii.hdr.dime.xyzt_units=2; 
nm=sprintf('%s%sOutlierMask.nii',trble,filesep);
save_nii(nii,nm);
clear nii H Y_dwi Ydata art_home atlases atlases_file g gr scan_info_text scheme out_warp_dir out_labels_dir out_dir OutlierMask nm 
clear T %_______________________________________________________________________________clear line
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
progress='Finished calculating DTI. Preparaing page-1 of 4 for QA report'

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
clear affMat Taff translation rotation xt yt zt xr yr zr%...................................................................................................clear line

%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%MASK IMAGE STATISTICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_Mask2_mask_nii=sprintf('%s%sMask2_mask.nii',tmp,filesep);
stat_mask=load_nii(name_Mask2_mask_nii);stat_mask=stat_mask.img;stat_mask=stat_mask(:);

clear bo_ref%...................................................................................................clear line
brain=(dti_mask==1 & stat_mask==1);
nobrain=(brain==0);
Mask_data_vec = reshape(Reg_Im,[],size(Reg_Im,4));   %each column represents entire volume
for j=1:Nt;
    Mask_data_vec(nobrain,j)=NaN;
end
Mask_dwi_vec=Mask_data_vec(:,1:end-1);
clear   M2 Mask_data_vec gives%...................................................................................................clear line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULATE CHI_2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progress='Calculating pixel-chi-squared.'
name_tensor_dt_nii=sprintf('%s%stensor_dt.nii',tmp,filesep);
T=load_nii(name_tensor_dt_nii); T=T.img;
name_md_nii=sprintf('%s%smd.nii',tmp,filesep); name_fa_nii=sprintf('%s%sfa.nii',tmp,filesep);
ADC=load_nii(name_md_nii); ADC=reshape(ADC.img,[],1);
FA=load_nii(name_fa_nii); FA=reshape(FA.img,[],1);

mADC=ADC(brain); mFA=FA(brain);
stdFA=std(FA(brain)); stdADC=std(ADC(brain));

[sd ModelData Errors]=calDTIrevC_norm(T,grad_file,bval_vec,Reg_Im,dti_mask);

clear T ADC FA mADC mFA stdFA stdADC dti_mask%______________________________________________________-clear line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %   MAKE MAP OF CHI-2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
map=zeros(Ng,Nz); 
brain_vol=reshape(brain,Nx,Ny,Nz);
sd=reshape(sd,Nx,Ny,Nz,Ng);
for volume=1:Ng
    for slice=1:Nz
        mask_slice=reshape(brain_vol(:,:,slice),1,[]);
        sss=reshape(sd(:,:,slice,volume),1,[]);
        map(volume,slice)=sum(sss(mask_slice))./sum(mask_slice); %normalize by number of voxels in each image
        map(volume,slice)=map(volume,slice)*Ng;  % should be in units of chi_sq_p
    end
end
clear sss sdf k ks brain_vol %..................................................................................................................clear line
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %  PLOT Projections of MAP
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

%%%%%%%%%%%%%%%%%%%%%%%%%%
 % PLOT RESTORE OUTLIERS
 %%%%%%%%%%%%%%%%%%%%%%%%
 axes(handles.axes40) %%gradient number

 plot(1:Ng,outs,'k','linewidth',2);
    set(gca,'Yaxislocation','right','box','off','fontweight','demi','TickDir','out','linewidth',2)
    fax=get(gca,'Position'); fax(1)=model(1); fax(3)=model(3); fax(4)=model(4);
    set(gca,'Position',fax);
    xlabel('gradient #')
    ylabel('% outliers','fontweight','demi')
    axis tight
   
  clear outs %____________________________________________________________________________________-clear line
 %%%%%%%%%%%%%%%%%%%%%%%%%%
 % MAKE HIST OF CHI_2
 %%%%%%%%%%%%%%%%%%%%%%%%
 
ChiSqHist
savename=sprintf('%s/voxel_wise_chi_q',trble);
save(savename,'sd','map')

clear MM chisq chi_sq_p sd xx h%..........................................%%%..............................................................clear line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % PLOT MAP-EXAMPLES  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
num=5;
div=1/num;
st=1;

mx=max(max(max(max(Reg_Im)))); RI=Reg_Im/mx;
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
for b=1:num
        t=1:num; t=flipdim(t,2);
        name=['handles.axes' num2str(t(b)*2-1)];
        en=round(div*b*Nz); %determine end value for sampling M
        M=map(:,st:en);
        mm=min(reshape(M,[],1));
        [a d]=find(M==mm);
        if isempty(d), continue, end
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
clear RR RI RR2 RRn mm value st en %____________________________________________________________________clear line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % PLOT MAP  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
map(isnan(map))=0.2; %  fill NaN with max value before drawing
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

% F=getframe(fig);
% close(fig)
% figure(1)
% imagesc(F.cdata);
% colormap(F.colormap);
set(gcf,'PaperPosition',[.25 .15 8.25 10.7]);
print('-dpsc2',testps,gcf);
close(gcf)
% clear F


clear stat_mask brain Mask_dwi_vec yh ym xxx slice slice_order slice_Ax mask_slice m labels RR2n labels_file hhh h ho csqp colmap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PAGE 2    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progress='Finished page 1. Preparing page-2 of 4 for QA report'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FA and MD boxplots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FA=load_nii(name_fa_nii); MD=load_nii(name_md_nii);
FA=FA.img; MD=MD.img;
name_fig2=sprintf('%s/homemade_mfiles/DTI_calp2.fig',QAmfiles_loc);
fig2 = openfig(name_fig2);
handles2=guihandles(fig2);
page2FAMD_KKI_array %_____________________________________________clear lines within code


clear grp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SHRINK DATA FOR SIMEX AND BOOT
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
progress='Running SIMEX'
R=load_nii(name_RegIm);R=R.img;
R=reshape(R,Nx*Ny*Nz,Ng+1);
sampleVox=R(NewMask,:);
sampleVox=permute(sampleVox,[2 1]);
clear R %_____________________________________________________________________clear line
fsimex=sprintf('%s/simex.mat',trble);
if exist(fsimex),
    load(fsimex)
else
    [FAsmx Bias]=simex(sampleVox,bval_vec,grad_file(:,1:end-1),sigmaEst,n_bo, numsim);
    clear FAsmx %____________________________________________________________________-clear line
    save(fsimex, 'Bias')
end
clear sampleVox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BOOT (not with RESTORE yet?)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progress='SIMEX finished. Running Bootstrap'
 ModelData=reshape(ModelData,Ng,Nx*Ny*Nz); 
sampleVox2=ModelData(:,NewMask); Errors=Errors(:,NewMask);
FAsample=FA(NewMask);
fboot=sprintf('%s/boot.mat',trble);
if exist(fboot),
    load(fboot)
else
    FAboot=boot(sampleVox2, Errors, bval_vec,grad_file(:,1:end-1),bootnum,FAsample');
    save(fboot, 'FAboot')
end
ffstd=nanstd(FAboot);
nm=sprintf('%s/ModelData',trble);
save(nm,'ModelData')
clear FAboot ModelData Errors;

page2SimAndBoot_KKI_array %______________________________________________________________clear lines within code


set(gcf,'PaperPosition',[.25 .15 8.25 10.7]);
print('-dpsc2','-append',testps,gcf);
close(gcf)
% clear F

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PAGE 3   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progress='Finished page-2. Preparing page-3 of 4 for QA report.'
name_fig3 = sprintf('%s/homemade_mfiles/DTI_calp3.fig',QAmfiles_loc);   
fig3 = openfig(name_fig3);
handles3 = guihandles(fig3);
page3pow
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALC DIFFUSION DIRCTION MAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V1=ones(Nx,Ny,Nz,3);
    for k=2:4
        dteig=sprintf('%s000%d.nii',name_dteig,k);
        temp=load_nii(dteig);
        V1(:,:,:,k-1)=temp.img;
    end
    name_v1_nii=sprintf('%s%sV1.nii',tmp,filesep);
    nii=make_nii(V1,resolution,[0 0 0],16); save_nii(nii,name_v1_nii); nii.hdr.dime.xyzt_units=2; 
    V1=abs(V1); 
    FA=load_nii(name_fa_nii); FA=FA.img;
    FA(FA>1)=1; FA(FA<0)=0;
    cc=V1.*repmat(FA,[1 1 1 3]); cc=cc(:); nobrain_b=([nobrain nobrain nobrain]==1); 
    cc(nobrain_b)=1; cc=reshape(cc,Nx,Ny,Nz,3);
    cc=flipdim(cc,2);
    cc=permute(cc,[2 1 3 4]);
    clear FA V1 nobrain_b dteig temp%____________________________________________clear line
    page3checkV1_v2 %______________________________________________________________clear lines in code
    clear cc %___________________________________________________________________clear line

set(gcf,'PaperPosition',[.2 .15 8.2 10.7]);
print('-dpsc2','-append',testps,gcf);
close(gcf)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PAGE 4  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progress='Finished page-3. Preparing page-4 of 4 for QA report.'
name_fig4 = sprintf('%s/homemade_mfiles/DTI_calp4.fig',QAmfiles_loc);   
fig4 = openfig(name_fig4);
handles4 = guihandles(fig4);

V1=load_nii(name_v1_nii); av=V1.img; fa=load_nii(name_fa_nii);
f=fa.img; f(f>1)=1; f(f<0)=0; 
faV=repmat(f,[1 1 1 3]);
dtmap=faV.*abs(av);
dtmap(dtmap==0)=1;
Nx=size(av,1); Ny=size(av,2); Nz=size(av,3);

name_Mask_mask_nii=sprintf('%s%sMask_mask.nii',tmp,filesep);
mask=load_nii(name_Mask_mask_nii);mask=mask.img;

makeVectorMap(round(Nz/2), round(Ny/2),dtmap,av,handles4.axes53,handles4.axes56,mask,tmp,1)
makeVectorMap(round(Nz/2), round(Ny/2),dtmap,av,handles4.axes52,handles4.axes55,mask,tmp,2)
makeVectorMap(round(Nz/2), round(Ny/2),dtmap,av,handles4.axes49,handles4.axes54,mask,tmp,3)

clear av dtmap

set(gcf,'PaperPosition',[.2 .15 8.2 10.7]);
print('-dpsc2','-append',testps,gcf);
close(gcf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  PUT OUTPUTS IN FOLDER %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
progress='Finished report. Storing Outputs.'
name=[output_folder filesep 'QA_maps'];
mkdir(name)

%save good stuff
cmmd=sprintf('!cp  %s %s %s %s %s%sQA_maps', name_dt_Bdouble, name_md_nii, name_fa_nii, name_v1_nii,output_folder,filesep); 
eval(cmmd)

%save QA report
write_file = [output_folder filesep Fname '_' strrep(P,filesep,'-') '.pdf' ];
ps2pdf('psfile',testps,'pdffile',write_file);

cmmd=sprintf('!cp %s%spowerData.mat %s',tmp,filesep,trble);
eval(cmmd)

close all
% 
cmmd=sprintf('!rm -r %s',tmp);
eval(cmmd)
progress='Pipeline Completed'
  
       
   

