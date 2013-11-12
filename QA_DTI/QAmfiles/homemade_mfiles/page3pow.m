%%% PLot mid slice of FA, MD, sigma_FA and bias FA
aa=.015; ay=.015; sy=.02;
 mi=round(Nz/2);
Fm=reshape(FA,Nx,Ny,Nz); Mm=reshape(MD,Nx,Ny,Nz);
axes(handles3.axes40)
imagesc(flipdim(Mm(:,:,mi),2)')
colorbar('fontweight','demi')
 set(gca,'XTick',[],'YTick',[],'clim',[0 3])
title('MD 10^{-3} mm^2/s','fontsize',12,'fontweight','demi')
set(gca,'Position',[ 0.0172-aa 0.8155-ay 0.2000 0.1423]);
cbar=findobj('tag','Colorbar'); cbar=cbar(1);
set(cbar(1,:),'Position',[.22-aa .8184-ay .0150 .1433-sy]);

clear Mm MD %_____________________________________________________clear line

 axes(handles3.axes41)
imagesc(flipdim(Fm(:,:,mi),2)')
colorbar('fontweight','demi')
 set(gca,'XTick',[],'YTick',[])
title('FA','fontsize',12,'fontweight','demi')
p=[.255-aa .8155-ay .1980 .1423];set(gca,'Position',p);
cbar=findobj('tag','Colorbar'); cbar=cbar(1);
set(cbar,'Position',[.455-aa .8184-ay .0150 .1453-sy]);
set(gca,'clim',[0 1]);
clear Fm FA %________________________________________________________________clear line

medianBootSim % get median values to represent entire region
axes(handles3.axes42)
imagesc(flipdim(Sm(:,:,mi),2)')
colorbar('fontweight','demi')
set(gca,'XTick',[],'YTick',[])
title('\sigma_{FA} from Bootstrap [3]','fontsize',12,'fontweight','demi')
p=[.505-aa .8155-ay .1980 .1423];set(gca,'Position',p);
cbar=findobj('tag','Colorbar'); cbar=cbar(1);
set(cbar,'Position',[.705-aa .8184-ay .0150 .1453-sy]);
lable=sprintf('median from %i%s sampling',round(percentage),'%');
xlabel(lable)
set(gca,'clim',[0 .05]);
clear Sm %__________________________________________________________________clear line

axes(handles3.axes43)
imagesc(flipdim(Bm(:,:,mi),2)')
colorbar('fontweight','demi')
 set(gca,'XTick',[],'YTick',[])
title('Bias of FA from SIMEX [4]','fontsize',12,'fontweight','demi')
p=[.76-aa .8155-ay .1980 .1423]; set(gca,'Position',p);
cbar=findobj('tag','Colorbar'); cbar=cbar(1);
set(cbar,'Position',[.96-aa .8184-ay .0150 .1453-sy]);
colormap('jet')
lable=sprintf('median from %i%s sampling',round(percentage),'%');
xlabel(lable)
set(gca,'clim',[0 .08]);
clear Bm %_____________________________________________________________________clear line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Calculate Power and Plot Power GM and WM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%NO bias, RHS and LHS GM
ES=-.1:.001:.1;
us=NewMask+maskROI(:,2)+maskROI(:,3); us=(us==2); %%us is voxels used for BOOt that are right hand say GM and LHS GM
sig=FAbootM(us);

bias=zeros(sum(us),1); % no bias here
pow=zeros(sum(us),length(ES));

n=5;
for esize=1:length(ES)
    pow(:,esize) = two_sided_t_pow(ES(esize),sig,bias,n,.05);
end
p=median(pow);
axes(handles3.axes44)
plot(ES,p,'k','linewidth',2)
hold on

n=15;
for esize=1:length(ES)
    pow(:,esize) = two_sided_t_pow(ES(esize),sig,bias,n,.05);
end
p=median(pow);
axes(handles3.axes44)
plot(ES,p,'r','linewidth',2)
hold on

n=30;
for esize=1:length(ES)
    pow(:,esize) = two_sided_t_pow(ES(esize),sig,bias,n,.05);
end
p=median(pow);
axes(handles3.axes44)
plot(ES,p,'b','linewidth',2)
hold on
ylabel('Power','fontsize',10,'fontweight','demi')
xlabel('effect size (FA)','fontsize',10,'fontweight','demi')
set(gca,'XTick',[-.1 -.05 0 .05],'XTickLabel',['-0.1|-.05|0|.05|'],'fontweight','demi');
set(gca,'linewidth',2,'box','off')
set(gca,'TickLength',[.02 .1])

%%%NO BIAS RHS and LHS White Matter 
ES=-.1:.001:.1;
us=NewMask+maskROI(:,4)+maskROI(:,5); us=(us==2);
sig=FAbootM(us);
bias=zeros(sum(us),1);
pow=zeros(sum(us),length(ES));

n=5;
for esize=1:length(ES)
    pow(:,esize) = two_sided_t_pow(ES(esize),sig,bias,n,.05);
end
p=median(pow);
axes(handles3.axes45)
plot(ES,p,'k','linewidth',2)
hold on

n=15;
for esize=1:length(ES)
    pow(:,esize) = two_sided_t_pow(ES(esize),sig,bias,n,.05);
end
p=median(pow);
axes(handles3.axes45)
plot(ES,p,'r','linewidth',2)
hold on

n=30;
for esize=1:length(ES)
    pow(:,esize) = two_sided_t_pow(ES(esize),sig,bias,n,.05);
end
p=median(pow);
axes(handles3.axes45)
plot(ES,p,'b','linewidth',2)
hold on
set(gca,'YTickLabel','');
set(gca,'XTick',[-.1 -.05 0 .05],'XTickLabel',['+/-0.1|-.05|0|.05|'],'fontweight','demi');
set(gca,'TickLength',[.02 .1],'Yaxislocation','right','linewidth',2,'box','off')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% WITH BIAS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%BIAS white matter
ES=-.1:.001:.1;
us=NewMask+maskROI(:,2)+maskROI(:,3); us=(us==2);
sig=FAbootM(us);
bias =Broi(us); 
pow=zeros(sum(us),length(ES));
n=5;
for esize=1:length(ES)
    pow(:,esize) = two_sided_t_pow(ES(esize),sig,bias,n,.05);
end
p=median(pow);
axes(handles3.axes46)
plot(ES,p,'k','linewidth',2)
hold on
set(gca,'XTick',[-.1 -.05 0 .05],'XTickLabel',['+/-0.1|-.05|0|.05|'],'YTickLabel','','fontweight','demi');

n=15;
for esize=1:length(ES)
    pow(:,esize) = two_sided_t_pow(ES(esize),sig,bias,n,.05);
end
p=median(pow);
axes(handles3.axes46)
plot(ES,p,'r','linewidth',2)
hold on

n=30;
for esize=1:length(ES)
    pow(:,esize) = two_sided_t_pow(ES(esize),sig,bias,n,.05);
end
p=median(pow);
axes(handles3.axes46)
plot(ES,p,'b','linewidth',2)
hold on
set(gca,'linewidth',2,'Box','off')
set(gca,'TickLength',[.02 .1])

%%%Bias WM
ES=-.1:.001:.1;
us=NewMask+maskROI(:,4)+maskROI(:,5); us=(us==2);
sig=FAbootM(us);
bias =Broi(us); 
pow=zeros(sum(us),length(ES));

n=5;
for esize=1:length(ES)
    pow(:,esize) = two_sided_t_pow(ES(esize),sig,bias,n,.05);
end
p=median(pow);
axes(handles3.axes47)
plot(ES,p,'k','linewidth',2)
hold on
set(gca,'YTickLabel','','linewidth',2)
set(gca,'XTick',[-.1 -.05 0 .05 .1],'XTickLabel',['+/-0.1|-.05|0|.05|0.1'],'fontweight','demi');

n=15;
for esize=1:length(ES)
    pow(:,esize) = two_sided_t_pow(ES(esize),sig,bias,n,.05);
end
p=median(pow);
axes(handles3.axes47)
plot(ES,p,'r','linewidth',2)
hold on

n=30;
for esize=1:length(ES)
    pow(:,esize) = two_sided_t_pow(ES(esize),sig,bias,n,.05);
end
p=median(pow);
axes(handles3.axes47)
plot(ES,p,'b','linewidth',2)
hold on
set(gca,'TickLength',[.02 .1],'Yaxislocation','right','linewidth',2,'box','off')
savename=sprintf('%s/powerData.mat',tmp);
save(savename,'NewMask','maskROI','FAbootM','Broi')
clear FAbootM esize pow ES p %_____________________________________________________clear line