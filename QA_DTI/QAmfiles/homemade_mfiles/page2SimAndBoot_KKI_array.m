
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%BOXPLOT BOOT-SIGMA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FAbootM=zeros(size(NewMask,1),1); %size brain 
FAbootM(NewMask,:)=permute(ffstd,[2 1]);  %% FAboot in proper location in brain
doubleMask=maskROI+repmat(NewMask,1,25); doubleMask=(doubleMask==2);
%%%%%put data in ROIs
FAbootROI=[];  
grp=[];
Fboot=zeros(1,25);
for k=1:25
    FAbootROI=[FAbootROI FAbootM(doubleMask(:,k))' Sigma_ds_KKI(:,k)'];
    Fboot(1,k)=median(FAbootM(doubleMask(:,k)));
    grp=[grp ones(1,sum(doubleMask(:,k)))*(2*k-1) ones(1,1000)*2*k];
end
axes(handles2.axes43)
h=boxplot(FAbootROI,grp,'orientation','horizontal');
nm=sprintf('%s/BoxplotsFAsigma',trble);
save(nm,'FAbootROI','grp')
clear Sigma_ds_KKI FAbootROI%____________________________________________________clear line
for k=1:25
    set(h(:,2*k-1),'color','b','linewidth',2);
    set(h(:,2*k),'color','k','linewidth',1);
end
q=cell2mat(get(h(3,:),'XData')); xll=max(max(q));
set(h(7,:),'visible','off')
set(h(1,:),'linestyle','-'); set(h(2,:),'linestyle','-');
fax3=get(gca,'Position'); fax3(2)=model(2); fax3(3)=model(3)*2; fax3(4)=model(4);
set(gca,'Position',fax3,'LineWidth',2,'Color','none','xlim',[0 xll+.1*xll],'fontsize',10,'fontweight','demi','Xaxislocation','top','YDir','reverse')

aT=get(gca,'XTickLabel'); 
set(gca,'YTickLabel',{' '})
%%%xlabel('\sigma_{FA}') title set in page2FAMD since text location will be
%%%stable for FA
aT=get(gca,'XTickLabel'); set(gca,'XTickLabel',[]);
bT=get(gca,'XTick'); 
text(bT,repmat(-1,length(bT),1),aT,'rotation',-90,'fontweight','demi')
set(gca,'YTickLabel',{' '})
clear h %_________________________________________________________________________________-clear line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%BOXPLOT SIMEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Broi=zeros(size(maskROI,1),1); %size brain
Broi(NewMask,:)=Bias; %put bias in proper place in brain;
Biasroi=[];
Broi_b=zeros(1,25);
for k=1:25
    Biasroi=[Biasroi Broi(doubleMask(:,k))' Bias_ds_KKI(:,k)'];
    Broi_b(1,k)=median(Broi(doubleMask(:,k)));
end
axes(handles2.axes51)
h=boxplot(Biasroi,grp,'orientation','horizontal');
nm=sprintf('%s/BoxplotsBias',trble);
save(nm,'Biasroi','grp')
clear Bias_ds_KKI Biasroi %___________________________________________________________clear line
for k=1:25
    set(h(:,2*k-1),'color','b','linewidth',2);
    set(h(:,2*k),'color','k','linewidth',1);
end
q=cell2mat(get(h(3,:),'XData')); xml=max(max(q)); 
q=cell2mat(get(h(4,:),'XData')); xll=min(min(q));
set(h(7,:),'visible','off')
set(h(1,:),'linestyle','-'); set(h(2,:),'linestyle','-');
fax3=get(gca,'Position'); fax3(2)=model(2); fax3(3)=model(3)*2; fax3(4)=model(4);
set(gca,'Position',fax3,'YDir','reverse','LineWidth',2,'Color','none','xlim',[xll-.1*xll xml+.1*xml],'fontsize',10,'fontweight','demi','Xaxislocation','top')

aT=get(gca,'XTickLabel'); 
set(gca,'YTickLabel',{' '})
%%%xlabel('\sigma_{FA}') title set in page2FAMD since text location will be
%%%stable for FA
set(gca,'XTickLabel',[]);
bT=get(gca,'XTick'); 
text(bT,repmat(-1,length(bT),1),aT,'rotation',-90,'fontweight','demi')
set(gca,'YTickLabel',{' '})
clear h %_______________________________________________________________________________clear line