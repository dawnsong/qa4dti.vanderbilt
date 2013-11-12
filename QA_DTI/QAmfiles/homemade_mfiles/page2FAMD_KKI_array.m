load BoxPlotKKI %% later can be replaced by intra-study subjects
FA=FA(:); MD=MD(:)*10^3;
FAroi=[]; MDroi=[]; grp=[];

for k=1:25
    FAroi=[FAroi FA(maskROI(:,k))' FA_ds_KKI(:,k)'];%%% values from this data
    MDroi=[MDroi MD(maskROI(:,k))' MD_ds_KKI(:,k)'*1000];
    grp=[grp ones(1,sum(maskROI(:,k)))*(2*k-1) ones(1,1000)*2*k];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PLOT JPEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tst=8
Y=imread('brain-labels_only.tif');
tst=9
axes(handles2.axes40)
    Y2=permute(Y,[2 1 3]);
    Y2=flipdim(Y2,2);
imagesc(Y2)
tst=10
set(gca,'XTick',[],'YTick',[])
model=get(gca,'Position');
clear Y Y2 %____________________________________________________________________________clear line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%BOXPLOT MD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles2.axes41)
set(gca,'YDir','reverse')
tst=11

h=boxplot(MDroi,grp,'orientation','horizontal');
tst=12
for k=1:25
    set(h(:,2*k-1),'color','b','linewidth',2);
    set(h(:,2*k),'color','k','linewidth',1);
end
set(h(7,:),'visible','off')
set(h(1,:),'linestyle','-'); set(h(2,:),'linestyle','-');
q=cell2mat(get(h(3,:),'XData')); xll=max(max(q));
fax2=get(gca,'Position'); fax2(2)=model(2); fax2(3)=model(3)*2; fax2(4)=model(4);
set(gca,'Position',fax2,'YDir','reverse','LineWidth',2,'Color','none','xlim',[0 xll+.1*xll],'fontsize',10,'fontweight','demi','Xaxislocation','top')
aT=get(gca,'XTickLabel'); set(gca,'XTickLabel',[]);
aTn=aT(1:end-1,:);
bT=get(gca,'XTick'); bTn=bT(1:end-1);
tst=13
text(bTn,repmat(-.2,length(bTn),1),aTn,'rotation',-90,'fontweight','demi')
set(gca,'YTickLabel',{' '}) %,'Ticklength',[0 0])
h=text(1.8,-2,'MD');
set(h,'rotation',-90,'fontsize',12,'fontweight','demi');
h=text(1.4,-2,'x10^-^3');
set(h,'rotation',-90,'fontsize',10,'fontweight','demi');
h=text(.7,-2.5,'mm^2/s');
set(h,'rotation',-90,'fontsize',10,'fontweight','demi');
tst=14
nm=sprintf('%s/BoxplotsMD',trble);
save(nm,'MDroi','grp')
clear MDroi h aT aTn bT MD_ds_KKI %_____________________________________________________________clear line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%BOXPLOT FA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles2.axes44)
tst=15
clear h
tst=16
h=boxplot(FAroi,grp,'orientation','horizontal');
tst=17
for k=1:25
    set(h(:,2*k-1),'color','b','linewidth',2);
    set(h(:,2*k),'color','k','linewidth',1);
end
set(h(7,:),'visible','off')
set(h(1,:),'linestyle','-'); set(h(2,:),'linestyle','-');
q=cell2mat(get(h(3,:),'XData')); xll=max(max(q));
fax=get(gca,'Position'); fax(2)=model(2); fax(3)=model(3)*2; fax(4)=model(4);
set(gca,'Position',fax,'LineWidth',2,'fontsize',10,'fontweight','demi','Xaxislocation','top','YDir','reverse')
h=text(.5,-2,'FA');
h2=text(2.0,-2,'\sigma_{FA}'); 
h3=text(3.2,-2,'bias');
set(h,'fontweight','demi');set(h2,'fontweight','demi');set(h3,'fontweight','demi');
aT=get(gca,'XTickLabel'); set(gca,'XTickLabel',[]);
bT=get(gca,'XTick'); 
text(bT,repmat(-.6,length(bT),1),aT,'rotation',-90,'fontweight','demi')
set(gca,'YTickLabel',{' '})
tst=18
nm=sprintf('%s/BoxplotsFA',trble);
save(nm,'FAroi','grp')
clear FAroi FA_ds_KKI h%_________________________________________________________clear line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PLOT DASHED LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles2.axes42)
tst=19
for k=1:24
    plot(0:10,ones(11,1)*k,'--k')
    hold on
    if mod(k,2)==0
    plot(0:10,ones(11,1)*k,'-k')
    end
end
fax2(3)=model(3)*8;
set(gca,'Color','none','Position',fax2,'YTick',[],'XTick',[],'YTicklabel','')
tst=20
