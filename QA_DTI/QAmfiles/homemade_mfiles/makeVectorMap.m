function makeVectorMap(slice,coron,dtmap, av, axial_handle, coronal_handle,mask,tmp,ks)
times=.75;
%%% AXIAL SLICE
%determine image resize for blow-up view
amask=mask(:,:,slice);
c=sum(amask,2)'; center_x=(sum((1:max(size(c))).*c))/sum(c);
top=max(find(c)); bottom=min(find(c));
c=sum(amask,1);center_y=(sum((1:max(size(c))).*c))/sum(c);
left=max(find(c)); right=min(find(c)); 
new_diamx=((center_x-bottom)+(top-center_x))*times;
new_diamy=((left-center_y)+(center_y-right))*times;
xsize=[center_x-new_diamx/2 center_x+new_diamx/2];
ysize=[center_y-new_diamy/2 center_y+new_diamy/2];
xsize=round(xsize); ysize=round(ysize);

dt=flipdim(dtmap,2); dt=permute(dt,[2 1 3 4]);
dt=flipdim(dt,1);
svname=sprintf('%s/faV1.tif',tmp);
imwrite(double(squeeze(dt(:,:,slice,:))),svname,'tif');
bg=imread(svname);
axes(axial_handle)
imagesc(bg);
hold on
set(gca,'YDir','normal')
set(gca,'XDir','normal')
ax=squeeze(av(:,:,slice,1:2));
for xi=1:size(bg,2) %%imshow flips y and j so x-is on 'y' in order to show correctly
    for yj=1:size(bg,1)
        m=ax(xi,yj,2)/ax(xi,yj,1);
        mid_y=yj; mid_x=xi;
        dim=sqrt(ax(xi,yj,2)^2+ax(xi,yj,1)^2); %percent in x,y plane vs z
        if abs(m)>1
            y1=yj-.5*dim; y2=yj+.5*dim;
            x1=(y1-mid_y+m*mid_x)/m;
            x2=(y2-mid_y+m*mid_x)/m;
        else
            x1=xi-.5*dim; x2=xi+.5*dim;
            y1=m*(x1-mid_x)+mid_y;
            y2=m*(x2-mid_x)+mid_y;
        end
        plot([x1 x2],[y1 y2],'w')
        hold on
    end
end
%set(gca,'xlim',[ceil((1+size(bg,2)*.2)) floor(size(bg,2)*.8)],'ylim',[ceil(1+size(bg,1)*.2) floor(size(bg,2)*.8)])
if ks==1
set(gca,'xlim',xsize,'ylim',ysize)
end
if ks==2
    ysize=[round((ysize(2)-ysize(1))/2)+ysize(1) ysize(2)];
    xsize=[xsize(1)+(xsize(2)-xsize(1))*.2 xsize(2)-(xsize(2)-xsize(1))*.2];
    set(gca,'xlim',xsize,'ylim',ysize)
end
if ks==3
    ysize=[ysize(1) round((ysize(2)-ysize(1))/2)+ysize(1)];
    xsize=[xsize(1)+(xsize(2)-xsize(1))*.2 xsize(2)-(xsize(2)-xsize(1))*.2];
    set(gca,'xlim',xsize,'ylim',ysize)
end
%%CORONAL SLICE

%determine image resize for blow-up view
cmask=squeeze(mask(:,coron,:));
c=sum(cmask,2)'; center_x=(sum((1:max(size(c))).*c))/sum(c);
top=max(find(c)); bottom=min(find(c));
c=sum(cmask,1);center_y=(sum((1:max(size(c))).*c))/sum(c);
left=max(find(c)); right=min(find(c)); 
new_diamx=((center_x-bottom)+(top-center_x))*times;
new_diamy=((left-center_y)+(center_y-right))*times;
xsize=[center_x-new_diamx/2 center_x+new_diamx/2];
ysize=[center_y-new_diamy/2 center_y+new_diamy/2];
xsize=round(xsize); ysize=round(ysize);


dt=squeeze(dtmap(:,coron,:,:));
%dt=flipdim(dt,2);
dt=permute(dt,[2 1 3 4]);
svname=sprintf('%s/faV1_coron.tif',tmp);
imwrite(double(dt),svname,'tif');
bg=imread(svname);
axes(coronal_handle)
imagesc(bg);
hold on
set(gca,'YDir','normal')
set(gca,'XDir','normal')
co = squeeze(av(:,coron,:,1:2:3));
for xi=1:size(bg,2)
    for zk=1:size(bg,1)
        m=co(xi,zk,2)/co(xi,zk,1);
        mid_z=zk; mid_x=xi;
        dim=sqrt(co(xi,zk,2)^2+co(xi,zk,1)^2); %percent in x,z plane vs y
        if abs(m)>1
            z1=zk-.5*dim; z2=zk+.5*dim;
            x1=(z1-mid_z+m*mid_x)/m;
            x2=(z2-mid_z+m*mid_x)/m;
        else
            x1=xi-.5*dim; x2=xi+.5*dim;
            z1=m*(x1-mid_x)+mid_z;
            z2=m*(x2-mid_x)+mid_z;
        end
        plot([x1 x2],[z1 z2],'w')
        hold on
    end
end

% set(gca,'xlim',[ceil((1+size(bg,2)*.2)) floor(size(bg,2)*.8)],'ylim',[ceil(1+size(bg,1)*.2) floor(size(bg,2)*.8)])

if ks==1
set(gca,'xlim',xsize,'ylim',ysize)
end
if ks==2
    ysize=[round((ysize(2)-ysize(1))/2)+ysize(1) ysize(2)];
    xsize=[xsize(1)+(xsize(2)-xsize(1))*.2 xsize(2)-(xsize(2)-xsize(1))*.2];
    set(gca,'xlim',xsize,'ylim',ysize)
end
if ks==3
     ysize=[ysize(1) round((ysize(2)-ysize(1))/2)+ysize(1)];
    xsize=[xsize(1)+(xsize(2)-xsize(1))*.2 xsize(2)-(xsize(2)-xsize(1))*.2];
    set(gca,'xlim',xsize,'ylim',ysize)
end
