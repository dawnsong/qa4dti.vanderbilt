   %%%%%%%%%%PLOT

    a_slice = squeeze(cc(:,:,round(Nz/2),:));
    c_slice = squeeze(cc(round(Nx/2),:,:,:)); c_slice=flipdim(c_slice,2); c_slice=permute(c_slice,[2 1 3]);
    s_slice = squeeze(cc(:,round(Ny/2),:,:)); s_slice=flipdim(s_slice,1); s_slice=flipdim(s_slice,2); s_slice=permute(s_slice,[2 1 3]);
    
    % Plot Coronal
    axes(handles3.axes4(3))
     imagesc(c_slice)
     hold on
    set(gca,'XTick',[],'YTick',[])
    xlabel('INFERIOR','fontweight','demi','fontsize',12)
    ylabel('LEFT','fontweight','demi','fontsize',12)
    labelTitle = sprintf('Coronal Slice, #%02d',round(Nx/2));
    title(labelTitle,'fontweight','demi','fontsize',12)
    pos=get(gca,'Position');
    
    %Plot axial slice
    axes(handles3.axes4(2))
    imagesc(a_slice)
    hold on
    set(gca,'XTick',[],'YTick',[])
    xlabel('POSTERIOR','fontweight','demi','fontsize',12)
    %ylabel('LEFT','fontweight','demi','fontsize',12)
    labelTitle = sprintf('Axial Slice, #%02d',round(Nz/2));
    title(labelTitle,'fontweight','demi','fontsize',12)
   p2=get(gca,'Position'); p2(2:end)=pos(2:end); set(gca,'Position',p2);
   
    % Sagittal Slice
     axes(handles3.axes4(1))
    imagesc(s_slice)
    hold on
    
    xlabel('INFERIOR','fontweight','demi','fontsize',12)
    ylabel('ANTERIOR','fontweight','demi','fontsize',12);
    labelTitle = sprintf('Sagittal Slice, #%02d',round(Ny/2));
    title(labelTitle,'fontweight','demi','fontsize',12)
    p2=get(gca,'Position'); p2(2:end)=pos(2:end); set(gca,'Position',p2);
      set(gca,'XTick',[],'YTick',[],'Yaxislocation','right')