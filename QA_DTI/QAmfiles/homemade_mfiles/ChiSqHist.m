chi_sq_p=sum(sd,4); chi_sq_p=chi_sq_p(:);csqp=chi_sq_p(brain);
bb=(csqp>.3); csqp(bb)=.3;
 xx=0:.005:10;  %with wide x, histnorm stays stable 
h=histnorm(csqp,xx);
 axes(handles.axes39)
 plot(xx,h,'linewidth',2);
 ym=get(gca,'ylim'); ym=ym(2);
 h=h./ym;
 plot(xx,h,'linewidth',2);
 set(gca,'xlim',[0 .3])
 set(gca,'xlim',[0 .3])
 set(gca,'fontweight','bold')
 mChi=median(csqp);
 stChi=std(csqp);
 xlabel('Chi Squared [1]')
 
 %make insert
 yh=get(gca,'ylim'); ys=(yh(2)*.15);
 ho=figure(200);
 plot(xx(25:61),h(25:61));
 limmy=get(gca,'ylim');
 close(ho)
 axes(handles.axes39)
 hold on
 plot(xx(21:61),(h(21:61)*yh(2)/limmy(2))+ys,'linewidth',2)
 set(gca,'ylim',yh)
  hhh=h(21:61); xxx=xx(21:61);h2=(hhh>limmy(2)*.55);
 plot(xxx(h2),(hhh(h2)*yh(2)/limmy(2))+ys,'w','linewidth',2)
 plot([xx(21),xx(61)],[ys,ys],'k','linewidth',2)
 plot([xx(21),xx(21)],[ys,(yh(2)*.66)+ys],'k','linewidth',2)
 set(gca,'box','off')
 h=text(.085,ys,'0');
 set(h,'fontweight','demi')
 hh=sprintf('%g',((yh(2)*.66-ys)*limmy(2))/yh(2));
 loc= yh(2)*.7+ys;
 h=text(.07,loc,hh);
 set(h,'fontweight','demi')
%  h=text(.1,ys-.033*yh(2),'0.1');
%  set(h,'fontweight','demi')
%   h=text(.29,ys-.033*yh(2),'0.3');
%  set(h,'fontweight','demi')
%  plot([xx(41) xx(41)],[round(ys*.8) round(ys*.3)],'linewidth',2)
%  plot([xx(41) xx(40)],[round(ys*.8) round(ys*.65)],'linewidth',2)
%  plot([xx(41) xx(42)],[round(ys*.8) round(ys*.65)],'linewidth',2)
 hh=sprintf('x %d',round(yh(2)/limmy(2)));
 h=text(.25,yh(2)*.5+ys,hh);
 set(h,'fontweight','demi')