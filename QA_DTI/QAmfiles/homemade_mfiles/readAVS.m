%%read avscale output textfile
function [translation rotation]=readAVS(avs)
b=(avs=='='); tt=(avs=='T'); ss=(avs=='S');
loc=find(b==1); loc2=find(tt==1); loc3=find(ss==1);
%rotation
rot=str2num(avs(loc(1)+1:loc2(2)-1));
rotation=rot*360/(2*pi);
%translation
translation=str2num(avs(loc(2)+1:loc3(1)-1));