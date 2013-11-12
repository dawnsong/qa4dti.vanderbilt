function [FAsmx Bias FAc]=simex(R,bval_vec,grad_file,sigmaEst,n_bo,numsim)

vox =size(R,2);
nosim = numsim;  % accuracy decreases with variance and increases with N  err(mean(X))~sigma/sqrt(n)
b = bval_vec(1:end-1);
lambda = [0 2 4 6 8];
lambda_sqrt = sqrt(lambda);
nana=length(lambda_sqrt);
dwiObs=R(1:end-1,:); boObs=R(end,:);
grads=grad_file';
gtable=[grads.^2  2*grads(:,1).*grads(:,2) 2*grads(:,1).*grads(:,3) 2*grads(:,3).*grads(:,2)];

matlabpool close
try
    matlabpool open
    parfor jj=1:nana   
        sig = lambda_sqrt(jj)*sigmaEst;
        sim=nosim(jj);
        [FAsnr{jj} ]=dawn_DTIfit(vox,sim,dwiObs,boObs,sig,b,gtable, n_bo);
    end
catch err,
end
matlabpool close


 FA=zeros(vox,nana);

for j=1:nana
        FA(1:vox,j) = FAsnr{j};
end
 FAobs=FA(:,1);
 FAc=FAobs;
% % % % % % % % % % % % % % % %     
    %Calculate FAsimex
% % % % % % % % % % % % 
FAsmx=zeros(vox,1);

    for v=1:vox
        metric=FA(v,:);
        datax=double(lambda); 
        datay=double(metric);
        param= polyfit(datax,datay,2);
        FAsmx(v)=polyval(param, -1);      
    end
    Bias = FAobs-FAsmx;
    

