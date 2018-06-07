function error = modelTestFitSimp_both(vars, fS, GcampZs, Gcamps, Rcamps)

% Gcamps=meanPSF(:,170:end);
% test=zeros(length(meanPSF),1);
% Gcamps=[Gcamps(:,1:84) test Gcamps(:,85:end)];

h1=[];
h2=[];

amp=vars(1);
sigma=vars(2);
for k1=1:fS*2+1
    for k2=1:fS*2+1
        n1=k1-(fS+1);
        n2=k2-(fS+1);
        h1(k1,k2)=exp(-(n1^2+n2^2)^2/(2*sigma^2));
    end
end

h1 = h1/sum(sum(h1));
h1=amp*h1;

amp2=vars(3);
sigma2=vars(4);
for k1=1:fS*2+1
    for k2=1:fS*2+1
        n1=k1-(fS+1);
        n2=k2-(fS+1);
        h2(k1,k2)=exp(-(n1^2+n2^2)^2/(2*sigma2^2));
    end
end

h2 = h2/sum(sum(h2));
h2=amp2*h2;
h2(fS+1, fS+1)=0;

h1=reshape(h1,1,169);
h2=reshape(h2,1,169);

% DoG=(h1-h2);
% h=DoG;
% h=reshape(h,1,169);



G_pred=Rcamps*h1'+Gcamps*h2';

err=GcampZs-G_pred;
err2=err.^2;
error=sum(err2);


R = min(min(corrcoef(GcampZs,G_pred)));
R1=1-R;
if vars(1)<0 || vars(3)<0 || vars(2)<0 || vars(4)<0
    R1=100000; error=inf;
end

vars
R1
1-R1
end