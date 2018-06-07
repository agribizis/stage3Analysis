function error = modelTestFitSimp_Ronly(vars, fS, GcampZs, Rcamps)

h=[];

amp=vars(1);
sigma=vars(2);
for k1=1:fS*2+1
    for k2=1:fS*2+1
        n1=k1-(fS+1);
        n2=k2-(fS+1);
        h(k1,k2)=exp(-(n1^2+n2^2)^2/(2*sigma^2));
    end
end

h = h/sum(sum(h));
h=amp*h;
h=reshape(h,1,169);



G_pred=Rcamps*h';

err=GcampZs-G_pred;
err=err.^2;
error=sum(err);


R = min(min(corrcoef(GcampZs,G_pred)));
R1=1-R;
if vars(2)<0
    R1=100000;
end

vars
R1
1-R1
end