function error = modelTestFitSimpBi_both(vars, fS, GcampZs, Gcamps, Rcamps)

%vars=[5 8 5 -4 6 1.5 5 8 5 -4 6 1.5];

mu1 = [vars(1) vars(2)];
x = 1:(fS*2)+1;
[X1,X2] = meshgrid(x,x);
Sigma1=[vars(3) vars(4); vars(4) vars(5)];
if all(eig(Sigma1) > eps)
    F = mvnpdf([X1(:) X2(:)],mu1,Sigma1);
    F = reshape(F,length(x),length(x));
    F = F*vars(6);
    
    h1=F;
    h1=reshape(h1,1,169);

    mu2 = [vars(7) vars(8)];
    Sigma2=[vars(9) vars(10); vars(10) vars(11)];
    if all(eig(Sigma2) > eps)
        F = mvnpdf([X1(:) X2(:)],mu2,Sigma2);
        F = reshape(F,length(x),length(x));
        F = F*vars(12);
        h2=F;
        h2(fS+1, fS+1)=0;
        h2=reshape(h2,1,169);
        
        
        G_pred=Rcamps*h1'+Gcamps*h2';

        err=GcampZs-G_pred;
        err=err.^2;
        error=sum(err);

        R = min(min(corrcoef(GcampZs,G_pred)));
        R1=1-R;
    else
        error=inf; R1=inf;
    end
    
else
    error=inf; R1=inf;
end


% if vars(2)<0
%     R1=100000;
% end
vars
R1
1-R1

end