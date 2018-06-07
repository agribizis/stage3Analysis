function error = modelTestFitSimpBi(vars, fS, GcampZs, Rcamps)

%vars=[5 8 5 -4 6 1.5];

mu = [vars(1) vars(2)];
x = 1:(fS*2)+1;
[X1,X2] = meshgrid(x,x);
%Sigma = [4 -2; -2 8];
Sigma=[vars(3) vars(4); vars(4) vars(5)];
if all(eig(Sigma) > eps)
    F = mvnpdf([X1(:) X2(:)],mu,Sigma);
    F = reshape(F,length(x),length(x));
    F = F*vars(6);
    %figure; imagesc(F)
    h=F;
    h=reshape(h,1,169);

    G_pred=Rcamps*h';

    err=GcampZs-G_pred;
    err=err.^2;
    error=sum(err);


    R = min(min(corrcoef(GcampZs,G_pred)));
    R1=1-R;
    
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