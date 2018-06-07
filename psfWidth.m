

psf1_1(psf1_1<0)=0;
psf1_2(psf1_2<0)=0;
psf1_3(psf1_3<0)=0;
psf1_4(psf1_4<0)=0;
psf1_5(psf1_5<0)=0;
psf1_6(psf1_6<0)=0;
psf1_7(psf1_7<0)=0;
psf1_8(psf1_8<0)=0;
psf2_1(psf2_1<0)=0;
psf2_2(psf2_2<0)=0;
psf2_3(psf2_3<0)=0;
psf2_4(psf2_4<0)=0;
psf2_5(psf2_5<0)=0;
psf2_6(psf2_6<0)=0;
psf2_7(psf2_7<0)=0;
psf2_8(psf2_8<0)=0;


% filt=psf2_8;
% A=sum(sum(filt));
% mult=1:(fS*2+1);
% x=sum(sum(filt).*mult)./A;
% y=sum(sum(filt').*mult)./A;
% x=round(x); y=round(y);
% 
% filt(filt>0)=1;
% sx28=sum(filt(x,:));
% sy28=sum(filt(:,y));


% filt=psf1_1;
% 
% maxF=max(filt(:));
% filt=filt./maxF;
% 
% A=sum(sum(filt));
% mult=1:(fS*2+1);
% x=sum(sum(filt).*mult)./A;
% y=sum(sum(filt').*mult)./A;
% x=round(x); y=round(y);
% 
% sx26=std(filt(x,:));
% sy26=std(filt(:,y));



filt=psf2_8;
A=sum(sum(filt));
mult=1:(fS*2+1);
x=sum(sum(filt).*mult)./A;
y=sum(sum(filt').*mult)./A;

wSum=0;
for i=1:(fS*2+1)
    for j=1:(fS*2+1)
        wSum=wSum+(i-x)^2*(filt(i,j)/A)+(j-y)^2*(filt(i,j)/A);
    end
end

wSum28=wSum;





group = [repmat({'S2'}, 8, 1); repmat({'S3'}, 8, 1); repmat({'beta2'}, 6, 1)];
x=meas;
gplotmatrix(x,[],group,[],'+xo')
[d, p, stats] = manova1(x,group);
c1 = stats.canon(:,1);
c2 = stats.canon(:,2);
figure()
gscatter(c2,c1,group,[],'oxs')