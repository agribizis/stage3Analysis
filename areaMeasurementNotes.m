load('meanPSF_b2_all_rOnly_250k.mat')
fS=6;
G1=GcampZs(1:250000,:);
G2=GcampZs(250001:500000,:);
G3=GcampZs(500001:750000,:);
G4=GcampZs(750001:1000000,:);
G5=GcampZs(1000001:1250000,:);
G6=GcampZs(1250001:1500000,:);
test=ones(length(meanPSF),1);
meanPSF=horzcat(meanPSF,test);
m1=meanPSF(1:250000,:);
m2=meanPSF(250001:500000,:);
m3=meanPSF(500001:750000,:);
m4=meanPSF(750001:1000000,:);
m5=meanPSF(1000001:1250000,:);
m6=meanPSF(1250001:1500000,:);
[meanPSF_G,bint,r,rint,stats]=regress(G1,m1);
sqSz=(fS*2+1)*(fS*2+1);
psf1=meanPSF_G(1:sqSz);
%psf2=vertcat(meanPSF_G(290:329),0,meanPSF_G(330:369));
psf2=vertcat(meanPSF_G((sqSz+1):(1.5*(sqSz)-.5)),0,meanPSF_G((1.5*(sqSz)-.5):length(meanPSF_G)-2));
psf1=reshape(psf1, fS*2+1,fS*2+1);
psf2=reshape(psf2, fS*2+1,fS*2+1);
figure; imagesc(psf1)
figure; imagesc(psf2)
psf1_1=psf1; psf2_1=psf2;
[meanPSF_G,bint,r,rint,stats]=regress(G2,m2);
sqSz=(fS*2+1)*(fS*2+1);
psf1=meanPSF_G(1:sqSz);
%psf2=vertcat(meanPSF_G(290:329),0,meanPSF_G(330:369));
psf2=vertcat(meanPSF_G((sqSz+1):(1.5*(sqSz)-.5)),0,meanPSF_G((1.5*(sqSz)-.5):length(meanPSF_G)-2));
psf1=reshape(psf1, fS*2+1,fS*2+1);
psf2=reshape(psf2, fS*2+1,fS*2+1);
figure; imagesc(psf1)
figure; imagesc(psf2)
psf1_2=psf1; psf2_2=psf2;
[meanPSF_G,bint,r,rint,stats]=regress(G3,m3);
sqSz=(fS*2+1)*(fS*2+1);
psf1=meanPSF_G(1:sqSz);
%psf2=vertcat(meanPSF_G(290:329),0,meanPSF_G(330:369));
psf2=vertcat(meanPSF_G((sqSz+1):(1.5*(sqSz)-.5)),0,meanPSF_G((1.5*(sqSz)-.5):length(meanPSF_G)-2));
psf1=reshape(psf1, fS*2+1,fS*2+1);
psf2=reshape(psf2, fS*2+1,fS*2+1);
figure; imagesc(psf1)
figure; imagesc(psf2)
psf1_3=psf1; psf2_3=psf2;
[meanPSF_G,bint,r,rint,stats]=regress(G4,m4);
sqSz=(fS*2+1)*(fS*2+1);
psf1=meanPSF_G(1:sqSz);
%psf2=vertcat(meanPSF_G(290:329),0,meanPSF_G(330:369));
psf2=vertcat(meanPSF_G((sqSz+1):(1.5*(sqSz)-.5)),0,meanPSF_G((1.5*(sqSz)-.5):length(meanPSF_G)-2));
psf1=reshape(psf1, fS*2+1,fS*2+1);
psf2=reshape(psf2, fS*2+1,fS*2+1);
figure; imagesc(psf1)
figure; imagesc(psf2)
psf1_4=psf1; psf2_4=psf2;
[meanPSF_G,bint,r,rint,stats]=regress(G5,m5);
sqSz=(fS*2+1)*(fS*2+1);
psf1=meanPSF_G(1:sqSz);
%psf2=vertcat(meanPSF_G(290:329),0,meanPSF_G(330:369));
psf2=vertcat(meanPSF_G((sqSz+1):(1.5*(sqSz)-.5)),0,meanPSF_G((1.5*(sqSz)-.5):length(meanPSF_G)-2));
psf1=reshape(psf1, fS*2+1,fS*2+1);
psf2=reshape(psf2, fS*2+1,fS*2+1);
figure; imagesc(psf1)
figure; imagesc(psf2)
psf1_5=psf1; psf2_5=psf2;
[meanPSF_G,bint,r,rint,stats]=regress(G6,m6);
sqSz=(fS*2+1)*(fS*2+1);
psf1=meanPSF_G(1:sqSz);
%psf2=vertcat(meanPSF_G(290:329),0,meanPSF_G(330:369));
psf2=vertcat(meanPSF_G((sqSz+1):(1.5*(sqSz)-.5)),0,meanPSF_G((1.5*(sqSz)-.5):length(meanPSF_G)-2));
psf1=reshape(psf1, fS*2+1,fS*2+1);
psf2=reshape(psf2, fS*2+1,fS*2+1);
figure; imagesc(psf1)
figure; imagesc(psf2)
psf1_6=psf1; psf2_6=psf2;
clear
load('b2bothFiltersRegressAll.mat')
psf1_1(psf1_1<0)=0;
psf1_2(psf1_2<0)=0;
psf1_3(psf1_3<0)=0;
psf1_4(psf1_4<0)=0;
psf1_5(psf1_5<0)=0;
psf1_6(psf1_6<0)=0;
sum1=sum(psf1_1(:));
sum2=sum(psf1_2(:));
sum3=sum(psf1_3(:));
sum4=sum(psf1_4(:));
sum5=sum(psf1_5(:));
sum6=sum(psf1_6(:));
psf2_1(psf2_1<0)=0;
psf2_2(psf2_2<0)=0;
psf2_3(psf2_3<0)=0;
psf2_4(psf2_4<0)=0;
psf2_5(psf2_5<0)=0;
psf2_6(psf2_6<0)=0;
sum21=sum(psf2_1(:));
sum22=sum(psf2_2(:));
sum23=sum(psf2_3(:));
sum24=sum(psf2_4(:));
sum25=sum(psf2_5(:));
sum26=sum(psf2_6(:));
A1=find(psf1_1>0);
A2=find(psf1_2>0);
A3=find(psf1_3>0);
A4=find(psf1_4>0);
A5=find(psf1_5>0);
A6=find(psf1_6>0);
preSums=[sum1 sum2 sum3 sum4 sum5 sum6 sum7 sum8];
postSums=[sum21 sum22 sum23 sum24 sum25 sum26 sum27 sum28];
