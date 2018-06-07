load '02_tmpMatsS3.mat'

Gcamp=zeros(sz(1)*sz(2),sz(3));
Gcamp(index,:)=Gcamp_tmp';
Rcamp=zeros(sz(1)*sz(2),sz(3));
Rcamp(index,:)=Rcamp_tmp';
Rcamp=reshape(Rcamp,sz);
Gcamp=reshape(Gcamp,sz);
%Gcamp=fliplr(Gcamp);
%Rcamp=fliplr(Rcamp);

Cors=zeros(sz(1)*2-1,sz(2)*2-1,sz(3));
parfor i=1:sz(3)
Cors(:,:,i)=xcorr2(Gcamp(:,:,i),Rcamp(:,:,i));
end
%figure; imagesc(mean(Cors,3))

[r,c]=find(mean(Cors,3)==max(max(mean(Cors,3))))
xshift=r-sz(1)
yshift=c-sz(2)

newRcamp=zeros(sz);

%newRcamp(xshift:end, yshift:end, :)=Rcamp(1:end-xshift+1, 1:end-yshift+1,:);
newRcamp(xshift:end, 1:end+yshift, :)=Rcamp(1:end-xshift+1, -yshift+1:end,:);
Rcamp=newRcamp;

Gcamp=imresize(Gcamp,0.625,'nearest');
Rcamp=imresize(Rcamp,0.625,'nearest');
Gcamp=imresize(Gcamp,0.5,'nearest');
Rcamp=imresize(Rcamp,0.5,'nearest');


meanG=mean(Gcamp,3);

index=find(meanG);

%constrain to pixels that are not near the edge of the image
index2=[];

[sROI] = ReadImageJROI('RoiSet4.zip');
mask = poly2mask(sROI{1,1}.mnCoordinates(:,1), sROI{1,1}.mnCoordinates(:,2),sz(1),640)+poly2mask(sROI{1,2}.mnCoordinates(:,1), sROI{1,2}.mnCoordinates(:,2),sz(1),640);
%mask = poly2mask(sROI{1,1}.mnCoordinates(:,1), sROI{1,1}.mnCoordinates(:,2)-270,sz(1),640)+poly2mask(sROI{1,2}.mnCoordinates(:,1), sROI{1,2}.mnCoordinates(:,2)-270,sz(1),640);
%changed sz(2) to 640, subtracted 270 if ROI is on the bottom...check this
%on next round
%mask = fliplr(mask); 
mask=imresize(mask,0.625,'nearest');
mask=imresize(mask,0.5,'nearest');
sz=size(Gcamp);
mask=mask(1:sz(1), 1:sz(2));
%added mask=mask(1:sz)...check to make sure this always works too

fS=6;
tic
for i=1:length(index)
    [r,c]=ind2sub(size(mask),index(i));
     if r-fS>0 && c-fS>0 && r+fS<sz(1) && c+fS<sz(2) && min(min(mask(r-fS:r+fS, c-fS:c+fS,:))) ~=0
         %r>fS && c>fS && r+fS<sz(1) && c+fS<sz(2) && mask(r-fS,c)>0 && mask(r, c-fS)>0 && mask(r+fS,c)<sz(1) && mask(r,c+fS)<sz(2)
        index2(end+1)=i;
     end
end
toc


count=0;
% meanPSF=zeros((fS*2+1*fS*2+1)*2,length(index2));
meanPSF=[];
GcampZs=[];

length(index2)

index2=datasample(index2,length(index2),'Replace',false);

tic
for i=1:length(index2)
    [r,c]=ind2sub(size(mask),index(index2(i)));
     %if r-100>0 && c-100>0 && r+100<sz(1) && c+100<sz(2)
        tmp=Rcamp(r-fS:r+fS, c-fS:c+fS,:);
        tmp=reshape(tmp, (fS*2+1)*(fS*2+1), sz(3));
        tmp2=Gcamp(r-fS:r+fS, c-fS:c+fS,:);
        tmp2=reshape(tmp2, (fS*2+1)*(fS*2+1), sz(3));
        locMean=(mean([tmp; tmp2]));
        %subtract local average, i.e. 'pre-whiten'
        %tmp=tmp-locMean;
        tmp=tmp';
        
        %for autocorrelative component, add this part:
        %remember to take out middle pixel      
        %tmp2=tmp2-locMean;
        tmp2(((fS*2+1)^2+1)/2,:)=0;
        tmp2=tmp2';
        
%         %subtract "G_0"
        GcampZ=(Gcamp(r,c,:));%-mean(tmp2);%meanG(index(index2(i))));
        GcampZ=reshape(GcampZ,1,sz(3));
        %GcampZ=GcampZ-locMean;
        GcampZ=GcampZ';
        GcampZs=vertcat(GcampZs, GcampZ);
        
        %meanPSF_G=regress(GcampZs,tmp);
        meanPSF_1=horzcat(tmp, tmp2);
        %meanPSF_1=tmp;
        meanPSF=vertcat(meanPSF,meanPSF_1);
               
%         if abs(mean(mean(meanPSF_G)))<1000
%             meanPSF(:,:,i)=meanPSF_G;
            count=count+1;
%         end
     %end
end
toc
meanPSF1=meanPSF;
GcampZs1=GcampZs;



load '01_tmpMatsS3.mat'

Gcamp=zeros(sz(1)*sz(2),sz(3));
Gcamp(index,:)=Gcamp_tmp';
Rcamp=zeros(sz(1)*sz(2),sz(3));
Rcamp(index,:)=Rcamp_tmp';
Rcamp=reshape(Rcamp,sz);
Gcamp=reshape(Gcamp,sz);
Gcamp=fliplr(Gcamp);
Rcamp=fliplr(Rcamp);

Cors=zeros(sz(1)*2-1,sz(2)*2-1,sz(3));
parfor i=1:sz(3)
Cors(:,:,i)=xcorr2(Gcamp(:,:,i),Rcamp(:,:,i));
end
%figure; imagesc(mean(Cors,3))

[r,c]=find(mean(Cors,3)==max(max(mean(Cors,3))))
xshift=r-sz(1)
yshift=c-sz(2)

newRcamp=zeros(sz);

newRcamp(xshift:end, 1:end+yshift, :)=Rcamp(1:end-xshift+1, -yshift+1:end,:);
Rcamp=newRcamp;

Gcamp=imresize(Gcamp,0.87,'nearest');
Rcamp=imresize(Rcamp,0.87,'nearest');
Gcamp=imresize(Gcamp,0.5,'nearest');
Rcamp=imresize(Rcamp,0.5,'nearest');


meanG=mean(Gcamp,3);

index=find(meanG);

%constrain to pixels that are not near the edge of the image
index2=[];

[sROI] = ReadImageJROI('RoiSet5.zip');
mask = poly2mask(sROI{1,1}.mnCoordinates(:,1), sROI{1,1}.mnCoordinates(:,2),sz(1),640)+poly2mask(sROI{1,2}.mnCoordinates(:,1), sROI{1,2}.mnCoordinates(:,2),sz(1),640);
mask = fliplr(mask); 
mask=imresize(mask,0.87,'nearest');
mask=imresize(mask,0.5,'nearest');
sz=size(Gcamp);
mask=mask(1:sz(1), 1:sz(2));

fS=6;
tic
for i=1:length(index)
    [r,c]=ind2sub(size(meanG),index(i));
     if r-fS>0 && c-fS>0 && r+fS<sz(1) && c+fS<sz(2) && min(min(mask(r-fS:r+fS, c-fS:c+fS,:))) ~=0
         %r>fS && c>fS && r+fS<sz(1) && c+fS<sz(2) && mask(r-fS,c)>0 && mask(r, c-fS)>0 && mask(r+fS,c)<sz(1) && mask(r,c+fS)<sz(2)
        index2(end+1)=i;
     end
end
toc


count=0;
% meanPSF=zeros((fS*2+1*fS*2+1)*2,length(index2));
meanPSF=[];
GcampZs=[];

length(index2)

index2=datasample(index2,2000,'Replace',false);

tic
parfor i=1:2000
    [r,c]=ind2sub(size(meanG),index(index2(i)));
     %if r-100>0 && c-100>0 && r+100<sz(1) && c+100<sz(2)
        tmp=Rcamp(r-fS:r+fS, c-fS:c+fS,:);
        tmp=reshape(tmp, (fS*2+1)*(fS*2+1), sz(3));
        tmp2=Gcamp(r-fS:r+fS, c-fS:c+fS,:);
        tmp2=reshape(tmp2, (fS*2+1)*(fS*2+1), sz(3));
        locMean=(mean([tmp; tmp2]));
        %subtract local average, i.e. 'pre-whiten'
        %tmp=tmp-locMean;
        tmp=tmp';
        
        %for autocorrelative component, add this part:
        %remember to take out middle pixel      
        %tmp2=tmp2-locMean;
        tmp2(((fS*2+1)^2+1)/2,:)=[];
        tmp2=tmp2';
        
%         %subtract "G_0"
        GcampZ=(Gcamp(r,c,:));%-mean(tmp2);%meanG(index(index2(i))));
        GcampZ=reshape(GcampZ,1,sz(3));
        %GcampZ=GcampZ-locMean;
        GcampZ=GcampZ';
        GcampZs=vertcat(GcampZs, GcampZ);
        
        %meanPSF_G=regress(GcampZs,tmp);
        meanPSF_1=horzcat(tmp, tmp2);
        %meanPSF_1=tmp;
        meanPSF=vertcat(meanPSF,meanPSF_1);
               
%         if abs(mean(mean(meanPSF_G)))<1000
%             meanPSF(:,:,i)=meanPSF_G;
            count=count+1;
%         end
     %end
end
toc

meanPSF2=meanPSF;
GcampZs2=GcampZs;





load '10_precor.mat'

Gcamp=zeros(sz(1)*sz(2),sz(3));
Gcamp(index,:)=Gcamp_tmp';
Rcamp=zeros(sz(1)*sz(2),sz(3));
Rcamp(index,:)=Rcamp_tmp';
Rcamp=reshape(Rcamp,sz);
Gcamp=reshape(Gcamp,sz);

Cors=zeros(sz(1)*2-1,sz(2)*2-1,sz(3));
parfor i=1:sz(3)
Cors(:,:,i)=xcorr2(Gcamp(:,:,i),Rcamp(:,:,i));
end
%figure; imagesc(mean(Cors,3))

[r,c]=find(mean(Cors,3)==max(max(mean(Cors,3))))
xshift=r-sz(1)
yshift=c-sz(2)

newRcamp=zeros(sz);

newRcamp(xshift:end, 1:end+yshift, :)=Rcamp(1:end-xshift+1, -yshift+1:end,:);
Rcamp=newRcamp;

Gcamp=imresize(Gcamp,0.714,'nearest');
Rcamp=imresize(Rcamp,0.714,'nearest');
Gcamp=imresize(Gcamp,0.5,'nearest');
Rcamp=imresize(Rcamp,0.5,'nearest');


meanG=mean(Gcamp,3);

index=find(meanG);

%constrain to pixels that are not near the edge of the image
index2=[];

[sROI] = ReadImageJROI('RoiSet6.zip');
mask = poly2mask(sROI{1,1}.mnCoordinates(:,1), sROI{1,1}.mnCoordinates(:,2),sz(1),640)+poly2mask(sROI{1,2}.mnCoordinates(:,1), sROI{1,2}.mnCoordinates(:,2),sz(1),640);
mask=imresize(mask,0.714,'nearest');
mask=imresize(mask,0.5,'nearest');
sz=size(Gcamp);
mask=mask(1:sz(1), 1:sz(2));

fS=6;
tic
for i=1:length(index)
    [r,c]=ind2sub(size(meanG),index(i));
     if r-fS>0 && c-fS>0 && r+fS<sz(1) && c+fS<sz(2) && min(min(mask(r-fS:r+fS, c-fS:c+fS,:))) ~=0
         %r>fS && c>fS && r+fS<sz(1) && c+fS<sz(2) && mask(r-fS,c)>0 && mask(r, c-fS)>0 && mask(r+fS,c)<sz(1) && mask(r,c+fS)<sz(2)
        index2(end+1)=i;
     end
end
toc


count=0;
% meanPSF=zeros((fS*2+1*fS*2+1)*2,length(index2));
meanPSF=[];
GcampZs=[];

length(index2)

index2=datasample(index2,2500,'Replace',false);

tic
parfor i=1:2500
    [r,c]=ind2sub(size(meanG),index(index2(i)));
     %if r-100>0 && c-100>0 && r+100<sz(1) && c+100<sz(2)
        tmp=Rcamp(r-fS:r+fS, c-fS:c+fS,:);
        tmp=reshape(tmp, (fS*2+1)*(fS*2+1), sz(3));
        tmp2=Gcamp(r-fS:r+fS, c-fS:c+fS,:);
        tmp2=reshape(tmp2, (fS*2+1)*(fS*2+1), sz(3));
        locMean=(mean([tmp; tmp2]));
        %subtract local average, i.e. 'pre-whiten'
        %tmp=tmp-locMean;
        tmp=tmp';
        
        %for autocorrelative component, add this part:
        %remember to take out middle pixel      
        %tmp2=tmp2-locMean;
        tmp2(((fS*2+1)^2+1)/2,:)=[];
        tmp2=tmp2';
        
%         %subtract "G_0"
        GcampZ=(Gcamp(r,c,:));%-mean(tmp2);%meanG(index(index2(i))));
        GcampZ=reshape(GcampZ,1,sz(3));
        %GcampZ=GcampZ-locMean;
        GcampZ=GcampZ';
        GcampZs=vertcat(GcampZs, GcampZ);
        
        %meanPSF_G=regress(GcampZs,tmp);
        meanPSF_1=horzcat(tmp, tmp2);
        %meanPSF_1=tmp;
        meanPSF=vertcat(meanPSF,meanPSF_1);
               
%         if abs(mean(mean(meanPSF_G)))<1000
%             meanPSF(:,:,i)=meanPSF_G;
            count=count+1;
%         end
     %end
end
toc

meanPSF3=meanPSF;
GcampZs3=GcampZs;


meanPSF=[meanPSF1; meanPSF2; meanPSF3];

test=ones(length(meanPSF),1);
meanPSF=horzcat(meanPSF,test);

GcampZs=[GcampZs1; GcampZs2; GcampZs3];


[meanPSF_G,bint,r,rint,stats]=regress(GcampZs,meanPSF);
%meanPSF_G=regress(GcampZs,meanPSF);

save('combinedS2S3_Ronly_250k.mat', 'meanPSF_G', 'bint', 'r', 'rint', 'stats');

sqSz=(fS*2+1)*(fS*2+1);
psf1=meanPSF_G(1:sqSz);
%psf2=vertcat(meanPSF_G(290:329),0,meanPSF_G(330:369));
psf2=vertcat(meanPSF_G((sqSz+1):(1.5*(sqSz)-.5)),0,meanPSF_G((1.5*(sqSz)-.5):length(meanPSF_G)-2));
psf1=reshape(psf1, fS*2+1,fS*2+1);
psf2=reshape(psf2, fS*2+1,fS*2+1);
figure; imagesc(psf1)
figure; imagesc(psf2)


index2=1:length(G1);
index2=datasample(index2,250000,'Replace',false);
G1=G1(index2); m1=m1(index2,:);
index2=1:length(G2);
index2=datasample(index2,250000,'Replace',false);
G2=G2(index2); m2=m2(index2,:);
index2=1:length(G3);
index2=datasample(index2,250000,'Replace',false);
G3=G3(index2); m3=m3(index2,:);
index2=1:length(G4);
index2=datasample(index2,250000,'Replace',false);
G4=G4(index2); m4=m4(index2,:);
index2=1:length(G5);
index2=datasample(index2,250000,'Replace',false);
G5=G5(index2); m5=m5(index2,:);
index2=1:length(G6);
index2=datasample(index2,250000,'Replace',false);
G6=G6(index2); m6=m6(index2,:);
index2=1:length(G7);
index2=datasample(index2,250000,'Replace',false);
G7=G7(index2); m7=m7(index2,:);
index2=1:length(G8);
index2=datasample(index2,250000,'Replace',false);
G8=G8(index2); m8=m8(index2,:);
index2=1:length(G9);
index2=datasample(index2,250000,'Replace',false);
G9=G9(index2); m9=m9(index2,:);
index2=1:length(G10);
index2=datasample(index2,250000,'Replace',false);
G10=G10(index2); m10=m10(index2,:);