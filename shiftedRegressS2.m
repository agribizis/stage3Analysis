load '02_tmpMats.mat'

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

newRcamp(xshift:end, yshift:end, :)=Rcamp(1:end-xshift+1, 1:end-yshift+1,:);
%newRcamp(xshift:end, 1:end+yshift, :)=Rcamp(1:end-xshift+1, -yshift+1:end,:);
Rcamp=newRcamp;

Gcamp=imresize(Gcamp,.5,'nearest');
Rcamp=imresize(Rcamp,.5,'nearest');


meanG=mean(Gcamp,3);

index=find(meanG);

%constrain to pixels that are not near the edge of the image
index2=[];

[sROI] = ReadImageJROI('RoiSet1.zip');
mask = poly2mask(sROI{1,1}.mnCoordinates(:,1), sROI{1,1}.mnCoordinates(:,2),sz(1),640)+poly2mask(sROI{1,2}.mnCoordinates(:,1), sROI{1,2}.mnCoordinates(:,2),sz(1),640);
%mask = poly2mask(sROI{1,1}.mnCoordinates(:,1), sROI{1,1}.mnCoordinates(:,2)-270,sz(1),640)+poly2mask(sROI{1,2}.mnCoordinates(:,1), sROI{1,2}.mnCoordinates(:,2)-270,sz(1),640);
%changed sz(2) to 640, subtracted 270 if ROI is on the bottom...check this
%on next round
%mask = fliplr(mask); 
mask=imresize(mask,0.5,'nearest');
sz=size(Gcamp);
mask=mask(1:sz(1), 1:sz(2));
%added mask=mask(1:sz)...check to make sure this always works too

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

index2=datasample(index2,4000,'Replace',false);

tic
parfor i=1:4000
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
meanPSF1=meanPSF;
GcampZs1=GcampZs;



load '04_tempMat.mat'

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

Gcamp=imresize(Gcamp,0.83,'nearest');
Rcamp=imresize(Rcamp,0.83,'nearest');
Gcamp=imresize(Gcamp,.5,'nearest');
Rcamp=imresize(Rcamp,.5,'nearest');

meanG=mean(Gcamp,3);

index=find(meanG);

%constrain to pixels that are not near the edge of the image
index2=[];

[sROI] = ReadImageJROI('RoiSet2.zip');
%mask = poly2mask(sROI{1,1}.mnCoordinates(:,1), sROI{1,1}.mnCoordinates(:,2),sz(1),640)+poly2mask(sROI{1,2}.mnCoordinates(:,1), sROI{1,2}.mnCoordinates(:,2),sz(1),640);
mask = poly2mask(sROI{1,1}.mnCoordinates(:,1), sROI{1,1}.mnCoordinates(:,2)-270,sz(1),640)+poly2mask(sROI{1,2}.mnCoordinates(:,1), sROI{1,2}.mnCoordinates(:,2)-270,sz(1),640);
mask = fliplr(mask); 
mask=imresize(mask,0.83,'nearest');
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

index2=datasample(index2,3000,'Replace',false);

tic
parfor i=1:3000
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





load '05_tmpMats.mat'

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

Gcamp=imresize(Gcamp,.5,'nearest');
Rcamp=imresize(Rcamp,.5,'nearest');

meanG=mean(Gcamp,3);

index=find(meanG);

%constrain to pixels that are not near the edge of the image
index2=[];

[sROI] = ReadImageJROI('RoiSet3.zip');
mask = poly2mask(sROI{1,1}.mnCoordinates(:,1), sROI{1,1}.mnCoordinates(:,2),sz(1),640)+poly2mask(sROI{1,2}.mnCoordinates(:,1), sROI{1,2}.mnCoordinates(:,2),sz(1),640);
%mask = poly2mask(sROI{1,1}.mnCoordinates(:,1), sROI{1,1}.mnCoordinates(:,2)-270,sz(1),640)+poly2mask(sROI{1,2}.mnCoordinates(:,1), sROI{1,2}.mnCoordinates(:,2)-270,sz(1),640);
mask = fliplr(mask); 
mask = imresize(mask,0.5,'nearest');
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

index2=datasample(index2,5000,'Replace',false);

tic
parfor i=1:5000
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

save('meanPSF_s206_fS6.mat', 'meanPSF', 'GcampZs','-v7.3');

sqSz=(fS*2+1)*(fS*2+1);
psf1=meanPSF_G(1:sqSz);
%psf2=vertcat(meanPSF_G(290:329),0,meanPSF_G(330:369));
psf2=vertcat(meanPSF_G((sqSz+1):(1.5*(sqSz)-.5)),0,meanPSF_G((1.5*(sqSz)-.5):length(meanPSF_G)-2));
psf1=reshape(psf1, fS*2+1,fS*2+1);
psf2=reshape(psf2, fS*2+1,fS*2+1);
figure; imagesc(psf1)
figure; imagesc(psf2)