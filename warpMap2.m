%this is the correct warp map to use

V1movie(isnan(V1movie))=0;
SCmovie(isnan(SCmovie))=0;
SCmask=(max(SCmovie,[],3)>0);
V1mask=(max(V1movie,[],3)>0);
indexSC=find(SCmask);
% nV1=(mat2gray(V1movie));
% nSC=(mat2gray(SCmovie));
%SCmovie=reshape(SCmovie,sz(1)*sz(2),sz(3));
%V1movie=reshape(V1movie,sz(1)*sz(2),sz(3));
% nV1=(zscore(V1movie(:)));
% nSC=(zscore(SCmovie(:)));
nSC=SCmovie;
nV1=V1movie;
nSC=reshape(nSC,sz);
nV1=reshape(nV1,sz);
%nV1=zscore(nV1(:)); nSC=zscore(nSC(:)); nV1=reshape(nV1,sz); nSC=reshape(nSC,sz);
% newV1=zeros(sz);
% % for i = 1:sz(3)
% %     newV1(:,:,i) = imwarp(V1movie(:,:,i),tform,'OutputView',imref2d([sz(1) sz(2)]));
% % end
% for i = 1:sz(3)
%     RegisteredImage = imwarp(nV1(:,:,i), movingRefObj, tform, 'OutputView', fixedRefObj, 'SmoothEdges', true);
%     newV1(:,:,i)=imwarp(RegisteredImage,D);
% end
% Cors=zeros(sz(1)*2-1,sz(2)*2-1,sz(3));
% parfor i=1:sz(3)
% Cors(:,:,i)=xcorr2(newV1(:,:,i),SCmovie(:,:,i));
% end
% %figure; imagesc(mean(Cors,3))
%
% [r,c]=find(mean(Cors,3)==max(max(mean(Cors,3))))
% xshift=r-sz(1)
% yshift=c-sz(2)
%
% V1shifted=zeros(sz);
%
% V1shifted(xshift:end, yshift:end, :)=newV1(1:end-xshift+1, 1:end-yshift+1,:);
% %V1shifted(xshift:end, 1:end+yshift, :)=newV1(1:end-xshift+1, -yshift+1:end,:);
% newV1=V1shifted;
fS=6;
index2=[];
tic
for i=1:length(indexSC)
    [r,c]=ind2sub([sz(1) sz(2)],indexSC(i));
    if r-fS>0 && c-fS>0 && r+fS<sz(1) && c+fS<sz(2) && min(min(SCmask(r-fS:r+fS, c-fS:c+fS,:))) ~=0
        %r>fS && c>fS && r+fS<sz(1) && c+fS<sz(2) && mask(r-fS,c)>0 && mask(r, c-fS)>0 && mask(r+fS,c)<sz(1) && mask(r,c+fS)<sz(2)
        index2(end+1)=i;
    end
end
toc
indexSC=indexSC(index2);
%SCpix=find(max(SCmovie,[],3)>0);
[r c]=ind2sub([sz(1) sz(2)],indexSC(:));
SCpix=[r c];
V1pix=[];
[V1pix(:,2),V1pix(:,1)]=transformPointsForward(tform,SCpix(:,2),SCpix(:,1));
V1pix=round(V1pix);
% V1ind=[];
% for i = 1:length(V1pix)
%     r=V1pix(i,1); c=V1pix(i,2);
%     if min(min(V1mask(r-fS:r+fS, c-fS:c+fS,:))) ~=0
%         V1ind(end+1)=i
%     end
% end
meanPSF=[];
GcampZs=[];
tic
for i=1:size(SCpix,1)
%[r c]=SCpix(i,:);
    tmp=nSC(SCpix(i,1)-fS:SCpix(i,1)+fS, SCpix(i,2)-fS:SCpix(i,2)+fS,:);
    tmp=reshape(tmp, (fS*2+1)*(fS*2+1), sz(3));
    %[r2 c2]=V1pix(i,:);
    tmp2=nV1(V1pix(i,1)-fS:V1pix(i,1)+fS, V1pix(i,2)-fS:V1pix(i,2)+fS,:);
    tmp2=reshape(tmp2, (fS*2+1)*(fS*2+1), sz(3));
    tmp=tmp';
    tmp2(((fS*2+1)^2+1)/2,:)=0;
    tmp2=tmp2';
    GcampZ=(nV1(V1pix(i,1),V1pix(i,2),:));%-mean(tmp2);%meanG(index(index2(i))));
    GcampZ=reshape(GcampZ,1,sz(3));
    GcampZ=GcampZ';
    GcampZs=vertcat(GcampZs, GcampZ);
    meanPSF_1=horzcat(tmp, tmp2);
    meanPSF=vertcat(meanPSF,meanPSF_1);
end
toc
test=ones(length(meanPSF),1);
meanPSF=horzcat(meanPSF,test);
[meanPSF_G,bint,r,rint,stats]=regress(GcampZs,meanPSF);
sqSz=(fS*2+1)*(fS*2+1);
psf1=meanPSF_G(1:sqSz);
psf1=reshape(psf1, fS*2+1,fS*2+1);
figure; imagesc(psf1)
psf2=meanPSF_G(sqSz+1:sqSz*2);
psf2=reshape(psf2, fS*2+1,fS*2+1);
figure; imagesc(psf2)

save('meanPSFdataMappedAffine4_s2_noBadV1.mat', 'meanPSF_G', 'meanPSF', 'GcampZs', 'bint', 'r', 'rint', 'stats', '-v7.3');