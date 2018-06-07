Aorig=A;
V1mask=imresize(V1mask,0.25,'nearest'); SCmask=imresize(SCmask,0.25,'nearest'); bothMasks=imresize(bothMasks,0.25,'nearest');
indexSC=find(SCmask);
indexV1=find(V1mask);
A(isnan(A))=0;
A = imresize(A,0.25,'bilinear');
%A=A(:,:,[1620:1710 3500:3590 5000:5080 6210:6320 7740:7850 8350:8450]);
%A=A(:,:,[1850:1990 4640:4800 8240:8430 9550:9690]);
sz=size(A);
A=imgaussfilt(A,.5);
A=reshape(A,sz(1)*sz(2), sz(3));
IMs=zeros(sz(1),sz(2),length(indexSC));
corx=[];
iSC=indexSC;
%iSC=indexSC(1:round(length(indexSC)/20):end);
%iSC=datasample(indexSC,50,'Replace',false);
tic
for i=1:length(iSC)
%     Iarr=reshape(Iarr,sz);
%     row=78+i;
%     %col=314+i*2;
%     %col=187+i;
%     col=90+i;
%     %145+i*2;
%     pix=squeeze(Iarr(90,col,:));
%     %pix=squeeze(Iarr(row,45,:));
%     %256
%     Iarr=reshape(Iarr, [sz(1)*sz(2) sz(3)]);
pix=A(iSC(i),:);
MovieCoef=zeros(sz(1)*sz(2),1);
for j=1:length(indexV1)
MovieCoef(indexV1(j))=min(min(corrcoef(A(indexV1(j),:),pix)));
end
IM=MovieCoef';
IM=reshape(IM, [sz(1) sz(2)]);
%figure; imagesc(IM)
IM(IM==1)=0;
IM=IM.*V1mask; corx(i)=max(IM(:));
%     [r c]=ind2sub([sz(1) sz(2)],find(IM==max(IM(:))));
%     [pix]=([round(mean(r)) round(mean(c))]);
%     IM=zeros(sz(1), sz(2)); IM(pix(1), pix(2))=1;
IM(iSC(i))=1;
IM(IM>prctile(IM(:),99.8))=1; IM(IM<prctile(IM(:),99.8))=0;
%comment above line out when resetting
%figure; imagesc(IM)
IMs(:,:,i)=IM*i;
end
toc
figure; imagesc(max(IMs,[],3))
mapIM=max(IMs,[],3);
mapIM = imresize(mapIM,4,'nearest');
sz=[540 640 sz(3)]
[sROI] = ReadImageJROI(filename2);
SCmask=poly2mask(sROI{1,1}.mnCoordinates(:,1), sROI{1,1}.mnCoordinates(:,2),sz(1),sz(2));
V1mask=poly2mask(sROI{1,2}.mnCoordinates(:,1), sROI{1,2}.mnCoordinates(:,2),sz(1),sz(2));
bothMasks=V1mask+SCmask;
index=find(bothMasks);
V1mask=imresize(V1mask,0.5,'nearest'); SCmask=imresize(SCmask,0.5,'nearest'); bothMasks=imresize(bothMasks,0.5,'nearest');
moving=mapIM(1:270,:).*V1mask; fixed=mapIM(1:270,:).*SCmask;
index=find(bothMasks);
%moving=mapIM(1:270,:).*V1mask; fixed=mapIM(1:270,:).*SCmask;
moving=mapIM(1:270,:).*SCmask; fixed=mapIM(1:270,:).*V1mask;
cpselect(moving,fixed)
figure; imagesc(moving+fixed)

tform = fitgeotrans(movingPoints,fixedPoints,'affine')
Jregistered = imwarp(moving,tform,'OutputView',imref2d(size(moving)));
figure; imagesc(Jregistered-fixed)
A=Aorig;
sz=size(A);

save('tformAffine4_handPicked.mat', 'tform', 'movingPoints', 'fixedPoints', 'sz');

%use this to look at warped movie
V1movie(isnan(V1movie))=0;
SCmovie(isnan(SCmovie))=0;

newV1=zeros(sz);
for i = 1:sz(3)
    newV1(:,:,i) = imwarp(SCmovie(:,:,i),tform,'OutputView',imref2d([sz(1) sz(2)]));
end
