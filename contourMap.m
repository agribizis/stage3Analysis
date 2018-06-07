Movie=reshape(Movie,sz(1)*sz(2),sz(3));
Movie = Movie ./ (mean(Movie,2) * ones(1,sz(3))) - 1;
Movie=reshape(Movie,sz);
%why does this work with the DF/F movie but not the bleach corrected Aorig
%movie??

newMovie=Movie(:,:,[110:160]);% 540:720 1870:1971 2680:2800 4210:4380 4190:4401]); 
%1620:1710 2920:3000 3500:3590 5000:5080 6210:6320 7740:7850 8350:8450
newMovie=newMovie.*bothMasks;
newMovie=imgaussfilt(newMovie,2);

%v1Map=zeros(sz(1),sz(2));
%scMap=zeros(sz(1),sz(2));
se = strel('disk',5);
imageMap=[];
for i=1:5:size(newMovie,3)-20
    IM=newMovie(:,:,i);
    %SCIM=IM.*SCmask; V1IM=IM.*V1mask;
    IM(IM>prctile(IM(:),99))=1; IM(IM<prctile(IM(:),99))=0;
    %IM = bwmorph(IM, 'close');
    %IM = bwmorph(IM, 'open');
    %IM = imdilate(IM,se);
    %IM = bwmorph(IM, 'thin',Inf);
    %%IM = edge(IM,'sobel');
    IM=IM.*i;
    %imageMap=imageMap+IM;
    imageMap(:,:,i)=IM;
end
figure; imagesc(max(imageMap,[],3))
%figure; im=imcontour(max(imageMap,[],3))
%figure; imagesc(imageMap)

imMap=max(imageMap,[],3);
imMap=imMap.*V1mask;
imMap=imrotate(imMap,180);
figure; imagesc(imMap)
Jregistered = imwarp(imMap,tform,'OutputView',imref2d(size(moving)));
figure; imagesc(Jregistered)

% [maxProj, Iarr] = timeColorMapProj(newMovie,1, size(newMovie,3), [1],[1 4]);
% maxProj=rgb2gray(maxProj);
% figure; im=imcontour(maxProj)


moving=max(rotV1,[],3); fixed=max(SCIM,[],3); cpselect(moving,fixed)

%this is no longer necessary
%
% rotV1=[];
% for i = 1:size(IMs,3)
%     rotV1(:,:,i)=imrotate(moving(:,:,i),90);
%     rotV1(:,:,i)=flipud(rotV1(:,:,i));    
% end
% shiftedV1 = zeros(size(rotV1));
% for i = 1:size(IMs,3)  
%     shiftedV1(31:80, 20:50,i) = rotV1(1:50,25:55,i);
% end


V1pix=[];
SCpix=[];
V1IM=IMs.*V1mask;
SCIM=IMs.*SCmask;
for i=1:length(corx)
%     col=90+i;
%     row=78+i;
%     SCpix(i,:)=[90 col];
    %SCpix(i,:)=[row 45];
    %V1pix(i,:)=[row 112];
    [r c]=ind2sub([sz(1) sz(2)],find(shiftedV1(:,:,i)));
    V1pix(i,:)=[r c];
    [r c]=ind2sub([sz(1) sz(2)],find(SCIM(:,:,i)));
    SCpix(i,:)=[r c];    
end
