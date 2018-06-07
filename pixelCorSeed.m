Movie = cat(3, A, B, C, D);
clear A B C D
Movie = double(Movie);
sz = size(Movie);
A=Movie; clear Movie
sz=[540 640 sz(3)];
[sROI] = ReadImageJROI(filename2);
V1mask=poly2mask(sROI{1,1}.mnCoordinates(:,1), sROI{1,1}.mnCoordinates(:,2),sz(1),sz(2));
SCmask=poly2mask(sROI{1,2}.mnCoordinates(:,1), sROI{1,2}.mnCoordinates(:,2),sz(1),sz(2));
SCmask=imresize(SCmask,0.25,'nearest');
V1mask=imresize(V1mask,0.25,'nearest');

SCmovie=A.*SCmask;
V1movie=A.*V1mask;

indexSC=find(SCmask);
indexV1=find(V1mask);

sz=size(A);

SCmovie=reshape(SCmovie,sz(1)*sz(2), sz(3));
V1movie=reshape(V1movie,sz(1)*sz(2), sz(3));
 
hat = 150;
se = strel('line', hat, 0);
flpA=horzcat(fliplr(SCmovie(:,1:hat)),SCmovie(:,:));
A_ln=zeros(sz(1)*sz(2), (sz(3))+hat);

for i = 1:length(indexSC)
A_ln(indexSC(i),:)=imtophat(flpA(indexSC(i),:), se);
end
SCmovie=A_ln(:,(hat+1):end); %A=reshape(A,sz);
%sigma=input('sigma = ');  %3px gauss smooth
clear A_ln flpA

flpA=horzcat(fliplr(V1movie(:,1:hat)),V1movie(:,:));
A_ln=zeros(sz(1)*sz(2), (sz(3))+hat);

for i = 1:length(indexV1)
A_ln(indexV1(i),:)=imtophat(flpA(indexV1(i),:), se);
end
V1movie=A_ln(:,(hat+1):end); %A=reshape(A,sz);
%sigma=input('sigma = ');  %3px gauss smooth
clear A_ln flpA

SCmovie = SCmovie ./ (mean(SCmovie,2) * ones(1,sz(3))) - 1;
V1movie = V1movie ./ (mean(V1movie,2) * ones(1,sz(3))) - 1;

% for i=1:200:length(indexSC)
%     %tempMov=vertcat(idxSCMov,idxSCMov(i,:));
%     MovieCoef=[];
%     for j=1:length(indexV1)
%         MovieCoef(j)=min(min(corrcoef(V1movie(indexV1(j),:),SCmovie(indexSC(i),:))));
%     end
%     %MovieCoef=corrcoef(tempMov');
%     %why do these all look the same regardless of which round it is?
%     %IM=MovieCoef(1,2:end);
%     IM=MovieCoef';
%     max(IM(:))
%     imagesCor=zeros(sz(1),sz(2));
%     imagesCor(indexV1)=IM;
%     figure; imagesc(imagesCor)
% end

IMs=zeros(sz(1),sz(2),round(length(indexSC)/50));
corx=[];
A=reshape(A, [sz(1)*sz(2) sz(3)]);
for i=5:50:length(indexSC)
    pix=squeeze(A(indexSC(i),:));
    MovieCoef=zeros(1,size(A,1));
    for j=1:length(indexV1)
        MovieCoef(indexV1(j))=min(min(corrcoef(A(indexV1(j),:),pix)));
    end
    IM=MovieCoef';
    IM=reshape(IM, [sz(1) sz(2)]);
    %figure; imagesc(IM)
    IM(IM==1)=0;
    IM=IM.*V1mask; corx(i)=max(IM(:));
    [r c]=ind2sub([sz(1) sz(2)],find(IM==max(IM(:))));
    [pix]=([round(mean(r)) round(mean(c))]);
    IM=zeros(sz(1), sz(2)); IM(pix(1), pix(2))=1;
    %IM(IM>prctile(IM(:),99))=1; IM(IM<prctile(IM(:),99))=0; 
    %comment above line out when resetting
    IMs(:,:,i)=IM*i;
end

figure; imagesc(max(IMs,[],3).*(SCmask+V1mask))
