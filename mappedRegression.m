


IMpoints=max(IMs,[],3);

BW=roipoly(IMpoints)
se = strel('disk',5);
BW2 = imdilate(IMpoints,se);
figure; imagesc(imgaussfilt(BW2,4))
BW3=imgaussfilt(BW2,4);


filename2 = 'RoiSetcropped.zip';
sz=[540 640]
[sROI] = ReadImageJROI(filename2);
SCmask=poly2mask(sROI{1,1}.mnCoordinates(:,1), sROI{1,1}.mnCoordinates(:,2),sz(1),sz(2));
V1mask=poly2mask(sROI{1,2}.mnCoordinates(:,1), sROI{1,2}.mnCoordinates(:,2),sz(1),sz(2));
bothMasks=V1mask+SCmask;
V1mask=imresize(V1mask,0.5,'nearest'); SCmask=imresize(SCmask,0.5,'nearest'); bothMasks=imresize(bothMasks,0.5,'nearest'); 
index=find(bothMasks);
%SCindex=find(SCmask);
%V1index=find(BW);
sz=size(SCmask);

V1pix=[];
SCpix=[];
V1IM=BW3.*BW;
Vvals=V1IM(V1index);
SCIM=BW3.*SCmask;
for i=1:round(length(SCindex)/2)-1
    [r c]=ind2sub([sz(1) sz(2)],(SCindex(i*2)));
    SCpix(i,:)=[r c];
    val=BW3(SCindex(i*2));
    [Vval index]=(min(abs(Vvals-val)));
    if find(V1IM==val-Vval)
        [r c]=ind2sub([sz(1) sz(2)],find(V1IM==val-Vval));
    else
        [r c]=ind2sub([sz(1) sz(2)],find(V1IM==Vval+val));
    end
    V1pix(i,:)=[r c];    
end