%used for thalamus analysis
function preProcessingAG_5 (filename1, filename2)
%pre-processing: create files.txt with names of .tif files
filename1 = input('filename.tif = ');
%filename2 = input('filename2.tif = ');
filename2 = 'RoiSet.zip';
A = openMovie(filename1);
A = double(A);
sz = size(A);


[sROI] = ReadImageJROI(filename2);
sz=[540 640 sz(3)];
mask1=poly2mask(sROI{1,1}.mnCoordinates(:,1), sROI{1,1}.mnCoordinates(:,2),sz(1),sz(2));
mask2=poly2mask(sROI{1,2}.mnCoordinates(:,1), sROI{1,2}.mnCoordinates(:,2),sz(1),sz(2));
%mask1=imresize(mask1, 0.5, 'nearest');
%mask2=imresize(mask2, 0.5, 'nearest');
sz = size(A);
%ROI = ReadImageJROI('RoiSet.zip');       
%bothMasks = poly2mask(ROI{1,1}.mnCoordinates(:,1), ROI{1,1}.mnCoordinates(:,2),sz(1),sz(2))+poly2mask(ROI{1,2}.mnCoordinates(:,1), ROI{1,2}.mnCoordinates(:,2),sz(1),sz(2));
bothMasks=mask1+mask2;
index=find(bothMasks);
figure; imagesc(mask1)
mask= input('mask name =');


%%
region.badPCs=badPCs;
A=mov;
A=reshape(A,sz(1)*sz(2),sz(3));
%index=1:length(A); %REMOVE THIS LATER
    

hat = 150;
se = strel('line', hat, 0);
flpA=horzcat(fliplr(A(:,1:hat)),A(:,:));
A_ln=zeros(sz(1)*sz(2), (sz(3))+hat);

for i = 1:length(index)
A_ln(index(i),:)=imtophat(flpA(index(i),:), se);
end
A=A_ln(:,(hat+1):end); %A=reshape(A,sz);
%sigma=input('sigma = ');  %3px gauss smooth

clear A_ln flpA
A = A ./ (mean(A,2) * ones(1,sz(3))) - 1;

A=reshape(A,sz(1)*sz(2),sz(3));
C=medfilt2(A);


%%

sigma=3;        
%thresh=input('thresh = ');  %threshold in no. of std dev
thresh=2;
B = zscore(C'); %covert to standardized zscores so that mean=0, and sd = 1;
%clear C 

%mask=1;
B=reshape(B',sz).*mask;
B(B<thresh)=0;
%V1Mov=B.*V1mask;
%SCMov=B(B>3).*SCmask;

A = reshape(A,sz);
%A2=false(size(A));
Iarr=zeros(size(A));
bothMasks3D = repmat(bothMasks,[1 1 sz(3)]);

tic
for fr = 1:sz(3); %option:parfor, do NOT use parfor it is 3x slower
	I = B(:,:,fr);
		I2 = gaussSmooth(I,sigma,'same');
	Iarr(:,:,fr) = I2;
end	
toc


bw = Iarr > thresh;
A2 = bothMasks3D & bw;
region.graythresh = thresh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==2==Detection============================
%[A3, CC, STATS] = wholeBrain_detectAG(A2,A,[4 3],[],fnm,region,hemisphereIndices);
CC = bwconncomp(A2);
%-----Get 3D connected components structure and some STATS available for n-D arrays-----
STATS = regionprops(CC,A,'Area','BoundingBox', 'Centroid', 'MaxIntensity', 'MinIntensity', 'MeanIntensity','MinorAxisLength');  %some of the properties in regionprops that work on n-D arrays
%Get measurements that will be used as inputs to kmeans---
% centr = vertcat(STATS.Centroid);
% centrZ = round(centr(:,3));
% roiArea=[STATS.Area];  %Scalar; the actual number of pixels in the region. (This value might differ slightly from the value returned by bwarea, which weights different patterns of pixels differently.
% roiMean=[STATS.MeanIntensity];
% roiMax=[STATS.MaxIntensity];
roiBoundingBox = zeros(length(STATS),6);
for i = 1:length(STATS)
    roiBoundingBox(i,:) = STATS(i).BoundingBox;
end
durations = roiBoundingBox(:,6);  %add to kmeans
diameters = mean([roiBoundingBox(:,4) roiBoundingBox(:,5)], 2);
badComponents = find(durations<2);  %if Removing only 1fr activation domains from the NoiseClusterIdx
ObjectIndices =  setdiff(1:length(durations),badComponents);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----Remake CC data structure with desired objects (deleting noise components)----
%now we remake the connected components (CC) data structure by keeping only the objects in the cluster of interest (the functional signal cluster)
newPixelIdxList = {};
count = 0;
tic
for i = 1:length(ObjectIndices)
    count = count+1;
    newPixelIdxList{count} = CC.PixelIdxList{ObjectIndices(i)};
end
toc
%CCorig=CC;
CC.PixelIdxList = newPixelIdxList;
CC.NumObjects = length(CC.PixelIdxList);
STATS = regionprops(CC,A,'Area','BoundingBox', 'Centroid', 'MaxIntensity', 'MinIntensity', 'MeanIntensity', 'Image', 'PixelIdxList', 'PixelList', 'SubarrayIdx'); %all the properties in regionprops that work on n-D arrays
%Now make a binary movie array based on the segmented functional signal domains
A3 = false(sz);
for i = 1:CC.NumObjects
    A3(CC.PixelIdxList{i}) = true;
end

% for i = 1:CC.NumObjects
%     if mask1(STATS(i).PixelList(1,2),STATS(i).PixelList(1,1))==1
%         location(i)=1;
%     else
%         location(i)=2;
%     end
% end
% region.location=location; 


M(size(A3,3)) = struct('cdata',[],'colormap',[]);
for fr=1:size(A3,3)
    [I2, map] = gray2ind(A3(:,:,fr), 8); %figure; imshow(I2,map)
    M(fr) = im2frame(I2,map);
end
writeMovie(M,['Movie' datestr(now,'yyyymmdd-HHMMSS') '.avi']);

% fnm2 = [fnm(1:length(fnm)-4) '_' datestr(now,'yyyymmdd-HHMMSS') '.mat'];
% save([fnm2(1:length(fnm2)-4) '_connComponents_BkgndSubtr60' '.mat'],'A2','A3','CC','STATS','-v7.3');  
% clear A2 A3;




% %==3==Format domain data structures============================
% domains = DomainSegmentationAssignmentAG(CC,STATS, 'false');  
% region.domainData.domains = domains;      
region.CC = CC;      
region.STATS = STATS;  

for i = 1:length(region.STATS)
    region.durations(i) = region.STATS(i).BoundingBox(6);
end
 
for i = 1:length(region.STATS)
    region.onsets(i)=region.STATS(i).PixelList(1,3);
region.offsets(i) = region.onsets(i)+region.STATS(i).BoundingBox(6)-1;
end

sig = [1.2];
onsets=region.onsets; offsets=region.offsets;
AVx=zeros(sz);
AVy=zeros(sz);
winSig = [3];
for fr = 1:sz(3)-1; %option:parfor
	img1 = A3(:,:,fr);
	img2 = A3(:,:,fr+1);
	[Vx,Vy,reliab] = optFlowLk( img1, img2, [], winSig, sig, 3e-6);
	AVx(:,:,fr) = Vx;
	AVy(:,:,fr) = Vy;
end
nDomains = numel(region.STATS);
Vsum(1:nDomains) = struct('theta',[],'rho',[]);
for k = 1:nDomains
	%if region.location(k) == 1
       [theta, rho]= cart2pol(-1*sum(AVx(region.STATS(k).PixelIdxList)), sum(AVy(region.STATS(k).PixelIdxList))); %this should reverse direction for the right hemisphere such that 0 degrees is toward ipsilateral (medial) for both eyes
    %else
        [theta, rho]= cart2pol(sum(AVx(region.STATS(k).PixelIdxList)), sum(AVy(region.STATS(k).PixelIdxList)));
    %end
    Vsum(k).theta = theta.*-1;
	Vsum(k).rho = rho;
end
mx = max([Vsum.rho]);

region.Vsum=Vsum;

clear CC STATS;




% % clims=[0 1];
% % arrowLength=200;
% % M(size(A3,3)) = struct('cdata',[],'colormap',[]);
% %     h=figure;
% %     imgAx=imagesc(A(:,:,1)); colormap(gray)
% %     for fr = 1:sz(3)-1
% % 
% %         delete(findobj(gca,'Type','line'))
% %         set(imgAx,'CData',A(:,:,fr)); %set image data directly instead of making new imagesc call
% % 
% %         for k = 1:nDomains
% %             if fr >= onsets(k) & fr <= offsets(k)
% %                 theta=Vsum(k).theta;
% %                 rho=Vsum(k).rho;
% %                 xpts = region.STATS(k).Centroid(1);
% %                 ypts = region.STATS(k).Centroid(2);
% %                 % [dxpts,dypts] = pol2cart(theta.* -1,rho);
% %                 [dxpts,dypts] = pol2cart(theta.*-1,(rho./mx).*arrowLength);
% %                 % amquiver(X, Y, Vx, Vy);
% %                 %arrow([xpts, ypts], [xpts + dxpts, ypts + dypts], 10, 15, 20, 1, 'EdgeColor', 'b', 'FaceColor', 'b', 'LineWidth', 2);
% %                 line(xpts,ypts,'LineStyle','none','Marker','o','MarkerSize',8,'Color','b')
% %                 line([xpts xpts+dxpts],[ypts ypts+dypts],'LineWidth',2,'Color','b')
% %             end
% %         end
% % 
% %         cdata=print(h,'-RGBImage'); %TODO: remove zbuffer dependency in these two lines and replace with 
% %         M(fr) = im2frame(cdata);  %grab figure data for frame without getframe screen draw issues
% %     end
% %     
% % fnm2=['opticflow' datestr(now,'yyyymmdd-HHMMSS') '.avi']
% % writeMovie(M,fnm2);





%A3 = reshape(A3,prod(sz(1:2)),sz(3)); %space x time matrix
clear A2 A3 AVx AVy B

save(['regionDataSVD2wVsum' datestr(now,'yyyymmdd-HHMMSS')],'region','-v7.3');     
clear region
clear Vsum durations onsets offsets
close all
end
