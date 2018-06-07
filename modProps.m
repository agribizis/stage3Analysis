function [onsets, distC, velC, pixIEIs, pixDurs, ObjectIndices, Areas, durs, thetas, freq]=modProps(region)

Areas=[region.STATS.Area];
rho=[region.Vsum.rho];
rho=rho(1:length(Areas));
badComponents = find(Areas<1000 | rho < 500);  %if Removing only 1fr activation domains from the NoiseClusterIdx
ObjectIndices =  setdiff(1:length(region.STATS),badComponents);
region.ObjectIndices=ObjectIndices;

Areas=Areas(ObjectIndices);
durs=region.durations(ObjectIndices);
freq=length(ObjectIndices)/5;
thetas=[region.Vsum.theta];
thetas=thetas(ObjectIndices);
onsets=region.onsets(ObjectIndices);

sz=region.CC.ImageSize;
A3 = false(sz);
for i = 1:length(ObjectIndices)
    A3(region.CC.PixelIdxList{ObjectIndices(i)}) = true;
end

implay(A3)

index=find(max(A3,[],3));
A3=reshape(A3, sz(1)*sz(2), sz(3));

pixs=A3(index,:);

pixIEIs=[];
pixDurs=[];

for i=1:size(pixs,1)
    t=find(pixs(i,:));
    diffT=diff(t);
    pixIEIs=cat(2,pixIEIs,diffT(find(diffT>1)));
    onset=1;
     for j=1:length(diffT)
         if diffT(j)>1
             pixDurs=cat(2, pixDurs, j-onset);
             onset=j;
         end
     end   
end

data=region.STATS;
%data=region.domainData.STATS;

badWaves=[];

for i = 1:length(data) 
    %if data(i).Area < 1000 | region.domainData.Vsum(i).rho < 1000 
    if data(i).Area < 1000 || region.Vsum(i).rho < 500
        badWaves(end+1)=i;
    end 
end
data(badWaves) = []; 

nDomains = numel(data);
distC = zeros(nDomains,1);
velC = zeros(nDomains,1);

duration=region.durations;


for k = 1:nDomains
    tmp=data(k).PixelList;
    mi=min(tmp(:,3));
    ma=max(tmp(:,3));
    
    c = zeros(nDomains,2);
    tempDist=zeros(ma-mi,1);
    
    for i= 1:(ma-mi+1)
        [a]=find(tmp(:,3)==mi+i-1);
    
        x = tmp(a,1);
        y = tmp(a,2);
   
        q = (sum(x(1:end-1).*y(2:end)) - sum(x(2:end).*y(1:end-1)))/2;
        cx = sum((x(1:end-1)+x(2:end)).*(x(1:end-1).*y(2:end)-x(2:end).*y(1:end-1)))/(6*q);
        cy = sum((y(1:end-1)+y(2:end)).*(x(1:end-1).*y(2:end)-x(2:end).*y(1:end-1)))/(6*q);
        
        c(i,:)=[cx cy];
    end
    
    for j= 1:(ma-mi)
        centMo=c(j+1,:)-c(j,:);
        [th,r]=cart2pol(centMo(1),centMo(2));
        tempDist(j)=tempDist(j)+r; 
    end
    
 distC(k)=sum(tempDist);
 velC(k)=sum(tempDist)./duration(k); 
 
end

save(['regionVelsDisPixIEIdurs'],'onsets','distC','velC','pixDurs','pixIEIs', 'Areas', 'durs', 'freq', 'thetas', '-v7.3');

end
