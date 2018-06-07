function [distC, velC]=distanceCummulative(region)

data=region.STATS;
%data=region.domainData.STATS;

badWaves=[];

for i = 1:length(data) 
    %if data(i).Area < 1000 | region.domainData.Vsum(i).rho < 1000 
    if data(i).Area < 1000 %| region.Vsum(i).rho < 1000 
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
end