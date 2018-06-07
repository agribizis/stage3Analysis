function [R2, er, error] = modelTest_gaussian(meanPSF_G, Gcamp, Rcamp)

fS=6;

sqSz=(fS*2+1)*(fS*2+1);
psf1=meanPSF_G(1:sqSz);
psf1=psf1;
%psf1=reshape(psf1, fS*2+1,fS*2+1);
%figure; imagesc(psf1)

psf2=meanPSF_G(sqSz+1:end);
psf2=psf2;
%psf2=reshape(psf2, fS*2+1,fS*2+1);
%figure; imagesc(psf2)

% load '01_s3_prePro_noBleach.mat'
% 
% Gcamp=zeros(sz(1)*sz(2),sz(3));
% Gcamp(index,:)=Gcamp_tmp';
% Rcamp=zeros(sz(1)*sz(2),sz(3));
% Rcamp(index,:)=Rcamp_tmp';
% Rcamp=reshape(Rcamp,sz);
% Gcamp=reshape(Gcamp,sz);
% Gcamp=fliplr(Gcamp);
% Rcamp=fliplr(Rcamp);
% Rcamp=Rcamp(:,:,1:200);
% Gcamp=Gcamp(:,:,1:200);

% xshift = 2
% yshift = -3
% newRcamp=zeros(sz);
% %newRcamp(xshift+1:end,yshift+1:end,:)=Rcamp(1:end-xshift,1:end-yshift,:);
% newRcamp(xshift+1:end,1:end+yshift,:)=Rcamp(1:end-xshift,-yshift+1:end,:);
% Rcamp=newRcamp;
% Gcamp=imresize(Gcamp,0.87,'nearest');
% Rcamp=imresize(Rcamp,0.87,'nearest');
% Gcamp=imresize(Gcamp,0.5,'nearest');
% Rcamp=imresize(Rcamp,0.5,'nearest');

% Rcamp=Rcamp(:,:,[220:410 900:1100]);
% Gcamp=Gcamp(:,:,[220:410 900:1100]);
% Rcamp=Rcamp(:,:,[220:410 900:1100 1440:1500]);
% Gcamp=Gcamp(:,:,[220:410 900:1100 1440:1500]);
% 


sz=size(Rcamp);

%psf1=reshape(psf1,1,(fS*2+1)*(fS*2+1));
%psf2=reshape(psf2,1,(fS*2+1)*(fS*2+1));

meanG=mean(Gcamp,3);
index=find(meanG);

G_pred=zeros(length(index),sz(3));

tic
parfor i=1:length(index)
    [r,c]= ind2sub([sz(1),sz(2)],index(i));
    if r-fS>0 && c-fS>0 && r+fS<sz(1) && c+fS<sz(2)        
        tmpR=Rcamp(r-fS:r+fS, c-fS:c+fS,:);
        tmpR=reshape(tmpR, (fS*2+1)*(fS*2+1),sz(3));
        tmpG=Gcamp(r-fS:r+fS, c-fS:c+fS,:);
        tmpG=reshape(tmpG, (fS*2+1)*(fS*2+1),sz(3));
        %G_pred(i,:)=myNNF5(tmpR');
        %for j=1:sz(3)
%             G_pred(i,j)=myNNF(tmpR(:,j));
        if length(meanPSF_G) > length(psf1)+1
            G_pred(i,:)=meanPSF_G(end)+psf1*tmpR+psf2*tmpG;
        else
            G_pred(i,:)=meanPSF_G(end)+psf1*tmpR;%+psf2*tmpG;
        end
    end
    i
end
toc



G_predicted=zeros(sz(1)*sz(2),sz(3));
G_predicted(index,:)=G_pred;

G_predicted=reshape(G_predicted,sz);

R = corrcoef(Gcamp,G_predicted);
R2=R.^2

err=Gcamp-G_predicted;
er=sum(err(:))
err2=err.^2;
error=sum(err2(:))

% % Rcamp=Rcamp./max(Rcamp(:));
% % Gcamp=Gcamp./max(Gcamp(:));
% % G_predicted=G_predicted./max(G_predicted(:));
% implay(G_predicted)
% Res=Gcamp-G_predicted;
% %residual=sum(Res(:))
% implay(Res)
% 
% meanRes=mean(Res,3);
% figure; imagesc(meanRes)
% maxRes=max(Res,[],3);
% figure; imagesc(maxRes)
% minRes=min(Res,[],3);
% figure; imagesc(minRes)
% 
% % [Iarr, ~] = gray2ind(G_predicted, 256);
% % Iarr2avi(Iarr,1,sz(3),'S2_04_predicted.avi');