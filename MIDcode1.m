%this is an attempt to preform maximally informative dimension (MID) analysis on our RF data

%part 1: load presented stimuli
%

filelist = readtext('files.txt',' '); %load file list

tmat=double(rgb2gray(imread(filelist{1})));
tmat=imresize(tmat, 0.25, 'nearest');
sz=size(tmat);
tmat=reshape(tmat,sz(1)*sz(2),1);
tmat=repmat(tmat,1,8);

A(:,:,1)=tmat;

sz(3)=2000;
for i = 2:sz(3)
    tmat=double(rgb2gray(imread(filelist{i})));
    tmat=imresize(tmat, 0.25, 'nearest');
    tmat=reshape(tmat,sz(1)*sz(2),1);
    tmat=repmat(tmat,1,8);
    A(:,:,i)=tmat;
end

%part 2: load calcium trace and multiply stimulus set by calcium trace
load('roiFrames_dF_one.mat');
load('roiFrames_dF_two.mat');
load('roiFrames_dF_three.mat');
load('roiFrames_dF_four.mat');
load('roiFrames_dF_five.mat');
load('roiFrames_dF_six.mat');
load('roiFrames_dF_seven.mat');
load('roiFrames_dF_eight.mat');
load('roiFrames_dF_nine.mat');
load('roiFrames_dF_ten.mat');

roiFrames=cat(3, roiFrames_dF_one, roiFrames_dF_two, roiFrames_dF_three, roiFrames_dF_four, roiFrames_dF_five, roiFrames_dF_six, roiFrames_dF_seven, roiFrames_dF_eight, roiFrames_dF_nine, roiFrames_dF_ten); 
sz=size(roiFrames);
roiFrames=reshape(roiFrames, sz(1), sz(2)*sz(3));
index=find(max(roiFrames,[],2)>.002);

roiFrames=reshape(roiFrames,sz);
topRoiFr=squeeze(roiFrames(index,:,:));
topRoiFr=reshape(topRoiFr,length(index),8*sz(3));
%roiFrames5=squeeze(roiFrames(index,5,:));

sz=size(A);
A=reshape(A,sz(1)*sz(2), sz(3));
tmat = repmat(topRoiFr(2,:),length(A),1);
A1 = tmat .* A;
%part 3: perform SVD, from wholeBrainAnalysis suite