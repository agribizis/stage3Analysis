
npix = sz(1)*sz(2);

A = reshape(A, npix, sz(3));
% Create covariance matrix
    % Perform SVD on temporal covariance
    c1 = (A'*A)/npix; movtm = mean(A,1); % Average over space
    covmat = c1 - movtm'*movtm; clear c1
    covtrace = trace(covmat) / npix;
    %[mixedsig, CovEvals, percentvar] = cellsort_svd(covmat, nPCs, sz(3), npix);
    

%------------
% Save the output data
% save(fnmat,'mixedfilters','CovEvals','mixedsig','movm','movtm','covtrace')
% fprintf(' CellsortPCA: saving data and exiting; ')

%-----------------------
% Perform SVD
nPCs=300;
opts.disp = 0;
opts.issym = 'true';
[mixedsig, CovEvals] = eigs(covmat, nPCs, 'LM', opts);  % pca_mixedsig are the temporal signals, mixedsig
CovEvals = diag(CovEvals);
% if nnz(CovEvals<=0)
%     nPCs = nPCs - nnz(CovEvals<=0);
%     fprintf(['Throwing out ',num2str(nnz(CovEvals<0)),' negative eigenvalues; new # of PCs = ',num2str(nPCs),'. \n']);
%     mixedsig = mixedsig(:,CovEvals>0);
%     CovEvals = CovEvals(CovEvals>0);
% end

mixedsig = mixedsig' * sz(3);
CovEvals = CovEvals / npix;

percentvar = 100*sum(CovEvals)/covtrace;
fprintf([' First ',num2str(nPCs),' PCs contain ',num2str(percentvar,3),'%% of the variance.\n'])


% Re-load movie data
Sinv = inv(diag(CovEvals.^(1/2)));

movtm = mean(A,1); % Average over space
movuse = A - ones(npix,1) * movtm;
mixedfilters = reshape(movuse * mixedsig' * Sinv, npix, nPCs);

mixedfilters = reshape(mixedfilters, sz(1),sz(2),nPCs);

implay(mixedfilters);

badPCs = input('bad PCs= (example "[1:2]")');  %***change these values***
szXY = sz(1:2); szZ = size(mixedsig,2);
PCuse=setdiff([1:300],badPCs);  %***change these values***
mixedfilters2 = reshape(mixedfilters(:,:,PCuse),npix,length(PCuse));  
mov = mixedfilters2 * diag(CovEvals(PCuse).^(1/2)) * mixedsig(PCuse,:);  
mov = zscore(reshape(mov,npix*szZ,1));
mov = reshape(mov, szXY(1), szXY(2), szZ); 