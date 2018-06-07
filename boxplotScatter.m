group = [repmat({'S2'}, length(S2), 1); repmat({'S3'}, length(S3), 1); repmat({'S3inh'}, length(S3i), 1)];

%figure; boxplot([S2S2reg S2S2gaus S2S3reg S2S3gaus S3S2reg S3S2gaus S3S3reg S3S3gaus],group)
figure; boxplot([S2' S3' S3i'],group)

%first parameter is your category groups, second parameter is a vector
%separating the groups- see matlab documentation on boxplot for more
%information on how to set this up
ydata=([S2; S3; S3i;]);
%same as first parameter above, but now set up into columns
[r, c] = size(ydata);

%you only need the bottom part if you want to add a scatter plot to your
%boxplot- if you have a lot of cells (which it looks like you do), this
%might end up looking too busy
xdata = [repmat(1, length(S2), 1); repmat(2, length(S3), 1); repmat(3, length(S3i), 1)];
hold on;
scatter(xdata(:), ydata(:), 'r.', 'jitter','on', 'jitterAmount', 0.05);
hold off;