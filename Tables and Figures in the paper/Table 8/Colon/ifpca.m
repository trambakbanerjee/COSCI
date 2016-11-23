function [label, stats, numselect,signals,ppp] = ifpca(Data, K, KSvalue, pvalcut, rep, kmeansrep, per)
%The function IPPCA gives an estimation of cluster labels with IF-PCA
%method according to Jin and Wang (2014).
%
%Function: [label, stats, numselect] = ifpca(Data, K, KSvalue, pvalcut, rep, kmeansrep, per)
%
%Inputs: 
%Data: p by n data matrix where p is the number of features and n is the 
%number of samples. Each column presents the observations for a 
%sample. 
%K: number of clusters
%KSvalue: simulated Kolmogorov-Smirnov statistics if possible, used to 
%estimate p-value for each feature. If left null, corresponding statistics
%will be genearated by the algorithm
%pvalcut: the threshold to elminate the effect of features as outliers.
%The default value is log(p)/p.
%rep: the number of Kolmogorov-Smirnov statistics to be simulated in the
%algorithm, defacult 100*p.
%kmeansrep: repetitions in kmeans algorithm for the last step of IF-PCA,
%defacult to be 30. 
%per: a number with 0 < per <= 1, the percentage of Kolmogorov-Smirnov 
%statistics that will be used in the normalization step, default to be 1. 
%When the data is highly skewed, this parameter can be specified, such as
%0.5.
%
%Output: 
%label: n by 1 vector, as the estimated labels for each sample
%stats: 4 by 1 struct, including the important statistics in the algorithm
%as following:
%  stats.KS: p by 1 vector shows normalized KS value for each feature;
%  stats.HC: p by 1 vector shows the HC value for each feature;
%  stats.pval: p by 1 vector shows the p-value for each feature;
%  stats.ranking: p by 1 vector shows the ranking for each feature
%    according to ranking with p-values
%numselect: number of selected features in IF-PCA
%
%Example:
% load('lungCancer.mat');
% Data = [lungCancer_test(1:149, 1:12533); lungCancertrain(:, 1:12533)];
% Data = Data';
% [p, n] = size(Data);
% [label, stats, L] = ifpca(Data, 2, [], [], 100*p, 30);
%
%Reference: 
%Jin and Wang (2014): Important Features PCA for High Dimensional
%Clustering. 

[p, n] = size(Data);

% Error checking
if (nargin<2 || isempty(K))
  error 'Please include the number of clusters'
end

if (nargin<3||isempty(KSvalue))
    nullsimu = true; 
else 
    nullsimu = false;
end
    
if (nargin<4||isempty(pvalcut))
    pvalcut = (log(p))/p;
end

if ((nargin<5||isempty(rep)) && nullsimu)
    rep = 100*p;
end

if (nargin<6||isempty(kmeansrep))
    kmeansrep = 30;
end

if (nargin<7||isempty(per))
    per = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  Main Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Normalize Data
gm = mean(Data'); gsd = std(Data');
Data = (Data - repmat(gm', 1, n))./repmat(gsd', 1, n);


%Simulate KS values
if(nullsimu)
KSvalue = zeros(rep,1); kk = (0:n)'/n;
for i = 1:rep
    
	x = randn(n,1); 
	z = (x - mean(x))/std(x);
	z = z/sqrt(1 - 1/n);
	pi = normcdf(z);
	pi = sort(pi);

	KSvalue(i) = max(max(abs(kk(1:n) - pi)), max(abs(kk(2:(n+1)) - pi)));
	end
KSvalue = KSvalue*sqrt(n);
clear x z pi kk;
end
KSmean = mean(KSvalue); KSstd = std(KSvalue);
if (per < 1)
    KSvalue = sort(KSvalue, 'ascend');
    KSmean = mean(KSvalue(1:round(rep*per)));
    KSstd = std(KSvalue(1:round(rep*per)));
end

%Calculate KS value for each feature in the data set
kk = (0:n)'/n;
KS = zeros(p, 1);
for j= 1:p
    pi = normcdf(Data(j,:)/sqrt(1 - 1/n));
    pi = sort(pi);
    KS(j) = sqrt(n)*max(max(abs(kk(1:n) - pi')), max(abs(kk(2:(n+1)) - pi')));
    clear pi;
end

% Standardize KS value according to Efron's idea
if (per == 1)
    KS = (KS - mean(KS))/std(KS)*KSstd + KSmean;
else 
    KS = sort(KS, 'ascend');
    KSm = mean(KS(1:round(per*p))); KSs = std(KS(1:round(per*p)));
    KS = (KS - KSm)/KSs*KSstd + KSmean;
end
% Calculate p-value with simulated KS values
pval = zeros(p,1);
for i = 1:p
    pval(i) = mean(KSvalue > KS(i));
end
[psort, ranking] = sort(pval, 'ascend');


% Calculate HC functional at each data point
kk = (1:p)'/(1 + p);
HCsort = sqrt(p)*(kk - psort)./sqrt(kk);
HCsort  =   HCsort./sqrt(max(sqrt(n)*(kk - psort)./kk, 0) + 1 );
HC = zeros(p,1);
HC(ranking) = HCsort;

% Decide the threshold
Ind = find(psort>pvalcut, 1, 'first');
ratio = HCsort;
ratio(1:Ind-1) = -Inf; ratio(round(p/2)+1:end)=-Inf;
L = find(ratio == max(ratio), 1, 'last');
numselect = L; 

% Record the statistics for every feature
stats.KS = KS; stats.HC = HC; stats.pval = pval; stats.ranking = ranking; 

% IF-PCA
data_select = Data(pval <= psort(L), :);
signals = (1:p);
signals = signals(pval<=psort(L));
G = data_select'*data_select;
[V, ~] = eigs(G, K - 1); 
label = kmeans(V, K, 'replicates', kmeansrep);
ppp = psort(L);


end