function U = initfcm_my(cluster_n, data,expo)
%INITFCM Generate initial fuzzy partition matrix for fuzzy c-means clustering.
%   U = INITFCM(CLUSTER_N, DATA_N) randomly generates a fuzzy partition
%   matrix U that is CLUSTER_N by DATA_N, where CLUSTER_N is number of
%   clusters and DATA_N is number of data points. The summation of each
%   column of the generated U is equal to unity, as required by fuzzy
%   c-means clustering.
%
%       See also DISTFCM, FCMDEMO, IRISFCM, STEPFCM, FCM.

%   Roger Jang, 12-1-94.
%   Copyright 1994-2000 The MathWorks, Inc. 
%   $Revision: 1.9 $  $Date: 2000/06/15 13:33:45 $

datamin=min(data);
datamax=max(data);
for ct=1:cluster_n
    center(ct,:)=ct*(datamax-datamin)/(cluster_n+1);
end;

dist = distfcm(center, data);       % fill the distance matrix
dist(find(dist==0))=1;
tmp = dist.^(-2/(expo-1));      % calculate new U, suppose expo != 1

U = tmp./(ones(cluster_n, 1)*sum(tmp));
