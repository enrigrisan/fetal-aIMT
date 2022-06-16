function [center, U, obj_fcn] = fcm_my(data, cluster_n, options)
% FCM Find clusters with fuzzy c-means clustering.
% FCM Fuzzy c-means clustering.
%   Synopsis
%   [center,U,obj_fcn] = fcm(data,cluster_n)
%   Description
%   [center, U, obj_fcn] = fcm(data, cluster_n) applies the fuzzy c-means
%   clustering method to a given data set.
%   The input arguments of this function are:
%   data: data set to be clustered; each row is a sample data point
%   cluster_n: number of clusters (greater than one)
%   The output arguments of this function are:
%   center: final cluster centers, where each row is a center
%   U: final fuzzy partition matrix (or membership function matrix)
%   obj_fcn: values of the objective function during iterations
%   fcm(data,cluster_n,options) uses an additional argument variable, options,
%   to control clustering parameters, introduce a stopping criteria, and/or set
%   the iteration information display:
%   options(1): exponent for the partition matrix U (default: 2.0)
%   options(2): maximum number of iterations (default: 100)
%   options(3): minimum amount of improvement (default: 1e-5)
%   options(4): info display during iteration (default: 1)
%   If any entry of options is NaN, the default value for that option is used
%   instead. The clustering process stops when the maximum number of iteration
%   is reached, or when the objective function improvement between two
%   consecutive iteration is less than the minimum amount of improvement
%   specified.
%   Example
%   data = rand(100, 2);
%   [center,U,obj_fcn] = fcm(data, 2);
%   plot(data(:,1), data(:,2),'o');
%   maxU = max(U);
%   index1 = find(U(1,:) == maxU);
%   index2 = find(U(2, :) == maxU);
%   line(data(index1,1),data(index1,2), ...
%   	'marker','*','color','g');
%   line(data(index2,1),data(index2,2), ...
%   	'marker','*','color','r');
%
%   See also FCMDEMO, INITFCM, IRISFCM, DISTFCM, and STEPFCM.

%   Roger Jang, 12-13-94.
%   Copyright 1994-2000 The MathWorks, Inc. 
%   $Revision: 1.9 $  $Date: 2000/06/15 13:33:43 $

if nargin ~= 2 & nargin ~= 3,
	error('Too many or too few input arguments!');
end

data_n = size(data, 1);
in_n = size(data, 2);

% Change the following to set default options
default_options = [2;	% exponent for the partition matrix U
		100;	% max. number of iteration
		1e-5;	% min. amount of improvement
		0];	% info display during iteration 

if nargin == 2,
	options = default_options;
else
	% If "options" is not fully specified, pad it with default values.
	if length(options) < 4,
		tmp = default_options;
		tmp(1:length(options)) = options;
		options = tmp;
	end
	% If some entries of "options" are nan's, replace them with defaults.
	nan_index = find(isnan(options)==1);
	options(nan_index) = default_options(nan_index);
	if options(1) <= 1,
		error('The exponent should be greater than 1!');
	end
end

expo = options(1);		% Exponent for U
max_iter = options(2);		% Max. iteration
min_impro = options(3);		% Min. improvement
display = options(4);		% Display info or not

obj_fcn = zeros(max_iter, 1);	% Array for objective function

U = initfcm_my(cluster_n, data,expo);			% Initial fuzzy partition
% Main loop
for i = 1:max_iter,
	[U, center, obj_fcn(i)] = stepfcm(data, U, cluster_n, expo);
	if display, 
		fprintf('Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
	end
	% check termination condition
	if i > 1,
		if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro, break; end,
	end
end

iter_n = i;	% Actual number of iterations 
obj_fcn(iter_n+1:max_iter) = [];
