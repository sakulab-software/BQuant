clear variables;

% parameters

%--------
% step 0
file = 'sample1.csv';    % Specify the data file that should be in the 'data' folder
dx   =  0.267;           % Observation interval
pp   =  10;              % Paraboloid parameter for background subtraction

%--------
% step 1
smooth_level = 1;        % Larger and smoother


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

%--------
% step 0
h = step0_backgrund_subtraction(file,dx,pp);

%--------
% step 1
N = step1_estimate_MT_location(h, smooth_level);

%--------
% step 2
for dataset = 1:N
  step2_gaussian_fitting(h, dataset, 'on');
  % Turn the last 'off' if you do not see the fitting
end

%--------
% step 3
step3_quantification(h, N);