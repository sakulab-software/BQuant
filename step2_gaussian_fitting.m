function step2_gaussian_fitting(header, dataset, visible)

close all;

file = sprintf('data_working/%s_step1_%02d.mat', header, dataset);
load(file,'data');

fprintf(1,'\nStep 2 : Data fitting by gaussians ...\n\n');

dx  = data.dx;             % scale  um/pixel
x = data.x;                % real scale coordinate
Y = data.smooth_I;         % smoothed fluorescent intencsity
g_idx = data.g_idx;        % estimated intensity index where MT exists
g_loc = x(g_idx);

N_gauss = length(g_idx);   % # of gauss funcs
N_data = length(x);        % # of data points

M = 100;                     % # of proposals in parameter estimation
max_reject_num = 2*N_gauss;  % for breaking loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The fitting error and movement of the Gaussian location are small, 
% and making it easy to adjust the Gaussian width and size.

%---
% likelihood hyper parameter (sd of error dist in likelihood)
%LL_sd = 5;
LL_sd = max(Y)*0.001;          % small error tolerance

%---
% Prior parameters

% mean of gaussian location
g_loc_mu  = g_loc';           % Estimation in step 1
% std of gaussian location
%g_loc_sd  = 0.01;
g_loc_sd  = 0.001;            % (um) Strong constraints

% mean of "gaussian width"
%g_sd_mu   = 2*dx;
g_sd_mu   = 1;                % 
% std of "gaussian sd"
%g_sd_sd   = 3*dx;
g_sd_sd   = 1;                %

% mean of "gaussian height"
g_hei_mu   = Y(g_idx)';       % Values at gaussian location
% sd of "gaussian height"
%g_hei_sd = 500;
g_hei_sd = max(Y)*0.1;        %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------
% Functions

LL = @(hei, loc, sd, base) ...
  sum( ...
        (repmat(Y',M,1) - model_calc(x, hei, loc, sd, base) ).^2 ...
        ./ LL_sd^2/2, 2 ...
     );         

P_loc = @(loc)  sum( (loc - g_loc_mu).^2, 2) ./ 2 ./ g_loc_sd^2;
P_sd = @(sd)  sum( (sd - g_sd_mu).^2, 2) ./ 2 ./ g_sd_sd^2;
P_hei = @(hei)  sum( (hei - g_hei_mu).^2, 2) ./ 2 ./ g_hei_sd^2;

Pos = @(hei, loc, sd, base)  LL(hei, loc, sd, base) + P_loc(loc) + P_sd(sd) + P_hei(hei);


%-------------------------------------------------
% initialize

hei = Y(g_idx)';
loc = g_loc';
sd  = g_sd_mu * ones(size(loc)) * 0.5;
base = 0;

E = 1e+100000;  % initial error

sweep_count = 1;
re = 1;
idx = 0;         % initial g-index to be fitted
direction = 1;   % local fitting moving direction 1: right, -1: left

total_reject = 0;

rng( sum(100*clock)+ feature('getpid'), 'twister');
scr = get(0, 'ScreenSize');
set(0, 'DefaultLineLineWidth', 1);
set(0, 'defaultAxesFontSize', 15);
set(0, 'defaultTextFontSize', 15);
fig = figure('Position', [100 20 scr(3)*0.8 scr(4)*0.8], ...
                'visible', visible);
             
%-------------------------------------------------
% fitting loop

while(1)
  
  % preparing prop_parameters
  test_hei = repmat(hei,M,1);
  test_loc = repmat(loc,M,1);
  test_sd  = repmat(sd,M,1);
  test_base   = base;

  % select successive part of the data points
  % for iteratively fitting data to the model.
  if(direction == 1)
    if(idx(end) + 1 <= N_gauss)
      idx = idx + 1;
    else
      direction = -1;
      idx = N_gauss-re+1:N_gauss;
    end
  else   % direction == -1
    if(idx(1) - 1 >= 1)
      idx = idx - 1;
    else
      direction = 1;
      sweep_count = sweep_count + 1;
      re = rem(sweep_count+1,2)+1;
      idx = 1:re;
    end
  end
  
  fprintf(1,sprintf('Intensity peak No. [%d:%d] ========\n',idx(1), idx(end)));
  plot_fig(header,dataset,x,Y,idx,g_idx,hei,loc,sd,base);
  
  k = length(idx);
  reject = 0;
  accept = 0;
  
  while(1) 
    % perturb parameters
    while(1)
      test_hei(:,idx) = repmat(hei(idx),M,1) .* normrnd(1,0.05, M, k);
      test_loc(:,idx) = repmat(loc(idx),M,1) .* normrnd(1,0.01, M, k);
      test_sd(:,idx) = repmat(sd(idx),M,1) .* normrnd(1,0.05, M, k);
      test_base = base + normrnd(0,0.0);
      if( all( abs(test_loc - repmat(g_loc',M,1)) < 10*dx)) break; end
    end
    
    test_E = Pos(test_hei, test_loc, test_sd, test_base);
    [min_E, min_E_idx] = min(test_E);
    
    if(min_E < E)
      hei = test_hei(min_E_idx,:);
      loc = test_loc(min_E_idx,:);
      sd = test_sd(min_E_idx,:);
      base = test_base;
      
      E = min_E;
      
      [gt, gx, gh] = plot_fig(header,dataset,x,Y,idx,g_idx,hei,loc,sd,base);
      
      fprintf(1,'%2.2f\n',E);
      
      accept = accept + 1;
      total_reject = 0;
      
      if( accept == 2)
        break;   % exit perturb parameters
      end
      
    else
      
      reject = reject + 1;

      if( reject == 5)
        total_reject = total_reject + 1;
        break;   % exit perturb parameters
      end
      
    end
    
  end % perturb parameters
  
  if( total_reject == max_reject_num)
    break;   % exit sliding window
  end
  
  
end  % sliding window


%-------------------------------------------------
% saving data

img_file = sprintf('data_working/%s_step2_%02d.png', header, dataset);
print(fig,'-dpng','-r300', img_file);

data_file = sprintf('data_working/%s_step2_%02d.mat', header, dataset);

data.E = E;        % fitting error (log posterior)
data.hei = hei;    % gaussian height vector
data.loc = loc;    % gaussian location vector
data.sd = sd;      % gaussian width vector
data.base = base;  % offset

data.Y  = Y;       % (smoothed) intensity vector
data.x  = x;       % coordinate vector
data.gt = gt;      % gaussian fitting vector
data.gx = gx;      % "gaussian coordinates" x "gaussian index" matrix
data.gh = gh;      % "each gaussian" x "gaussian index" matrix

save(data_file, 'data');


fprintf(1,'done\n');
fprintf(1,'Saved data and image in "data_working" folder.\n\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Z = model_calc(x, hei, loc, sd, base) 

N_gauss = size(hei,2);
N_data  = size(x,1);
M       = size(hei,1);   % # of proposal

X = repmat(x,1,N_gauss);
Z = zeros(M,N_data);

for i = 1:M
  
  HEI = repmat(hei(i,:),N_data,1);
  LOC = repmat(loc(i,:),N_data,1);
  SD  = repmat(sd(i,:),N_data,1);
  
  z = sum( HEI .* exp(-(X - LOC).^2 ./ SD.^2 / 2), 2);
  Z(i,:) = z' + base;
  
end

end


function  [gt, gx, gh] = plot_fig(header,dataset,x,Y,idx,g_idx,hei,loc,sd,base)

N_gauss  = length(g_idx);   % # of gauss funcs
N_data   = length(x);       % # of data points

% Sum of Gaussians (model of intensity data, 1D)
f =  @(hei, loc, sd, base) ...
  sum( repmat(hei, N_data, 1) .* ...
       exp( ...
           - (repmat(x,1,N_gauss)  - repmat(loc,N_data,1)).^2 ...
            ./ repmat(sd, N_data, 1).^2/2 ...
          ) , 2 ...
     ) + base;

%------

bh = max(Y)*0.1;
   
clf;
subplot(2,1,1);

gt = f(hei, loc, sd, base);

p1 = plot(x,Y,'k-');
hold on;
p2 = plot(x,gt,'r-');
p3 = plot( [loc(1); loc(1)], [0; bh], 'k--');
plot( [loc; loc], [zeros(1,length(loc)); bh*ones(1,length(loc))], 'k--');
p4 = plot( [loc(idx(1)); loc(idx(1))], [0; bh], 'b-', 'LineWidth', 3);
plot( [loc(idx); loc(idx)], [zeros(1,length(idx)); bh*ones(1,length(idx))], 'b-', 'LineWidth', 3);

xlabel('Location (um)'); ylabel('Intensity');
legend([p1 p2 p3 p4], 'Measurement', 'Model', 'Gaussian location', 'Adjustment target', 'Location', 'best');
xlim([0 x(end)]);
title(sprintf('dataset %d in "%s"',dataset,header));

subplot(2,1,2);

x_center  = repmat(loc,61,1);
x_diff = repmat([-30:30]',1,N_gauss) .* repmat(5*sd/30,61,1);
gx = x_center + x_diff;
gh  = repmat(hei,61,1).* exp(-x_diff.^2./repmat(sd,61,1)/2);

h1 = plot(gx, gh, 'k-');
hold on;
h2 = plot(x, gt,'r-');
legend([h1(1), h2],{'Each Gaussian', 'Sum Gaussain'}, 'Location', 'best');
xlabel('Location (um)'); ylabel('Intensity');
xlim([0 x(end)]);


drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
