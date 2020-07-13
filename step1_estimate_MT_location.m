function N = step1_estimate_MT_location(header,smooth_level)

%smooth_level = 1;      % SD of Gaussian filter

scr = get(0, 'ScreenSize'); 
set(0, 'DefaultLineLineWidth', 1);
set(0, 'defaultAxesFontSize', 15);
set(0, 'defaultTextFontSize', 15);

% load data
da = dir(sprintf('data_working/%s_step0_*.mat',header));

fprintf(1,'\nStep 1 : Data smoothing and Peak estimation ...');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(da,1);      % number of data set

for series_idx = 1:N     % index of fluorescence measurement
  
  file = da(series_idx).name;
  load(sprintf('data_working/%s', file), 'data');
  % x, y, bg_bias dx, header;
  
  org_I = data.y - data.bg_bias;
  data.org_I = org_I;
  dx = data.dx;
  
  %---------------------------------------------
  % smoothing line and finding peaks
  
  % Gaussian smoothing to remove small measurement noise
  kern = normpdf(-3*smooth_level:3*smooth_level, 0, smooth_level);
  kern = kern / sum(kern);
  smooth_I = conv(org_I,kern, 'same');
  org_I = org_I - min(org_I);            % remove background baseline
  %smooth_I = smooth_I - min(smooth_I);
  data.smooth_I = smooth_I;
  
  % find peaks
  [pos_idx, neg_idx] = find_peaks(smooth_I);
  
  % If there are negative peaks at both ends, truncate the ends for fitting
  init_idx = 1; end_idx = length(org_I);
  if(neg_idx(1) < pos_idx(1))
    init_idx = neg_idx(1);
  end
  if(neg_idx(end) > pos_idx(end))
    end_idx = neg_idx(end);
  end
  org_I = org_I(init_idx:end_idx);
  smooth_I = smooth_I(init_idx:end_idx);
  x = dx * (1:size(org_I,1));
  
  [pos_idx, ~] = find_peaks(smooth_I);
  
  %---------------------------------------------
  % find additional fluorescence peaks
  
  dY = diff(smooth_I);
  
  % estimate intensity index where MT exists as inflection point
  [dpos_idx, dneg_idx] = find_peaks(dY);
  
  z = [dpos_idx ones(size(dpos_idx)); dneg_idx -ones(size(dneg_idx))];
  [~, idx] = sort(z(:,1));
  z = z(idx,:);
  z = [z,  dY(z(:,1)),   sign(z(:,2) .* dY(z(:,1)))];
  
  % remove small inflection point
  e = 3;
  for i = 1:length(z(:,1))
    
    if( z(i,4) == -1)
      if(i == 1)
        if( abs(dY(z(i,1)) - dY(z(i+1,1))) < e)
          z(i,4) = 1;
        end
      elseif(i == length(z(:,1)))
        if(abs(dY(z(i,1)) - dY(z(i-1,1))) < e)
          z(i,4) = 1;
        end
      else
        if(abs(dY(z(i,1)) - dY(z(i-1,1))) < e || ...
            abs(dY(z(i,1)) - dY(z(i+1,1))) < e)
          z(i,4) = 1;
        end
      end
    end
  end
  
  % inflection index
  infl_idx = z(z(:,4) == -1, 1);
  
  g_idx = [pos_idx; infl_idx];
  g_idx = sort(g_idx);
  data.g_idx = g_idx;
  
  %---------------------------------------------
  % Draw
  
  fig = figure('Position', [100 20 scr(3)*0.9 scr(4)*0.8], ...
                'visible', 'off');
  clf;
  subplot(2,1,1);
  h(1) = plot(x, org_I, 'k-', 'LineWidth', 0.5); hold on;
  h(2) = plot(x, smooth_I,'r-');
  for i = 1:length(pos_idx)
    plot([pos_idx(i) pos_idx(i)]*dx, [0 smooth_I(pos_idx(i))], ...
      'r--', 'LineWidth', 0.5);
  end
  h(3) = plot(pos_idx*dx, smooth_I(pos_idx), 'ro', ...
    'MarkerSize', 5, 'MarkerFaceColor', 'r');
  for i = 1:length(infl_idx)
    plot([infl_idx(i) infl_idx(i)]*dx, [0 smooth_I(infl_idx(i))], ...
      'b--', 'LineWidth', 1);
  end
  h(4) = plot(infl_idx*dx, smooth_I(infl_idx), 'bo', ...
    'MarkerSize', 5, 'MarkerFaceColor', 'b');
  legend(h(1:4), {'Measured', 'Smoothed', 'Peak', 'Additional peak'});
  xlabel('Location'); ylabel('Intensity');
  title('Comparison between measured and smoothed intensity');
  
  subplot(2,1,2);
    h(1) = plot(x(1:end-1), dY,'b-'); hold on;
  h(2) = plot([0 x(end)], [0 0], 'k--');
  for i = 1:length(infl_idx)
    plot([infl_idx(i) infl_idx(i)]*dx, [0 dY(infl_idx(i))], ...
      'k-', 'LineWidth', 2);
  end
  h(3) = plot(infl_idx*dx, dY(infl_idx),'ko', 'MarkerSize', 5);
  legend(h([1 3]), {'Derivative', 'Detected inflection point'})
  xlabel('Location'); ylabel('Intensity (1/um)');
  title('Derivative of smoothed intensity');
      
  
  %---------------------------------------------
  % saving data
  
  img_file = sprintf('data_working/%s_step1_%02d.png',header, series_idx);
  print(fig,'-dpng','-r300', img_file);
    
  data_file = sprintf('data_working/%s_step1_%02d.mat', header, series_idx');
  save(data_file, 'data');

  % dx       : scale  um/pixel
  % x        : real scale coordinate
  % org_I    : original fluorescent intencsity
  % smooth_I : smoothed fluorescent intencsity
  % g_idx    : estimated intensity index where MT exists
  
end


fprintf(1,'done\n');
fprintf(1,'Saved data and image in "data_working" folder.\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos_idx, neg_idx] = find_peaks(y)

dy = sign(diff(y));
tmp = [0; dy(1:end-1)];
flag = dy + tmp;
  
peak_idx = find(flag ==0);
idx = 1:length(peak_idx);

if(dy(1) > 0)
  pos_idx = peak_idx(find(mod(idx,2)));
  neg_idx = peak_idx(find(~mod(idx,2)));
else
  neg_idx = peak_idx(find(mod(idx,2)));
  pos_idx = peak_idx(find(~mod(idx,2)));
end
  