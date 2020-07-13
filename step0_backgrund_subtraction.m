function header = step0_backgrund_subtraction(data_file, dx, pp)

% dx        = 0.267; 
% data_file = 'sample1.csv';
% pp        = 10;             % paraboloid parameter

da = csvread(sprintf('data/%s',data_file));  
header = erase(data_file, '.csv');
da(da==0) = nan;

% make working data directry
if(~exist('data_working', 'dir'))
  mkdir('data_working');
end

scr = get(0, 'ScreenSize'); 
set(0, 'DefaultLineLineWidth', 1);
set(0, 'defaultAxesFontSize', 15);
set(0, 'defaultTextFontSize', 15);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------

fprintf(1,'\nStep 0 : Background subtraction ...');


for series_idx = 1:size(da,2)
  
  y = rmmissing(da(:,series_idx));
  x = [1:length(y)]'*dx;
  
  yup = (max(y) - min(y))/1000;
  y0 = zeros(length(x),1);
  
  % sliding paraboloid
  for i = 1:length(x)

    y0_tmp = min(y)*0.5;
    while(1)
      f = -(x-x(i)).^2/pp + y0_tmp; % paraboloid
      d = y - f;
      
      if(min(d) > 0)
        y0_tmp = y0_tmp + yup;
      else
        y0(i) = y0_tmp - yup;
        break;
      end
    end
  end
  
  bg_bias = zeros(length(x),1);
  for i = 1:length(x)
    ff = -(x-x(i)).^2/pp + y0;
    bg_bias(i) = max(ff);
  end
  
  % Draw
  fig = figure('Position', [100 20 scr(3)*0.3 scr(4)*0.5], ...
                'visible', 'off');
  clf;
  subplot(2,1,1);
  plot(x, y);
  hold on;
  plot(x, bg_bias);
  legend('Raw', 'Background');
  xlabel('Location');  ylabel('Intencity');
  title('Measured data');
  
  subplot(2,1,2);
  plot(x, y-bg_bias);
  xlabel('Location');  ylabel('Intencity');
  title('Subtracted data');
  
  % Save
  file = sprintf('data_working/%s_step0_%02d.mat', header, series_idx);
  data.x = x;
  data.y = y;
  data.bg_bias = bg_bias;
  data.dx = dx;
  data.header = header;
  
  save(file, 'data');

  img_file = sprintf('data_working/%s_step0_%02d.png', header, series_idx);
  print(fig,'-dpng','-r300', img_file);
  
end


fprintf(1,'done\n');
fprintf(1,'Saved data and image in "data_working" folder.\n\n');





