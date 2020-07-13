function step3_quantification(header, N_data)

fprintf(1,'\nStep 3 : Statistics of height,width, ans area ...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir = 'data_stat';
if(~exist(dir)) mkdir(dir); end

H = nan(N_data,200);   % height
S = nan(N_data,200);   % width
A = nan(N_data,200);   % area

%--------------------------------------------------------------
% Data storage


for dataset = 1:N_data
  
  file = sprintf('data_working/%s_step2_%02d.mat', header, dataset);
  load(file,'data');
  
  hei = data.hei;    % gaussian height vector
  loc = data.loc;    % gaussian location vector
  sd  = data.sd;     % gaussian width vector
  
  H(dataset,1:length(hei)) = hei;
  S(dataset,1:length(sd)) = sd;
  A(dataset,1:length(sd)) = sqrt(2*pi) * sd .* hei;
end

%--------------------------------------------------------------
% Draw histogram 

scr = get(0, 'ScreenSize'); 
fig = figure('Position', [100 100 scr(3)*0.7 scr(4)*0.3]); clf; 

subplot(1,3,1); 

bin = max(max(H(:)))/10;
x = H(:); x = x(~isnan(x)); p = gamfit(x);
h = histogram(x, 'Normalization','pdf','BinWidth', bin);
hold on;
t = 0:max(x)/100:max(x);
pdf = gampdf(t, p(1), p(2));
plot(t, pdf, 'r-');
xlabel('Value'); ylabel('p.d.f.');
title('Height');
legend('Measurement', 'Gamma', 'Location', 'best');

subplot(1,3,2); 

bin = max(max(S(:)))/10;
x = S(:); x = x(~isnan(x)); p = gamfit(x);
h = histogram(x, 'Normalization','pdf','BinWidth', bin);
hold on;
t = 0:max(x)/100:max(x);
pdf = gampdf(t, p(1), p(2));
plot(t, pdf, 'r-');
xlabel('Value'); ylabel('p.d.f.');
title('Width');
legend('Measurement', 'Gamma', 'Location', 'best');

subplot(1,3,3); 

bin = max(max(A(:)))/10;
x = A(:); x = x(~isnan(x)); p = gamfit(x);
h = histogram(x, 'Normalization','pdf','BinWidth', bin);
hold on;
t = 0:max(x)/100:max(x);
pdf = gampdf(t, p(1), p(2));
plot(t, pdf, 'r-');
xlabel('Value'); ylabel('p.d.f.');
title('Area');
legend('Measurement', 'Gamma', 'Location', 'best');


% save data

img_file = sprintf('data_stat/%s_hist.png', header); 
print(fig,'-dpng','-r300', img_file);

outfile = sprintf('data_stat/%s_height_stat.csv', header); 
csvwrite(outfile, H);
outfile = sprintf('data_stat/%s_width_stat.csv', header); 
csvwrite(outfile, S);
outfile = sprintf('data_stat/%s_area_stat.csv', header); 
csvwrite(outfile, A);

fprintf(1,'done.\n\n');
