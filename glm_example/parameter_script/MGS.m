% VSWMDA Randomization

% Four types of images
% Four quadrants
% Four order positions
% The total number of trials has to be multiples of 16!!

clear
close
clc

for ses = 1%:5
    
runs = [1:2] + (ses-1)*2;

for rr = runs % two runs
    ttl = ['MGS_run' num2str(rr)]; % figure title
n = 16;

%% Trial length - TR = 1000 ms
tr = 1;

srt_t = 0;%0.5;
fix_t = 1.5;%tr*2-srt_t;

stim_t = 0.5;
%blank_t = 0.15;
delay_t = [8 10];
%cue_t = 0.5;
resp_t = tr;
%feed_t = 0.5;
ITI = tr*[8:10]'; % 9 s on average


trial_len = stim_t + mean(delay_t) + resp_t + mean(ITI) % feed_t
%trial_len*48/60 % 20 mins

run_len = (trial_len*16 + 25)/60

%% Now let's assign different locations

% % Behavior
% d = 67.7; % cm
% h = 30.5;
% w = 38;
% x_res = 1280;
% y_res = 1024;

% fMRI
d = 138.9; % cm
h = 35;
w = 43;
x_res = 1024;
y_res = 768;

%% Size
jit_ecc = 0.5;
jit_ang = 5;

fix_size = 0.5;
dot_size = 0.5;
ROI_size = 0.8; % image size

back_color = {'154,154,154,255'};
dot_color = {'0,0,0,255'};
fix_color = {'0,0,0,255'};
feed_color = {'50,205,50,255'};

ecc_all = [3 5.5]; ecc = ecc_all(mod(rr+1,2)+1);
ang = [15:20:75];%11.25:22.5:360;
ang = [ang ang+90 ang+180 ang+270];

%% Plot
figure(1);clf
scatter(6*cosd(ang),6*sind(ang))
line([-6 6],[0 0])
line([0 0],[-6 6])
%saveas(figure(1),'location16.jpg')

% %% Run optseq to get trial sequence
% pause % run optseq and import data
% quad = table2array(MGS1(1:2:end,2)); % quadrant number
% iti_all = table2array(MGS1(2:2:end,3)); % all the ITIs
% 
% save('MGS2500_001_optseq.mat','MGS1','iti_all','quad')

%% quad & iti
% each run, 16 trials, all 16 locations
quad = repmat([1:4]',n/4,1); quad = Shuffle(quad);
iti_all = datasample(ITI,n,1);

%% Assign different angles to each trial and make sure there are no 
% repetitive locations and images in sequential trials
dif = 1;ang_all = [];img_all = [];img_sub = []; img_sub_all = [];
while dif
for ii = 1:4
ang_quad{ii} = repmat(ang((1:4)+4*(ii-1)),1,n/16);%[repmat(ang((1:4)+ii-1),1,3),ang(randperm(4,3)+ii)];
if rr == 1|3|5|7|9|11|13|15
img_quad{ii} = [ones(1,n/4/4),2*ones(1,n/4/4),3*ones(1,n/4/4),4*ones(1,n/4/4)];
elseif rr == 2|4|6|8|10|12|14|16
    img_quad{ii} = [3*ones(1,n/4/4),4*ones(1,n/4/4),1*ones(1,n/4/4),2*ones(1,n/4/4)];
end
img_sub{ii} = (mod(rr+1,2)+1)*ones(n/4,1);%repmat([1;2],n/4/2,1);
ind = find(quad == ii);
ind_ang = randperm(n/4,n/4);
ang_all(ind) = ang_quad{ii}(ind_ang);
img_all(ind) = img_quad{ii}(ind_ang);
img_sub_all(ind) = img_sub{ii}(ind_ang);
end
dif = any(~diff(ang_all)) & any(~diff(img_all));
end
% %% Assign different images to each trial
% img_sub = datasample(1:2,n);
img_sub = img_sub_all;

%% delay
% delay duration
delay = zeros(n,1);
delay_t_io = repmat(delay_t',2,1);
for jj = 1:4 % each quadrant
       temp = find(ismember(ang_all, ang_quad{jj})); % always 4 trials % ret_ord
       delay(temp,1) = Shuffle(delay_t_io);
end

%% Construct the full experiment
run = rr*ones(n,1);
trial = [1:n]';

% actual image name
img_name = [];
for ii = 1:n % trial
img_name{ii,1} = [num2str(img_all(ii)) '_' num2str(img_sub(ii)) '.png'];
end

% fixation location
fix_x = x_res/2*ones(n,1);
fix_y = y_res/2*ones(n,1);

% time
iti_t_all = Shuffle([ITI(1)*ones(5,1);ITI(2)*ones(5,1);ITI(3)*ones(6,1)]); % iti_t_all(1) = 5.5; % datasample(ITI,n)
fix_t_all = fix_t*ones(n,1) + iti_t_all;
stim_t_all = fix_t_all + stim_t*ones(n,1);
delay_t_all = stim_t_all + datasample(delay_t,n)';%delay_t*ones(n,1);
resp_t_all = delay_t_all + resp_t*ones(n,1);

% size
fix_size_pix = tand(fix_size).*d./w.*x_res;
fix_size_pix_all = fix_size_pix*ones(n,1);
dot_size_pix = tand(dot_size).*d./w.*x_res;
dot_size_pix_all = dot_size_pix*ones(n,1);
roi_size_pix = tand(ROI_size).*d./w.*x_res;
roi_size_pix_all = roi_size_pix*ones(n,1);

% color
dot_color_all = repmat(dot_color,n,1);
back_color_all = repmat(back_color,n,1);
fix_color_all = repmat(fix_color,n,1);
feed_color_all = repmat(feed_color,n,1);

% targets actual location
jit = linspace(-1,1,16); 
jit = jit(randperm(16)); jit = reshape(jit,n,1);
ecc_rand = ecc + jit_ecc*jit;%(2*rand(n,4)-1);

jit = linspace(-1,1,16); 
jit = jit(randperm(16)); jit = reshape(jit,n,1);
ang_rand = ang_all' + jit_ang*jit;%(2*rand(n,4)-1);

% ecc_rand = ecc + jit_ecc*(2*rand(n,1)-1); 
% ang_rand = ang_all' + jit_ang*(2*rand(n,1)-1);
x = tand(ecc_rand).*d.*cosd(ang_rand)./w.*x_res + x_res/2;
y = y_res/2-tand(ecc_rand).*d.*sind(ang_rand)./h.*y_res;

% full table
mltable = table(run,trial,quad,img_name,...
    ecc_rand,ang_rand,x,y,fix_x,fix_y,...
   fix_t_all,stim_t_all,delay_t_all,resp_t_all,iti_t_all,...
    fix_size_pix_all,dot_size_pix_all,roi_size_pix_all,ang_all',ecc*ones(n,1),...
    'VariableNames',{'run_ind','trial_ind','stim_quad','stim_img','stim_ecc','stim_ang','stim_x','stim_y','fix_x','fix_y', ...
    'fix_t','stim_t','delay_t','resp_t','iti_t','fix_size','dot_size','roi_size','stim_ang_orig','ecc'});

% adjust baseline timing
mltable(17,:) = mltable(16,:);
mltable.iti_t(17) = 15;
mltable{1,11:15} = mltable{1,11:15}-mltable.iti_t(1)+5;
mltable.iti_t(1) = 5;
mltable.trial_ind = [1:17]';

% mltable(2:17,:) = mltable(1:16,:);
% mltable(18,:) = mltable(17,:);
% mltable.trial_ind = [1:18]';

%% save the table
% Write into a csv file
writetable(mltable,[ttl '.csv'])

%% Calculate measurements
% % tm = mltable(2:17,[11:16]);
% tm = mltable(:,[11:15]);
% tm_all = sum(table2array(tm(:,4)));
% tm_all_min = tm_all/60;
% n_vol = tm_all/tr;
% n_vol_runs(rr) = n_vol;

 % calculate measurements
 tm = mltable(:,[11:15]);
tm_all = sum(table2array(tm(1:16,4))) + table2array(tm(17,end));
tm_all_min = tm_all/60;
n_vol = tm_all/1
n_vol_runs(rr) = n_vol;

%% plot actual stimuli locations
temp = [reshape(x,[],1) reshape(y,[],1)];
figure(3);clf
%for ii = 1:12
scatter(temp(:,1),y_res-temp(:,2),5)
%hold on
%end
line([0 x_res],[y_res/2 y_res/2],'Color','k')
line([x_res/2 x_res/2],[0 y_res],'Color','k')
xlim([0 x_res])
ylim([0 y_res])
title('Actual Stimuli Locations')
saveas(figure(3),[ttl,'.jpg'])
end
end
%writematrix(n_vol_runs,'MGS_TR1000_n_vol.csv')
% %% 24 sets of locations
% % Let's find trials that probe the first quadrant
% ind = find(ret_quad == 1); % 32 trials
% % We want 8 trials per location in that quadrant (so each location probed
% % equally), AND we want each location appears at the same chance
% % Let's construct possible sets of locations (assymetrical & at least 45
% % degrees apart)
% ang_proto_1 = [];
% for ii = 1:4 % the first quadrant, four angles
%     ang_1 = ang(ii);
%     % Let's calculate possible second quadrant locations
%     temp = ang(5:8);
%     ang_2_all = temp((temp + ang_1)/2 ~= 90);
%     for jj = 1:length(ang_2_all) % for each second-quadrant locations
%         ang_2 = ang_2_all(jj);
%         % Possible third quadrant locations
%         temp = ang(9:12);
%         ang_3_all = temp((temp + ang_2)/2 ~= 180 & (temp - ang_1) ~= 180);
%         for kk = 1:length(ang_3_all) % for each third-quadrant locations
%             ang_3 = ang_3_all(kk);
%             temp = ang(13:16);
%             ang_4 = temp((temp+ang_3)/2~= 270 & temp-(ang_2-90)~=270 & temp+ang_1~=360);
% 
%             % store all the possible sets of locations
%             ang_proto_1 = [ang_proto_1;ang_1 ang_2 ang_3 ang_4];
%         end
%     end
% end
% 
% % each location appears equally?
% temp = histcounts(reshape(ang_proto_1,1,[]),[0 ang+2]);
% 
% % number of trials per order per quadrant
% ord_temp = ret_ord(ind); % 8 trials
% % for each order, randomly select 2 sets of location with each 4 locations
% % in that quadrant
% ang_proto_1(1:6,:)= ang_proto_1(randperm(6,6),:);
% ang_proto_1([1:6]+6,:) = ang_proto_1(randperm(6,6)+6,:);
% ang_proto_1([1:6]+12,:) = ang_proto_1(randperm(6,6)+12,:);
% ang_proto_1([1:6]+18,:) = ang_proto_1(randperm(6,6)+18,:);
% 
% %%
% ang_proto_1 = [];
% for ii = 1:4 % the first quadrant, four angles
%     ang_1 = ang(ii);
%     % Let's calculate possible second quadrant locations
%     temp = ang(5:8);
%     ang_2_all = temp((temp + ang_1)/2 ~= 90);
%     for jj = 1:length(ang_2_all) % for each second-quadrant locations
%         ang_2 = ang_2_all(jj);
%         % Possible third quadrant locations
%         temp = ang(9:12);
%         ang_3_all = temp((temp + ang_2)/2 ~= 180 | (temp - ang_1) ~= 180);
%         for kk = 1:length(ang_3_all) % for each third-quadrant locations
%             ang_3 = ang_3_all(kk);
%             temp = ang(13:16);
%             ang_4 = temp((temp+ang_3)/2~= 270 | temp-(ang_2-90)~=270 | temp+ang_1~=360);
% for ll = 1:length(ang_4)
%             % store all the possible sets of locations
%             ang_proto_1 = [ang_proto_1;ang_1 ang_2 ang_3 ang_4(ll)];
% end
%         end
%     end
% end
% % randomly pick 32 out of 60 sets of locations, with equal probability in
% % probing each location
% 
% %         % Now we start to combine order, quadrant, and images
% %         for oo = 1:4 % for each specific image-order association
% %             % We have four possible choices of image-quadrant associations
% %             for qq = 1:4
% %                 % Let's store those 4*4 = 16 combinations
% %                 order_spec(4*(oo-1)+qq,:) = order_temp(oo,:); % order
% %                 quad_spec(4*(oo-1)+qq,:) = quad_temp(qq,:); % quadrant
% %             end
% %         end
% %         % But only 4 of those 4*4 = 16 combinations are valid, the
% %         % combination of those 4 make sure that each image is
% %         % probed equally in each order & each quadrant
% %         dif = order_spec - quad_spec; dif = sort(dif,2);
% %         uniq_dif = unique(dif,'rows');
% %         for kk = 1:length(uniq_dif)
% %             occurence(kk) = sum(ismember(uniq_dif(kk,:),dif));
% %         end
%         
% %         % check: each quadrant appears in each order equally
% %         a = arrayfun(@(y)(sum(order_rs(quad_rs == 4) == y)),1:4);
% %         quad_rs = reshape(quad_spec(:,1),[],1);
% %         order_rs = reshape(order_spec(:,1),[],1);
% %         a = arrayfun(@(y)(sum(order_rs == y)),1:4);
