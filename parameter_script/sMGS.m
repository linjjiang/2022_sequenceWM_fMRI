% VSWMDA Randomization
% smgs, optseq version
% 10/28/2022
% finish a full set of randomization within run
% potential issues: each run is fixed at one eccentricity

% Four types of images
% Four quadrants
% Four order positions
% The total number of trials has to be multiples of 16!!

% This is just for one session (4 runs). To make multiple sessions, you
% need to change the run number below
clear
close
clc

for ses = 1%:4
    
runs = [1:4] + (ses-1)*4;

%% 1. Load locations
ses = ceil(runs(1)/4); % session number
% if session = 1, load set 1 & 2
ang_set1 = load(['../location/set',num2str(2*(ses-1)+1),'.mat']); 
ang_set2 = load(['../location/set',num2str(2*(ses-1)+2),'.mat']);
ang32 = [ang_set1.ang_set(randperm(16,16),:);ang_set1.ang_set(randperm(16,16)+16,:);...
    ang_set2.ang_set(randperm(16,16),:);ang_set2.ang_set(randperm(16,16)+16,:)];
% ang32 = [ang_set1.ang_set;ang_set1.ang_set;...
%     ang_set2.ang_set;ang_set2.ang_set];

% %% 2. Load optseq timing - we will adjust the delay and Iti timing when we shuffle the trials at last
% for ii = runs % we want to load run 1-16
% data = importfile(['../optseq/smgs/run',num2str(ii),'-001.par']);
% optseq(:,ii) = table2array(data(2:2:end,3))-15; % how many sec should be added to the end of delay or ITI
% end
% 
% % set up delay & ITI durations
% iti_add = zeros(size(optseq)); dly_add = iti_add;
% iti_add (optseq==0 | optseq==1) = optseq(optseq==0 | optseq==1);
% dly_add (optseq==0 | optseq==1) = 0;
% 
% ind = find(optseq > 1); 
% % we will use a psuedorandom number to indicate whether we need to add
% % 2-sec delay
% long_dly = ind(rand(1,length(ind))>0.5);
% dly_add(long_dly) = 2; 
% iti_add(ind) = optseq(ind)-dly_add(ind);

%% 2. Construct parameters in each trial
n = 16; % We have 16 trials per run
for rr = runs % No. of runs

ttl = ['sMGS_run' num2str(rr)]; % figure/file title

%% 2.1 Size of the stimuli
% % locations
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

jit_ecc = 0.5; % we want to have plus minus 0.5 dva jittering in eccentricity
jit_ang = 5; % we want to have plus minus 5 degrees of arc jittering in polar angle

fix_size = 0.5; % fixation cross size, in dva
dot_size = 0.5; % dot size
ROI_size = 0.8; % image size

back_color = {'154,154,154,255'}; % background color
dot_color = {'0,0,0,255'}; % dot color
fix_color = {'0,0,0,255'}; % fixation color
%feed_color = {'50,205,50,255'};

%% 2.2 Timing of the trial (for TR = 1000 ms)
tr = 1; % TR

srt_t = 0;% Start of the trial
fix_t = 1.5;%tr*2-srt_t; Fixation

stim_t = 0.5; % Stimulus duration
blank_t = 0.5; % Gap duration
%fix_t = tr*2-stim_t*4-blank_t*3

delay_t = [tr*8 tr*10]'; % possible delays, we will take care of this later

cue_t = 0.5; % cue
resp_t = 0.5; %tr; % response

%feed_t = 0.5;
ITI = tr*[9:12]'; % possible ITIs, we will take care of this later

trial_len = fix_t + stim_t*4 + blank_t*3 + mean(delay_t) + cue_t + resp_t + mean(ITI); % feed_t
%trial_len*48/60 % 20 mins

run_len = (trial_len*16  + 25)/60;

%% quadrant & image orders
% Let's construct three runs first
% order of the quadrants: we will always use a z-shape pattern
quad_subgrp = {[2 3 1 4;1 4 2 3;3 2 4 1;4 1 3 2];[2 3 1 4;1 4 2 3;3 2 4 1;4 1 3 2];...
    [2 3 1 4;1 4 2 3;3 2 4 1;4 1 3 2];[2 3 1 4;1 4 2 3;3 2 4 1;4 1 3 2]};
% ord_subgrp = {[4 1 3 2;3 2 4 1;1 4 2 3;2 3 1 4];[2 3 1 4;1 4 2 3;3 2 4 1;4 1 3 2];...
%     [3 2 4 1;4 1 3 2;2 3 1 4;1 4 2 3];[1 4 2 3;2 3 1 4;4 1 3 2;3 2 4 1]};

% order of the images: we will always use faces first or houses first
ord_subgrp = {[1 2 3 4;1 2 3 4;1 2 3 4;1 2 3 4];[2 1 4 3;2 1 4 3;2 1 4 3;2 1 4 3];...
    [3 4 1 2;3 4 1 2;3 4 1 2;3 4 1 2];[4 3 2 1;4 3 2 1;4 3 2 1;4 3 2 1]};

% we want to repeat the orders and assign them to four runs
img_in_order = repmat(ord_subgrp{mod(rr-1,4)+1},4,1);
quad_in_order = repmat(quad_subgrp{mod(rr-1,4)+1},4,1);

% %% Check quadrant & image orders and their pairing - comment this if you already tested it before
% temp_quad = [];
% for ii = 1:4
%     temp_quad = [temp_quad; quad_subgrp{ii}];
% end
% temp_quad = unique(temp_quad,'rows','stable');
% 
% temp_img = [];
% for ii = 1:4
%     temp_img = [temp_img; ord_subgrp{ii}];
% end
% temp_img = unique(temp_img,'rows','stable');
% 
% % test pairing between quadrant and image
% temp_pair = zeros(16,4); % column - quadrant; values - images
% temp_pair_ct = []; temp_pair_ct_cat = [];
% for ii = 1:16
%     temp_pair(ii,temp_quad(ii,:)) = temp_img(ii,:);
% end
% for ii = 1:4
% temp_pair_ct(ii,:) = histcounts(temp_pair(:,ii),0.5:4.5); % rows - quadrant; column - images
% end
% temp_pair_ct_cat = [sum(temp_pair_ct(:,1:2),2) sum(temp_pair_ct(:,3:4),2)];

%% Next thing is to determine the retrieval order

% Retrieve each order position equally
ret_ord = []; % which ordered position to retrieve
for ii = 1:4 % repeat 4 times
ret_ord = [ret_ord;ii*ones(4,1)]; % 4 quadrant/image orders
end

% retrieve image
for ii = 1:length(ret_ord)
ret_img(ii,1) = img_in_order(ii,ret_ord(ii));
end

% retrieve quadrant
for ii = 1:length(ret_ord)
ret_quad(ii,1) = quad_in_order(ii,ret_ord(ii));
end

%% Save the data first
%save([ttl '/order.mat'],'quad_in_order','img_in_order','ord_subgrp','quad_subgrp','ret_ord');

%% Now let's assign different locations

ecc_all = [3 5.5]; ecc = ecc_all(mod(rr+1,2)+1);
ang = [15:20:75];%11.25:22.5:360;
ang = [ang ang+90 ang+180 ang+270];

%%
ang_all = ang32((((mod(rr-1,4)+1)-1)*16+1) : ((mod(rr-1,4)+1)*16),:);

%% construct variables for manual shuffle
% 1. all the angles -> ang_all
% 2. all the delays -> delay
% 3. all the presentation orders:
% 3.1. quadrant order: col number -> order, col value -> quadrant
% quad_in_order: full order
% use a number to label each order: 1-24
quad_in_order_uniq = unique(quad_in_order,'rows','sorted');
for ii = 1:n
    [~, index]=ismember(quad_in_order_uniq,quad_in_order(ii,:),'rows');
    quad_in_order_lb(ii,1) = find(index);
end
% 3.2. image order: col number -> order, col value -> image
% order_final_reverse
img_in_order_uniq = unique(img_in_order,'rows','sorted');
for ii = 1:n
    [~, index]=ismember(img_in_order_uniq,img_in_order(ii,:),'rows');
    img_in_order_lb(ii,1) = find(index);
end
% % assign each image randomly: 1_1.jpg -> 4_2.jpg
% % first number -> image category: female face, male face, house, scene
% % second number -> actual image in each category, 2 per category
% for ii = 1:4 % each image
%     for jj = 1:4 % each order
%        temp = find(ret_img == ii & ret_ord == jj);
%       ret_img_act(temp,1) = randi(2,length(temp),1);
%     end
% end
% % Other than retrieval image, other images can be assigned randomly
% img_act = zeros(size(img_in_order,1),1);
% for ii = 1:n
%     img_act(ii) = ret_img_act(ii);
% end

if rr == 1|rr==3
    ret_img_act = ones(16,1);
    img_act = ret_img_act;
elseif rr == 2|rr==4
    ret_img_act = 2*ones(16,1);
    img_act = ret_img_act;
end

% 4. retrieval
% 4.1. retrieval order: ret_ord
% 4.2. retrieval quadrant: ret_quad
% 4.3. retrieval image: ret_img
% 4.4. retrieval location
for ii = 1:n
ret_ang(ii,1) = ang_all(ii,ret_quad(ii));
end

%% Next, delay & iti duration
% For each retrieved image category and each order,
% We would like the same average delay and iti duration

% Example:
% Run 1, for each quadrant order (e.g.,1432), 
% cue first, second, third and fourth -->>
% delay1, delay1, delay2, delay2; iti1, iti2, iti3, iti4

% Run 2, for each quadrant order,
% cue first, second, third and fourth -->>
% delay1, delay1, delay2, delay2; iti2, iti1, iti4, iti3

% Run 3, for each quadrant order,
% cue first, second, third and fourth -->>
% delay2, delay2, delay1, delay1; iti3, iti4, iti1, iti2

% Run 4, for each quadrant order,
% cue first, second, third and fourth -->>
% delay2, delay2, delay1, delay1; iti4, iti3, iti2, iti1

% delay duration
delay = zeros(n,1);
if rr == 1 | rr == 2 | rr ==5
delay_t_io = [delay_t(1);delay_t(1);delay_t(2);delay_t(2)];%repmat(delay_t,2,1);
elseif rr == 3 | rr == 4 | rr == 6
    delay_t_io = [delay_t(2);delay_t(2);delay_t(1);delay_t(1)];%repmat(delay_t,2,1);
end
switch mod(rr-1,4)+1
    case 1
        iti_t_io = ITI;
    case 2
        iti_t_io = [ITI(2);ITI(1);ITI(4);ITI(3)];
    case 3
        iti_t_io = [ITI(3);ITI(4);ITI(1);ITI(2)];
    case 4
        iti_t_io = [ITI(4);ITI(3);ITI(2);ITI(1)];
end

iti = zeros(n,1);
iti_t_io = ITI;
for jj = 1:4 % each presentation quadrant order
    for kk = 1:4 % each retrieval order
       temp = find(quad_in_order_lb == jj & ret_ord == kk); % always 4 trials % ret_ord
       delay(temp,1) = delay_t_io(kk);%Shuffle(delay_t_io);
       iti(temp,1) = iti_t_io(kk);%Shuffle(iti_t_io);
end
end
%% manual shuffle

% before shuffling, let's load our manually set-up order of the events

%         if ses <=4
%     data = importfile(['../optseq1028/run',num2str(rr+4*(ses-1)),'-001.par']);
%         elseif ses == 5
%     data = importfile(['../optseq1028/run',num2str(4),'-001.par']);
%         end
%         
%     optord = data{:,2}; % order of the events (not quadrant or image order, but trial order)
%     optord(optord==0) = []; % get rid of all null events (ord==0)
  

switch rr
    case 1
        optord = [1
2
3
4
2
1
4
3
4
3
1
2
3
4
2
1];
        cueod = [1
2
3
4
3
2
1
4
2
1
3
4
2
3
1
4];
    case 2
        optord = [3
1
4
2
1
3
2
4
2
4
3
1
4
2
1
3];
                cueod = [2
4
1
3
2
1
4
3
1
2
3
1
4
2
3
4];
    case 3
        optord = [2
4
1
3
4
2
3
1
3
1
2
4
1
3
4
2];
                cueod = [4
3
1
2
1
3
3
4
1
3
2
4
2
4
2
1];
    case 4
        optord = [4
3
2
1
3
4
1
2
1
2
4
3
2
1
3
4];
         cueod = [3
4
1
2
1
4
1
3
4
2
1
2
4
3
3
2];

end
    
%n = size(data,1); % how many trials in total
rep_ind = zeros(10,1); % repetitive trial index
% ang_all(randperm(16),:) 
data = [];
data_left = [];
df = [];

while ~isempty(rep_ind)
    data = [delay img_in_order_lb quad_in_order_lb ret_img ret_quad ret_ord img_act iti]; % concatenate all data together
    data_sff = [];
    rep_ind = [];
    
    %data_left = data; % data left to choose

            % we want to create four subtables, each table correspond to the trials
        % belonging to a specific quadrant order (1423,...)
        for ii = 1:4 % for each quadrant order
            for jj = 1:4 % for each cue
            data_left{ii,jj} = data(find(data(:,3)==ii & data(:,6) == jj),:); % each table should contain 16 rows - 16 trials
            end
        end

    % for trial from 1 to n
    for ii = 1:n

        % when it is the first trial of that run
        if ii == 1
                % randomly select a trial from the corresponding pool
                % data_left{ord(1)}
                df = data_left{optord(ii),cueod(ii)}; % current data_left for that trial
                tr = 1; tr_1 = randi(size(df,1)); % randomly select the first trial
                data_sff(ii,:) = df(tr_1,:); % assign tr_1th trial to the first trial of that run
                df(tr_1,:) = []; % remove that trial from data
                data_left{optord(ii),cueod(ii)} = df; % assign df back to data_left

        else % if it is not the first trial
                % find a trial with a different location from the previous trial
                % from the corresponding pool data_left{ord(ii)}
                df = data_left{optord(ii),cueod(ii)}; % current data_left for that trial
                
%                 % first calculate difference between the rest of the trial and the
%                 % previous trial
%                 % we are making sure that there is no repetitive polar angle of the dots                
%                 dif = all(df(:,[1:4]) - data_sff(ii-1,[1:4]),2); % if all the numbers are different (all the difference numbers should be non-zero)
%             ind = find(dif==1); % if all the numbers are different (find dif == 1)
% 
%             if isempty(ind) % if no trial left that has no repetition
%                 % proceed with the rest of the trial
%                 ind_rand = randi(size(data_left,1));
%                 % record the repetitive trial number
%                 rep_ind = [rep_ind ii];
%             else
%                 ind_rand = ind(randi(length(ind))); % pick one trial from those
%             end
ind_rand = 1;
                data_sff(ii,:) = df(ind_rand,:); % assign the trial to data_sff
                df(ind_rand,:) = []; % remove that trial from the data_left
                data_left{optord(ii),cueod(ii)} = df; % assign df back to data_left
        end
    end
end

%% add the locations
data = [];
data_left = [];
df = [];
%n = size(data,1); % how many trials in total
rep_ind = zeros(10,1); % repetitive trial index
% ang_all(randperm(16),:) 
while ~isempty(rep_ind)
    data = ang_all(randperm(16),:);
data_left = data;
    rep_ind = [];
    ang_sff = [];
    for ii = 1:n
        if ii == 1
                tr_1 = randi(size(data_left,1)); % randomly select the first trial
                ang_sff(ii,:) = data_left(tr_1,:); % assign tr_1th trial to the first trial of that run
                data_left(tr_1,:) = []; % remove that trial from data
        else
            dif = all(data_left(:,[1:4]) - ang_sff(ii-1,[1:4]),2); % if all the numbers are different (all the difference numbers should be non-zero)
            ind = find(dif==1); % if all the numbers are different (find dif == 1)

            if isempty(ind) % if no trial left that has no repetition
                % proceed with the rest of the trial
                ind_rand = randi(size(data_left,1));
                % record the repetitive trial number
                rep_ind = [rep_ind ii];
            else
                ind_rand = ind(randi(length(ind))); % pick one trial from those
            end
                ang_sff(ii,:) = data_left(ind_rand,:); % assign the trial to data_sff
                data_left(ind_rand,:) = []; % remove that trial from the data_left
        end

        end
end

%%
    data_sff_new = [ang_sff data_sff];
    data_sff = data_sff_new;

%% Construct the whole trial

for ii = 1:size(data_sff,1)
    temp = quad_in_order_uniq(data_sff(ii,7),:); % full order for that trial
    data_sff(ii,1:4) = [data_sff(ii,temp(1)),data_sff(ii,temp(2)),...
        data_sff(ii,temp(3)),data_sff(ii,temp(4))];
end

%% write to an excel file to double-check the shuffled results
%writematrix(data_sff,[ttl '/' ttl '.csv'])

%% Construct the full experiment
run = [rr*ones(16,1)];
trial = [1:n]';

% quadrant order (1-4, first-> 4th quadrant)
quad_ord = data_sff(:,7); % order condition (1-24)
quad_ord_full = quad_in_order_uniq(data_sff(:,7),:); % full quadrant order

% image category order (1-4)
img_cat_ord = data_sff(:,6);
img_cat_ord_full = img_in_order_uniq(data_sff(:,6),:); % full image order

% image sub-category (1-8, 8 different images under each category)
img_subcat = data_sff(:,11);

% actual image name
for ii = 1:n % trial
    for jj = 1:4 % for each order
img_name{ii,jj} = [num2str(img_cat_ord_full(ii,jj)) '_' num2str(img_subcat(ii)) '.png'];
    end
end

% fixation location
fix_x = x_res/2*ones(n,1);
fix_y = y_res/2*ones(n,1);

% time
iti_t_all = data_sff(:,end);%Shuffle([ITI(1)*ones(6,1);ITI(2)*ones(6,1);ITI(3)*ones(4,1)]); %iti_t_all(1) = 5.5; % datasample(ITI,n)
fix_t_all = fix_t*ones(n,1) + iti_t_all;
stim_t_all = fix_t_all + repmat([stim_t,stim_t*2+blank_t,stim_t*3+blank_t*2,stim_t*4+blank_t*3],n,1);
blank_t_all = fix_t_all + repmat([stim_t+blank_t,stim_t*2+blank_t*2,stim_t*3+blank_t*3],n,1);
delay_t_all = data_sff(:,5)+stim_t_all(:,4);
cue_t_all = delay_t_all + cue_t*ones(n,1);
resp_t_all = cue_t_all + resp_t*ones(n,1);

% iti_t_run = repmat(ITI,3,1); 
% iti_t_all = [Shuffle(iti_t_run);Shuffle(iti_t_run);Shuffle(iti_t_run);Shuffle(iti_t_run)];
% srt_t_all = srt_t*ones(n,1);

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
%feed_color_all = repmat(feed_color,n,1);

% targets actual location
% % what type of jitter do I want?
% % each quadrant will be repeated 16 times in a run
% % each location will be repeated 4 times in a run (there are 4 locations in
% % each quadrant)
% % note that here location = polar angle, because the eccentricity is fixed
% % in one run
% % here, what we want to do -- for each location, we want to psuedo-randomly jitter
% % the eccentricity and the polar angle and make sure that the mean is at
% % the center
% ang_run = data_sff(:,1:4); ang_rand = zeros(16,4); ecc_rand = zeros(16,4);
% for kk = ang
%     ind = find(ang_run==kk);
%     jit = linspace(-1,1,4); jit = jit(randperm(4))';
%     ang_rand(ind) = ang_run(ind) + jit_ang*jit;
%     jit = linspace(-1,1,4); jit = jit(randperm(4))';
%     ecc_rand(ind) = ecc + jit_ecc*jit;
% end

jit = linspace(-1,1,64); 
jit = jit(randperm(64)); jit = reshape(jit,n,4);
ecc_rand = ecc + jit_ecc*jit;%(2*rand(n,4)-1);

jit = linspace(-1,1,64); 
jit = jit(randperm(64)); jit = reshape(jit,n,4);
ang_rand = data_sff(:,1:4) + jit_ang*jit;%(2*rand(n,4)-1);

x = tand(ecc_rand).*d.*cosd(ang_rand)./w.*x_res + x_res/2;
y = y_res/2-tand(ecc_rand).*d.*sind(ang_rand)./h.*y_res;

% %% plot the jittering results:
% figure(1);clf
% scatter(reshape(x,1,[]),reshape(y,1,[]),10)
% xlim([0 1600])
% ylim([0 900])
%%
% retrieve order, image category and quadrant
ret_ord = data_sff(:,10);
ret_img = data_sff(:,8);
ret_quad = data_sff(:,9);

% retrieve location
% for report one
for ii = 1:n
    ret_ecc(ii,1) = ecc_rand(ii,ret_ord(ii));
    ret_ang(ii,1) = ang_rand(ii,ret_ord(ii));
end
ret_x = tand(ret_ecc).*d.*cosd(ret_ang)./w.*x_res + x_res/2;
ret_y = y_res/2-tand(ret_ecc).*d.*sind(ret_ang)./h.*y_res;

% cue
% cue image
% cue quadrant (also an image)
% cue order (an image)
for ii = 1:n % trial
cue_img{ii,1} = [num2str(ret_img(ii)) '_' num2str(img_subcat(ii)) '.png'];
cue_quad{ii,1} = [num2str(ret_quad(ii)) '.png'];
cue_ord{ii,1} = [num2str(ret_ord(ii)) '.png'];
end

% full table
mltable = table(run,trial,quad_ord,...
    quad_ord_full(:,1),quad_ord_full(:,2),quad_ord_full(:,3),quad_ord_full(:,4),...
    img_cat_ord, img_cat_ord_full(:,1),img_cat_ord_full(:,2),img_cat_ord_full(:,3),img_cat_ord_full(:,4),...
    img_subcat,img_name(:,1),img_name(:,2),img_name(:,3),img_name(:,4),...
    ecc_rand(:,1),ecc_rand(:,2),ecc_rand(:,3),ecc_rand(:,4),...
    ang_rand(:,1),ang_rand(:,2),ang_rand(:,3),ang_rand(:,4),...
    x(:,1),x(:,2),x(:,3),x(:,4),...
    y(:,1),y(:,2),y(:,3),y(:,4),...
    ret_ord,ret_quad,ret_img,cue_img,cue_quad,cue_ord, ...
    ret_ecc,ret_ang,ret_x,ret_y,fix_x,fix_y,...
    fix_t_all,stim_t_all(:,1),blank_t_all(:,1),stim_t_all(:,2),blank_t_all(:,2),...
    stim_t_all(:,3),blank_t_all(:,3),stim_t_all(:,4),...
    delay_t_all,cue_t_all,resp_t_all,iti_t_all,fix_size_pix_all,dot_size_pix_all,roi_size_pix_all,...
    data_sff(:,1),data_sff(:,2),data_sff(:,3),data_sff(:,4),ecc*ones(n,1),...
    'VariableNames',{'run_ind','trial_ind','quad_order_ind',...
    'stim1_quadrant','stim2_quadrant','stim3_quadrant','stim4_quadrant',...
    'img_cat_ind','stim1_cat','stim2_cat','stim3_cat','stim4_cat',...
    'img_subcat_ind','stim1_img','stim2_img','stim3_img','stim4_img',...
    'stim1_ecc','stim2_ecc','stim3_ecc','stim4_ecc',...
    'stim1_ang','stim2_ang','stim3_ang','stim4_ang',...
    'stim1_x','stim2_x','stim3_x','stim4_x',...
    'stim1_y','stim2_y','stim3_y','stim4_y',...
    'probe_order','probe_quad','probe_img_category','cue_img','cue_quad','cue_ord', ...
    'probe_ecc','probe_ang','probe_x','probe_y','fix_x','fix_y',...
    'fix_t','stim1_t','blank1_t','stim2_t','blank2_t','stim3_t','blank3_t','stim4_t',...
    'delay_t','cue_t','resp_t','iti_t','fix_size','dot_size','roi_size',...
    'stim1_ang_orig','stim2_ang_orig','stim3_ang_orig','stim4_ang_orig','ecc'});

% mltable(2:17,:) = mltable(1:16,:);
% mltable(18,:) = mltable(17,:);
% mltable.trial_ind = [1:18]';

%     mltable(17,:) = mltable(16,:);
%     mltable.trial_ind = [1:17]';
%     mltable.iti_t(17) = 14.5;
%     mltable{1,46:56} = mltable{1,46:56}-mltable.iti_t(1)+5;
%     temp.iti_t(1) = 5; 

%% change iti and other timing
        
        % copy the 16th trial to the 17th
        mltable(17,:) = mltable(16,:); mltable.trial_ind(17) = 17;
          
        % change the iti of the first trial (baseline) and subsequent timing of
        % all events in the first trial
        mltable{1,46:56} = mltable{1,46:56}-mltable.iti_t(1)+5;
        mltable.iti_t(1) = 5;
        
        % calculate how long the run is
        tm = mltable(:,[46:57]); tm = table2array(tm);
        tm_all = sum(tm(1:16,11));
        
        % set up the iti of the last trial(baseline) such that the duration
        % of the run is 419
        mltable.iti_t(17) = 419 - tm_all;


%% Calculate measurements
%tm = mltable(2:17,[45:53,55:58]);
tm = mltable(:,[46:54,56:57]); tm = table2array(tm); 
tm_all = sum(tm(1:16,10))+tm(17,end);
tm_all_min = tm_all/60;
n_vol = tm_all/tr
n_vol_runs(rr) = n_vol;

%  % calculate measurements
% tm = temp(:,[46:57]); tm = table2array(tm); 
% tm_all = sum(tm(1:16,11))+tm(17,end);
% tm_all_min = tm_all/60;
% n_vol = tm_all/1
% n_vol_runs(rr) = n_vol;

%% save the table
% Write into a csv file
% writetable(mltable,[ttl '/' ttl '_all.csv'])
writetable(mltable,[ttl '.csv'])

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
%saveas(figure(3),[ttl '.jpg'])

end
%writematrix(n_vol_runs,'sMGS_TR1000_n_vol.csv')
end