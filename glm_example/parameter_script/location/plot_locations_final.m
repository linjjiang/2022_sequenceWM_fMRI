%%
clear all
close all
clc
addpath(genpath('C:\Users\linjing\Desktop\VSWM-DA\datasource_script\'))

%% prototype locations
% locations
d = 138.9; % cm
h = 35;
w = 43;
x_res = 1024;
y_res = 768;

% maximum range of half of the screen screen
w_dva = atand(w/2/d) % horizontal half
h_dva = atand(h/2/d) % vertical half

% pix
ecc_all = [3 5.5];
ang = [15:20:75];%11.25:22.5:360; %[22.5:15:67.5]
ang = [ang ang+90 ang+180 ang+270];

X = [ecc_all(1)*cosd(ang) ecc_all(2)*cosd(ang)];
Y = [ecc_all(1)*sind(ang) ecc_all(2)*sind(ang)];

figure(1);clf
scatter(X,Y,100,'LineWidth',3)
xlim([-w_dva w_dva])
ylim([-h_dva h_dva])
line([0 0],[-h_dva h_dva],'LineStyle','--','LineWidth',3)
line([-w_dva w_dva],[0 0],'LineStyle','--','LineWidth',3);
set(gca,'PlotBoxAspectRatio',[w h 1])
axis off
saveas(figure(1),'32_locations.jpg')

%% How many groups of 4 locations are possible

% construct different combinations
data = combntns(ang,4);

% get rid of angles in the same quadrant
quad = zeros(size(data,1),1);
for ii = 1:size(data,1)
    v = quadrant(data(ii,:));
    if length(v) == length(unique(v)) % if in different quadrants
        quad(ii) = 0;
    else
        quad(ii) = 1;
    end
end
data_quad = data(quad==0,:);

% include fully asymmetric + only one axial symmetric + only one central symmetric case
s = zeros(size(data_quad,1),1);
for ii = 1:size(data_quad,1)
    s(ii) = symmetric5(data_quad(ii,:));
end
data_quad_sym = data_quad(s==0,:);

% If you want to check symmetric type, uncommend the following section
% check symmetric type
% include fully asymmetric + only one axial symmetric + only one central symmetric case
st = zeros(size(data_quad,1),1);
for ii = 1:size(data_quad,1)
    st(ii) = symmetric6(data_quad(ii,:));
end
% plot different types of location combinations
% st == 1: axial symmetric
% st == 2: central symmetry
% st == 3: fully assymetric
ind = find(st==1); ind = ind(1);
x = ecc_all(2)*cosd(data_quad(ind,:));
y = ecc_all(2)*sind(data_quad(ind,:));
X = ecc_all(2)*cosd(ang);
Y = ecc_all(2)*sind(ang);
figure(1);clf
scatter(X,Y,200,'LineWidth',3)
hold on
scatter(x,y,200,'filled','LineWidth',3)
xlim([-w_dva w_dva])
ylim([-h_dva h_dva])
line([0 0],[-h_dva h_dva],'LineStyle','--','LineWidth',3)
line([-w_dva w_dva],[0 0],'LineStyle','--','LineWidth',3);
set(gca,'PlotBoxAspectRatio',[w h 1])
axis off


% % further define fully asymmetric dataset
% s = zeros(size(data_quad_sym,1),1);
% for ii = 1:size(data_quad_sym,1)
%     s(ii) = symmetric3(data_quad_sym(ii,:));
% end
% data_quad_sym_sym = data_quad_sym(s==0,:);

% % get rid of location sets with all angles > 45 degress apart (70 degrees
% % apart)
% temp = [data_quad_sym data_quad_sym(:,1)+360];
% dif = diff(temp'); dif = dif';
% ind = all(dif>50,2);
% data_quad_sym_apt = data_quad_sym(~ind,:);
% 
% % get rid of all location sets that are too close (10 degrees) to the
% % axes
% s = zeros(size(data_quad_sym_apt,1),1);
% for ii = 1:size(data_quad_sym_apt,1)
%     s(ii) = near_axes(data_quad_sym_apt(ii,:),10);
% end
% data_quad_sym_apt_far = data_quad_sym_apt(s==0,:);

% %% Then I want to pick 128 out of 168 groups that are most dissimilar
% % And then I will divide them into 2 dissimilar groups

%% Dot simulation animation (uncommend the following section if you do not need animation anymore)

load set1.mat
X = ecc_all(2)*cosd(ang);
Y = ecc_all(2)*sind(ang);
figure(1);clf
upper = ang_set(1:16,1:2);
lower = ang_set(1:16,3:4);
% temp = [];
% for ii = 1:4 % four locations of the first quadrant
%     for jj = 1:4 % four locations of the second quadrant
%         temp = [temp;ang(ii) ang(jj+4)];
%     end
% end

% plot all 16 locations
for ii = 1:16 % 16 different combinations of two upper dots
    subplot(4,4,ii)
    scatter(X,Y,50,'LineWidth',1)
    xlim([-w_dva w_dva])
ylim([-h_dva h_dva])
line([0 0],[-h_dva h_dva],'LineStyle','--','LineWidth',3)
line([-w_dva w_dva],[0 0],'LineStyle','--','LineWidth',3);
%set(gca,'PlotBoxAspectRatio',[w h 1])
    axis off
box off
end

% plot 16 sets of upper 2 dots
for ii = 1:16
    subplot(4,4,ii)
x = ecc_all(2)*cosd(upper(ii,:));
y = ecc_all(2)*sind(upper(ii,:));
hold on
pause(0.5)
scatter(x,y,50,[0.8500 0.3250 0.0980],'filled','LineWidth',1)
print(['Frame ' num2str(ii)], '-dpng', '-r150');
end

% plot 16 sets of lower 2 dots
for ii = 1:16
    subplot(4,4,ii)
x = ecc_all(2)*cosd(lower(ii,:));
y = ecc_all(2)*sind(lower(ii,:));
hold on
pause(0.5)
scatter(x,y,50,[0.4940 0.1840 0.5560],'filled','LineWidth',1)

print(['Frame ' num2str(ii+16)], '-dpng', '-r150');

end

% stack images into gif animation
GifName = 'loc_select_animation.gif';
delay = 0.5;    % Delay between frames (s)
for ii = 1:32
    [A, ~] = imread(['Frame ' num2str(ii) '.png']);
    [X, map] = rgb2ind(A, 256);
    if ii == 1
        imwrite(X, map, GifName, 'gif', 'LoopCount', inf, 'DelayTime', delay)
    else
        imwrite(X, map, GifName, 'gif', 'WriteMode', 'append', 'DelayTime', delay)
    end
end

% delete intermediate files
delete Frame*


%% Check the dissimilar index
load 1000_set.mat
temp = [];
% check the uniqueness of the locations
for ii = 1:1000
    temp(ii,:) = reshape(q(ii).q1234,1,[]);
end
uniq_temp = unique(temp,'rows');

% check min, max, mean and the largest 2 dissimilar indexes
m = mean(d_mean_srt,2);
min(m)
max(m)
mean(m)
uniq_m = unique(m);
uniq_m(end)
uniq_m(end-1)


%% I simulated 1000 groups of 4 locations - see the last half of the script
% Re-select the sets

% copy "1000_set.mat" to the folder and load the mat file
load 1000_set.mat

% ind_new: index of different sets of locations, sorted by disimilarity
% each row has two sets of locations (two runs), each has 16 groups of 4 dots
% these two sets of locations do not overlap - in total, we have 32
% different groups of 4 dots
%ind_best = ind_new(end,:)
% locs: we want to pick 8 rows of ind_new, such that we have in total
% 32*8 = 256 trials
d_mean_uniq = unique(d_mean_srt,'rows'); % unique values of disimilarity
locs = [0 0];
for jj = [1000 999 998 997 996 995 994 993 992 991]
ind = find(ismember(d_mean_srt,d_mean_uniq(jj,:),'rows')==1); % find rows with that unique values of disimilarity
rows = ind_new(ind,:) % find location index of those rows
rows_uniq = unique(rows(:,1));
for ii = 1:length(rows_uniq)
    locs2 = rows(rows(:,1)==rows_uniq(ii),2);
    locs2_ind = ~ismember(locs2,locs(:,2));
    locs = [locs;rows_uniq(ii),locs2(locs2_ind(randi(length(locs2_ind))))];
end
end
%% check how many patterns overlap in the 2 datasets with highest dissimilar
% indexes

intersect(q(locs(2,1)).q1234,q(locs(2,2)).q1234,'rows')
% the output of the intersect should be empty, meaning that there are no
% overlapping rows/groups of lcoations

% combine 2 datasets into 1 larger datasets containing 32 non-overlapping
% locations
comb_q1 = [q(locs(2,1)).q1234;q(locs(2,2)).q1234]; % 2 datasets with the highest dissimilarity
comb_q2 = [q(locs(3,1)).q1234;q(locs(3,2)).q1234]; % 2 datasets with the second highest dissimilarity
% check how many patterns are overlapped
overlap = intersect(comb_q1,comb_q2,'rows')
num_nonoverlap = 64 - size(overlap,1)

%% Plot each set in a more pretty way

d = 138.9; % cm
h = 35;
w = 43;
x_res = 1024;
y_res = 768;

% maximum range of half of the screen screen
w_dva = atand(w/2/d) % horizontal half
h_dva = atand(h/2/d) % vertical half

% pix
ecc_all = [3 5.5];
ang = [15:20:75];%11.25:22.5:360; %[22.5:15:67.5]
ang = [ang ang+90 ang+180 ang+270];


figure(1);clf

ang_set1 = load('set1.mat');
ang_set2 = load('set2.mat');
ang_set = [ang_set1.ang_set;ang_set2.ang_set];
% plot all locations
for ii = 1:64 % 16 different combinations of four dots
    subplot(4,16,ii)
    X = ecc_all(2)*cosd(ang);
Y = ecc_all(2)*sind(ang);
    scatter(X,Y,20,'LineWidth',1)
    
    x = ecc_all(2)*cosd(ang_set(ii,:));
    y = ecc_all(2)*sind(ang_set(ii,:));
hold on
scatter(x,y,20,[0.8500 0.3250 0.0980],'filled','LineWidth',1)

    xlim([-w_dva w_dva])
ylim([-h_dva h_dva])
line([0 0],[-h_dva h_dva],'LineStyle','--','LineWidth',1)
line([-w_dva w_dva],[0 0],'LineStyle','--','LineWidth',1);
%set(gca,'PlotBoxAspectRatio',[w h 1])
    axis off
box off
end
saveas(figure(1),'set1_2.jpg')
%% set 2
figure(1);clf
% plot all locations
for ii = 33:64 % 16 different combinations of four dots
    subplot(4,8,ii-32)
    X = ecc_all(2)*cosd(ang);
Y = ecc_all(2)*sind(ang);
    scatter(X,Y,50,'LineWidth',1)
    
    x = ecc_all(2)*cosd(ang_set(ii,:));
    y = ecc_all(2)*sind(ang_set(ii,:));
hold on
scatter(x,y,50,[0.8500 0.3250 0.0980],'filled','LineWidth',1)

    xlim([-w_dva w_dva])
ylim([-h_dva h_dva])
line([0 0],[-h_dva h_dva],'LineStyle','--','LineWidth',1)
line([-w_dva w_dva],[0 0],'LineStyle','--','LineWidth',1);
%set(gca,'PlotBoxAspectRatio',[w h 1])
    axis off
box off
end
saveas(figure(1),'set2.jpg')



%% how many unique groups of locations in total
locs_all = [];
for ii = 2:11 % 9 sets of locations
    for jj = 1:2
locs_all = [locs_all;q(locs(ii,jj)).q1234];
    end
end
locs_uniq = unique(locs_all,'rows'); length(locs_uniq)
sum(locs_all==1)
%% save those sets
for kk = 2:11 % 9 sets of locations
% take a look at those locations
x_res = 1024; y_res = 768;
    figure(1);clf
    %set(gcf,'position',[0,0,2000,600])
    ang1234 = [q(locs(kk,1)).ang1234;q(locs(kk,2)).ang1234];
%         temp1 = ang1234(32,3:4);
%     temp2 = ang1234(19,3:4);
%     ang1234(32,3:4) = temp2; ang1234(19,3:4) = temp1;

for ii = 1:32
    subplot(4,8,ii)

temp1 = ang1234(ii,:);
temp_x = x_res/2+300.*cosd(temp1);
temp_y = y_res/2+300.*sind(temp1);
scatter(temp_x,temp_y,15)
xlim([0 x_res])
ylim([0 y_res])
hold on
scatter(x_res/2,y_res/2,15,'+')
set(gca,'PlotBoxAspectRatio',[w h 1])

%line([x_res/2 x_res/2],[0 y_res])
%line([0 x_res],[y_res/2,y_res/2]);
title(num2str(ii))
end
ang_set = ang1234;

save(['set',num2str(kk-1),'.mat'],'ang_set')
saveas(figure(1),['set',num2str(kk-1),'.jpg'])
end


%% Prototype of 4 locations
clear q
% str = {'set1','set2','set3','set4','set5','set6','set7','set8','set9','set10',...
%     'set11','set12','set13','set14','set15','set16','set17','set18','set19','set20'};
%str = {'set3','set4'};
for kk = 1:100%1:4 % for each set of randomization
%    mkdir(str{kk})
%    cd(str{kk})
q1 = [1:4];q2 = [5:8];
q12 = [];

for ii = 1:4
    for jj = 1:4
        q12 = [q12;q1(ii) q2(jj)];
    end
end

q34 = q12+8;
s = 1;
while any(s) 
q1234 = [q12 q34(randperm(16),:)];
ang1234 = ang(q1234);
for ii = 1:size(q1234,1)
    s(ii) = symmetric4(ang1234(ii,:)); % s = 0, asymetrical
end
end

% for ii = 1:16
%     %subplot(4,8,ii)
%     figure(ii);clf
%     set(gcf,'position',[0,0,x_res,y_res])
% temp1 = ang1234(ii,:);
% temp_x = x_res/2+200.*cosd(temp1);
% temp_y = y_res/2+200.*sind(temp1);
% scatter(temp_x,temp_y,300)
% xlim([0 x_res])
% ylim([0 y_res])
% line([x_res/2 x_res/2],[0 y_res])
% line([0 x_res],[y_res/2,y_res/2]);
% title(num2str(ii))
% saveas(figure(ii),[num2str(ii),'.jpg'])
% end
q(kk).q1234 = q1234;
q(kk).ang1234 = ang1234;
%save([str{kk},'.mat'],'ang1234','q1234')
%cd ..
end

%%
% str = {'set1','set2','set3','set4','set5','set6','set7','set8','set9','set10',...
%     'set11','set12','set13','set14','set15','set16','set17','set18','set19','set20'};
% 
% for ii = 1:20
%     cd(str{ii})
% q(ii) = load([str{ii},'.mat']);
% cd ..
% end

%%
save('1000_set.mat','q')


%% overlap
ind = [];
for ii = 1:1000
    for jj = 1:1000
ist{ii,jj} = intersect(q(ii).q1234,q(jj).q1234,'rows');
if isempty(ist{ii,jj}) %size(ist{ii,jj},1) <= 1
    ind = [ind;ii jj];
end
    end
end

%%
save('1000_set.mat','q','ist','ind')

%% procrustes similarity analysis
d_mean = []; d = [];
for kk = 1:size(ind,1)
    ang_ii = q(ind(kk,1)).ang1234;
    ang_jj = q(ind(kk,1)).ang1234;
    d_wt = [];
    d_bt = [];
    for hh = 1:16
        X_ii(hh,:) = 200.*cosd(ang_ii(hh,:));
        Y_ii(hh,:) = 200.*sind(ang_ii(hh,:));
        
        X_jj(hh,:) = 200.*cosd(ang_jj(hh,:));
        Y_jj(hh,:) = 200.*sind(ang_jj(hh,:));
    end
    
    for hh = 1:16
        tp = 1:16; tp(hh) = [];
        for ll = tp
     [temp,~,~] = procrustes([X_ii(hh,:)' Y_ii(hh,:)'],[X_ii(ll,:)' Y_ii(ll,:)']); d_wt = [d_wt;temp];
        end
    end
    
        for hh = 1:16
        tp = 1:16; tp(hh) = [];
        for ll = tp
     [temp,~,~] = procrustes([X_jj(hh,:)' Y_jj(hh,:)'],[X_jj(ll,:)' Y_jj(ll,:)']); d_wt = [d_wt;temp];
        end
    end
    
        for hh = 1:16
        for ll = 1:16
     [temp,~,~] = procrustes([X_ii(hh,:)' Y_ii(hh,:)'],[X_jj(ll,:)' Y_jj(ll,:)']); d_bt = [d_bt;temp];
        end
        end
    d(kk).d_wt = d_wt; d(kk).d_bt = d_bt;
    d_mean(kk,1) = mean(d_wt); d_mean(kk,2) = mean(d_bt);
end
%[d,Z,tr] = procrustes(X,Y);

%%
save('1000_set.mat','q','ist','ind','d_mean')

%%
[d_mean_srt,ind_srt] = sortrows(d_mean);
ind_new = ind(ind_srt,:);
save('1000_set.mat','q','ist','ind','d_mean_srt','ind_new','d_mean')


%%
ind_best = ind_new(end,:)
intersect(q(ind_best(1)).q1234,q(ind_best(2)).q1234,'rows')

%%
ind_best = ind_new(end-1,:)
intersect(q(ind_best(1)).q1234,q(ind_best(2)).q1234,'rows')

%%
%save('final_set.mat')
intersect(q(370).q1234,q(997).q1234,'rows')

%%
x_res = 1600; y_res = 900;
    figure(1);clf
    set(gcf,'position',[0,0,1600,600])
    ang1234 = [q(370).ang1234;q(997).ang1234];
        temp1 = ang1234(32,3:4);
    temp2 = ang1234(19,3:4);
    ang1234(32,3:4) = temp2; ang1234(19,3:4) = temp1;

for ii = 1:32
    subplot(4,8,ii)

temp1 = ang1234(ii,:);
temp_x = x_res/2+300.*cosd(temp1);
temp_y = y_res/2+300.*sind(temp1);
scatter(temp_x,temp_y,15)
xlim([0 x_res])
ylim([0 y_res])
hold on
scatter(x_res/2,y_res/2,15,'+')
%line([x_res/2 x_res/2],[0 y_res])
%line([0 x_res],[y_res/2,y_res/2]);
title(num2str(ii))
end
ang_set1 = ang1234;
save('set1.mat','ang_set1')
saveas(figure(1),'set1.jpg')

%%
    figure(1);clf
    set(gcf,'position',[0,0,1600,600])
    ang1234 = [q(328).ang1234;q(543).ang1234];
            temp1 = ang1234(2,3:4);
    temp2 = ang1234(11,3:4);
    ang1234(2,3:4) = temp2; ang1234(11,3:4) = temp1;

                temp1 = ang1234(17,3:4);
    temp2 = ang1234(32,3:4);
    ang1234(17,3:4) = temp2; ang1234(32,3:4) = temp1;

for ii = 1:32
    subplot(4,8,ii)

temp1 = ang1234(ii,:);
temp_x = x_res/2+300.*cosd(temp1);
temp_y = y_res/2+300.*sind(temp1);
scatter(temp_x,temp_y,15)
xlim([0 x_res])
ylim([0 y_res])
hold on
scatter(x_res/2,y_res/2,15,'+')
%line([x_res/2 x_res/2],[0 y_res])
%line([0 x_res],[y_res/2,y_res/2]);
title(num2str(ii))
end
ang_set2 = ang1234;
save('set2.mat','ang_set2')
saveas(figure(1),'set2.jpg')

%% 
intersect([q(328).ang1234;q(543).ang1234],[q(370).ang1234;q(997).ang1234],'rows');

%% set 3
x_res = 1600; y_res = 900;
    figure(1);clf
    set(gcf,'position',[0,0,1600,600])
    ang1234 = [q(370).ang1234;q(996).ang1234];
%         temp1 = ang1234(32,3:4);
%     temp2 = ang1234(19,3:4);
%     ang1234(32,3:4) = temp2; ang1234(19,3:4) = temp1;
% 
for ii = 1:32
    subplot(4,8,ii)

temp1 = ang1234(ii,:);
temp_x = x_res/2+300.*cosd(temp1);
temp_y = y_res/2+300.*sind(temp1);
scatter(temp_x,temp_y,15)
xlim([0 x_res])
ylim([0 y_res])
hold on
scatter(x_res/2,y_res/2,15,'+')
%line([x_res/2 x_res/2],[0 y_res])
%line([0 x_res],[y_res/2,y_res/2]);
title(num2str(ii))
end
ang_set3 = ang1234;
save('set3.mat','ang_set3')
saveas(figure(1),'set3.jpg')

%% set 4
    figure(1);clf
    set(gcf,'position',[0,0,1600,600])
    ang1234 = [q(328).ang1234;q(990).ang1234];
            temp1 = ang1234(2,3:4);
    temp2 = ang1234(11,3:4);
    ang1234(2,3:4) = temp2; ang1234(11,3:4) = temp1;

%                 temp1 = ang1234(17,3:4);
%     temp2 = ang1234(32,3:4);
%     ang1234(17,3:4) = temp2; ang1234(32,3:4) = temp1;

for ii = 1:32
    subplot(4,8,ii)

temp1 = ang1234(ii,:);
temp_x = x_res/2+300.*cosd(temp1);
temp_y = y_res/2+300.*sind(temp1);
scatter(temp_x,temp_y,15)
xlim([0 x_res])
ylim([0 y_res])
hold on
scatter(x_res/2,y_res/2,15,'+')
%line([x_res/2 x_res/2],[0 y_res])
%line([0 x_res],[y_res/2,y_res/2]);
title(num2str(ii))
end
ang_set4 = ang1234;
save('set4.mat','ang_set4')
saveas(figure(1),'set4.jpg')

%% set 5
x_res = 1600; y_res = 900;

    figure(1);clf
    set(gcf,'position',[0,0,1600,600])
    ang1234 = [q(90).ang1234;q(893).ang1234];
            temp1 = ang1234(2,3:4);
    temp2 = ang1234(11,3:4);
    ang1234(2,3:4) = temp2; ang1234(11,3:4) = temp1;

%                 temp1 = ang1234(17,3:4);
%     temp2 = ang1234(32,3:4);
%     ang1234(17,3:4) = temp2; ang1234(32,3:4) = temp1;

for ii = 1:32
    subplot(4,8,ii)

temp1 = ang1234(ii,:);
temp_x = x_res/2+300.*cosd(temp1);
temp_y = y_res/2+300.*sind(temp1);
scatter(temp_x,temp_y,15)
xlim([0 x_res])
ylim([0 y_res])
hold on
scatter(x_res/2,y_res/2,15,'+')
%line([x_res/2 x_res/2],[0 y_res])
%line([0 x_res],[y_res/2,y_res/2]);
title(num2str(ii))
end
ang_set5 = ang1234;
save('set5.mat','ang_set5')
saveas(figure(1),'set5.jpg')
% 
% 
% %% Manually shuffle some of the pairs
% %     ang1234 = [q(20).ang1234;q(15).ang1234];   
% %     temp1 = ang1234(2,3:4);
% %     temp2 = ang1234(11,3:4);
% %     ang1234(2,3:4) = temp2; ang1234(11,3:4) = temp1;
% 
%     ang1234 = [q(7).ang1234;q(2).ang1234];    
%     temp1 = ang1234(3,3:4);
%     temp2 = ang1234(14,3:4);
%     ang1234(3,3:4) = temp2; ang1234(14,3:4) = temp1;
%     
% %         temp1 = ang1234(6,3:4);
% %     temp2 = ang1234(8,3:4);
% %     ang1234(6,3:4) = temp2; ang1234(8,3:4) = temp1;
% % 
% %             temp1 = ang1234(26,3:4);
% %     temp2 = ang1234(28,3:4);
% %     ang1234(26,3:4) = temp2; ang1234(28,3:4) = temp1;
% % 
% 
% x_res = 1600; y_res = 900;
%     figure(1);clf
%     set(gcf,'position',[0,0,1600,600])
%     
% for ii = 1:32
%     subplot(4,8,ii)
% 
% temp1 = ang1234(ii,:);
% temp_x = x_res/2+300.*cosd(temp1);
% temp_y = y_res/2+300.*sind(temp1);
% scatter(temp_x,temp_y,15)
% xlim([0 x_res])
% ylim([0 y_res])
% hold on
% scatter(x_res/2,y_res/2,15,'+')
% %line([x_res/2 x_res/2],[0 y_res])
% %line([0 x_res],[y_res/2,y_res/2]);
% title(num2str(ii))
% end
% %saveas(figure(1),'final_set_2.jpg')
% 
% %%
% %%
% indd = find(ind_new(:,1)==370);
% inddd = find(ismember(ind_new(:,1),ind_new(indd,2))==1);
% d_mean_srt_ddd = d_mean_srt(inddd,1);
% %%
% % ang_set1 = ang1234;
% % save('set1.mat','ang_set1')
% 
% ang_set2 = ang1234;
% save('set2.mat','ang_set2')
% 
% 
% %% 
