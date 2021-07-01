
% 
% for i = 1:10
%     
%     fig=figure;
%     subplot(2,1,1)
%     imagesc(df.Pr{i}');
%     axis xy
%     colormap magma
%     hold on;
%     subplot(2,1,2)
%     map = imresize(df.Pr{i}(~df.inactive_bins{i},:),[28,28]);
%     imagesc(map');
%     axis xy
%     colormap magma
%     hold on;
%     
% end

% map = zeros(length(df.Pr),28,28);
% for i = 1:length(df.Pr)
%     map(i,:,:) = imresize(df.Pr{i}(~df.inactive_bins{i},:),[28,28]);
% end

% writeNPY(map, 'F:\Projects\PAE_PlaceCell\analysis\replay\maps.npy')



% mdl = fitlm(x,y,'Weights',w);
% df.r2(i) = mdl.Rsquared.Ordinary;

% function [x,y,score] = clean_trajectory(Pr,inactive_bins)

i = find(df.ripple_number==26422)

map = Pr.Pr{i};
cur_inactive_bins = inactive_bins{i};

[w,y] = max(map,[],2);

x = (1:size(map,1))';

y(inactive_bins{i}) = [];
x(inactive_bins{i}) = [];
w(inactive_bins{i}) = [];

tbl = table(y,x,x.^2,x.^3,'VariableNames',{'y','x','x2','x3'});

lm_obs = fitlm(tbl,'y~x','Weights',w);

obs = sum(lm_obs.Residuals.Raw.^2);

fig=figure;
imagesc(map');
axis xy
colormap magma
hold on;
    
plot(x,y,'.r')
hold on
plot(x,lm_obs.Fitted,'r')


clear temp_map
for i_shuff = 1:1000
    % column cycle shuff
    r = floor(1 + (size(map,2)-1).*rand(size(map,1),1));
    for s = 1:size(map,1)
        temp_map(s,:) = circshift(map(s,:), r(s));
    end
    
    [w,y] = max(temp_map,[],2);

    x = (1:size(temp_map,1))';

    y(cur_inactive_bins) = [];
    x(cur_inactive_bins) = [];
    w(cur_inactive_bins) = [];
    
    tbl = table(y,x,x.^2,x.^3,'VariableNames',{'y','x','x2','x3'});

    lm = fitlm(tbl,'y~x+x2+x3','Weights',w);
    ssr(i,1) = sum(lm.Residuals.Raw.^2);
        
    plot(x,lm.Coefficients.Estimate(1) +...
    x * lm.Coefficients.Estimate(2) +...
    x.^2 * lm.Coefficients.Estimate(3) +...
    x.^3 * lm.Coefficients.Estimate(4),'color',[1,1,1,.1])
end
    
plot(x,lm_obs.Coefficients.Estimate(1) +...
    x * lm_obs.Coefficients.Estimate(2) +...
    x.^2 * lm_obs.Coefficients.Estimate(3) +...
    x.^3 * lm_obs.Coefficients.Estimate(4),'r')

%%
z = get_rZ(obs,ssr)
p = get_pvalue(obs,ssr)
%%
function z = get_rZ(obs,null)
z = (abs(obs) - abs(mean(null))) / abs(std(null));
end

function p = get_pvalue(obs,null)
p = (sum(abs(null) >= abs(obs)) + 1) / (length(null) + 1);
end


%%
df = readtable('F:\Projects\PAE_PlaceCell\replay\processed\replay_df.csv');
Pr = load('F:\Projects\PAE_PlaceCell\replay\processed\replay_Pr.mat','Pr');
load('F:\Projects\PAE_PlaceCell\replay\processed\replay_inactive_bins.mat','inactive_bins');

df.Pr = Pr.Pr;
df.inactive_bins = inactive_bins;
alpha = 0.2
for i = 1:length(df.session)
    
    fig=figure;
    imagesc(df.Pr{i}');
    axis xy
    colormap magma
    hold on;

    % get max per time step
    [w,y] = max(df.Pr{i},[],2);
    % each spatial bin is 3cm
    y = y * 3; 
    % time bins
    x = 1:size(df.Pr{i},1);
    % mdl = fitlm(x,y,'Weights',w);
    % df.r2(i) = mdl.Rsquared.Ordinary;
    % remove bins without activity
    y(df.inactive_bins{i}) = [];
    x(df.inactive_bins{i}) = [];
    % median smooth to clean up trajectory
    y = smoothdata(y,'movmedian',3);
    
    plot(x,y/3,'w','linewidth',3)
    
    % get spatial difference between bins
    dy = abs(diff(y));
    % get cumulative distance 
    df.traj_dist(i) = sum(dy);
    % calculate avg speed of trajectory (dist(cm) / time(sec))
    df.traj_speed(i) = df.traj_dist(i) / (((max(x) - min(x)) * 10) / 1000);
    % get mean step size 
    df.traj_step(i) = mean(dy);
    
    % save pngs to folder (for final version, save to svg)
    sig = df.pvalue_cellID_shuf(i) <= alpha & df.pvalue_circular_shuf(i) <= alpha & df.pvalue_column_cycle(i) <= alpha;
    title(['sig ',num2str(sig)])
    
    id = [df.session{i},'_',...
        df.ep_type{i},'_',...
        num2str(df.Var1(i))];
    saveas(fig,['F:\Projects\PAE_PlaceCell\analysis\replay\trajectory_test_figs\',id,'.png'])
    pause(0.00001)
    close all
end

alpha = 0.05
idx = df.pvalue_cellID_shuf <= alpha & df.pvalue_circular_shuf <= alpha & df.pvalue_column_cycle <= alpha;

figure;histogram(df.traj_step(idx));
%%

%%
i = find(df.ripple_number==32855)
i = find(df.ripple_number==11)
i=2873

figure;
imagesc(Pr.Pr{i}');
axis xy
colormap magma
hold on;
[~,y] = max(Pr.Pr{i},[],2);

original_y = y;

x = 1:size(Pr.Pr{i},1);

plot(x,y,'linewidth',3)

y(inactive_bins{i}) = [];
x(inactive_bins{i}) = [];

plot(x,y,'linewidth',3)

% bad_idx = find(diff(y*3) > 50)' + 1;
% y(bad_idx) = [];
% x(bad_idx) = [];
% 
% plot(x,y,'linewidth',3)
plot(x,smoothdata(y,'movmedian',3),'linewidth',3)

final_y = smoothdata(y,'movmedian',3);

% plot(rescale(data.frames(idx,1),1,size(Pr.Pr{i},1)),data.frames(idx,2),'r')
title([df.ep_type(i), num2str(df.bayesLinearWeighted(i)), df.replay_type{i}])


legend({'original','inactive bins removed','median smooth'},'Location','best')

corr(final_y,original_y(ismember(1:length(original_y),x)))
%%
for i = 1:length(Pr.Pr)
    [x{i},y{i},score(i)] = clean_trajectory(Pr.Pr{i},inactive_bins{i});
end
%%
function [x,y,score] = clean_trajectory(Pr,inactive_bins)
[~,y] = max(Pr,[],2);

original_y = y;

x = 1:size(Pr,1);


y(inactive_bins) = [];
x(inactive_bins) = [];


% bad_idx = find(diff(y*3) > 50)' + 1;
% y(bad_idx) = [];
% x(bad_idx) = [];


y = smoothdata(y,'movmedian',3);
score = corr(y,original_y(ismember(1:length(original_y),x)));


end




%%
