function Out = Testing_Resampling_Proedures()

% Code to verify the well-foundness of different bootrapping procedures.
% This code was written to justify the soundness of the simulations in
% Lorca-Puls et al (2018, Neuropsychologia).
%
% Note, there is a document called "Verifying bootstrap resampling
% procedure used in Lorca Puls" that writes up the results of these
% simulations. This document should end up with the folder focussed on
% the Lorca-Puls work, although it may also be on a current laptop
% folder called, small samples.
%
% The key property we are seeking to verify is the following: 1) for a given
% "base-sample" (for Diego a large set of patients selected, according to
% some criteria from Cathy's database; for us here, a large set of numbers
% sampled from a Gaussian), sampling directly with replacement
% (direct_wi_rep) generates the same uncertainty (ie. sdts) as sampling
% indirectly with replacement (indirect_wi_rep).
%
% Figure ZZZZ is the main output. The standard deviations being estimated
% in that figure, are of sets of parameter estimates (Means or Cohen's ds), each
% of which was generated from a resampling. So, this mirrors what happens in
% Diego’s analysis, where bootstrap samples are generated, effect sizes of
% each are calculated and the results are put into a distribution. The Std’s
% shown in figure ZZZZ are standard deviations of such sets. So, we repeat
% Diego’s procedure for different bootstrap procedures.
%    Critically, in figure ZZZZ, the direct_wi_rep and indirect_wi_rep lines track each other.
% direct_wi_rep is what we did in Diego’s paper, i.e. bootstrap smaller and smaller
% samples from single (eventually much bigger) data sets. Consequently, the
% smaller samples will have less items in common between samples,
% than the larger samples, and one might think this is what causes the
% increased variability we see as samples get smaller.
%   The indirect_wi_rep line counters this, since it first samples an
% “intermediate” sample from the large, base, sample, which is the size of
% the bootstrap samples then generated from it, e.g. the indirect_wi_rep 60
% data point is generated from bootstrap samples of size 60 from a sample
% that was of size 60. Thus, at each sample size, bootstrap samples were
% generated from samples of the same size as the bootstrap samples.
%     Importantly, this does not change the pattern we observe of sdts.
% My explanation of why the indirect_wi_rep curve is not different from the
% direct_wi_rep one is the following. Even if, on the direct_wi_rep curve,
% you mapped subsamples to sets and then took the intersection, you would
% find less overlap, as subsamples get smaller, that does not automatically
% mean the procedure is broken. That is, the variability in bootstrapping
% is generated from the varying number of times (including zero times) that
% the same items appear in a set. The combinatorics of the variability
% generated in this way is so great that it swamps any decreased overlap
% that would be apparent when going to underlying sets (i.e. removing repetitions).

base_samp_size = 360; % NORMALLY 360; Size of base sample
no_iterations = 10000; % NORMALLY 10000; No of samples taken of a base sample
store_Stat_base = []; % Accumulates the statistic (means or Cohen's ds) of
                      % base samples, to plot at end

No_subsamp_sizes = 11; % This is the number of subsamples sizes tested. Keep
%  at this value, otherwise the code may break.
No_full_sims = 16; % Number of times we run what corresponds to Diego's full
% analysis. This enables us to obtain a more accurate estimate of the means
% and stds of the equivalant of Diego's distributions.

store_All = zeros(5,2,No_subsamp_sizes,No_full_sims);
                % Saves Means and Standard Devs across many replications
% of Diego's procedure - rows are the 5 resampling methods (including
% reference method); columns (2nd dim) are Mean
% and STD; 3rd dim is sizes of subsamples, from big to small; 4th dim is
% replications of Diego's procedure.

% Set statistic to run resampling over. If true statistic is the mean.
statistic_is_mean = false;  % If false, statistic is (Cohen's d) effect size.

% Outer loop that runs the equivalent of Diego's full analysis many times,
% in order that we can obtain estimates of the stds of his
% distributions across many replications.
for r=1:No_full_sims
    
% Generate base sample, which is fixed for each subsample size
base_sample = randn(base_samp_size,1);

% Statistics of base sample
Mean_base = mean(base_sample);
std_base = std(base_sample);
Cohens_d = Mean_base/std_base;

for j=1:No_subsamp_sizes % Cycles through sizes of subsamples

disp('repetitions of Diegos full procedure reached');
r
                     
disp('repetitions of subsample sizes reached');
j

% Step down from biggest to smallest sub sample sizes
sub_samp_size = base_samp_size - (j*30);

% Determine how many splits for disjoint sampling
k_for_disjoint_split = fix(base_samp_size/sub_samp_size);
i_disj_split = 1; % Initialize index for stepping through folds in disjoint sampling.

% Reset data structs to accummulate relevant statistic (either means or
% Cohen's ds) for different resampling/sampling procedures
store_Stat_SubS_direct_wi_rep = []; % Accumulates stats, when sampling with replacement
store_Stat_SubS_direct_wo_rep = []; % Same, but sampling without replacement
store_Stat_SubS_indirect_wi_rep = []; % Same, but sampling with replacement from
                       % a single sample of sub_samp_size from the base_sample
store_Stat_SubS_disjoint = []; % Disjoint sampling, i.e. into
                               % k_for_disjoint_split disjoint sets
store_Stat_Ref_samp = []; % A reference sampling directly from gaussian

% We take a single subsample, which indirect sampling will work from
Subs_for_indirect = datasample(base_sample,sub_samp_size,'Replace',false);

for i=1:no_iterations % Cycle through subsamples of a single base sample
                      % and sub sample size

% disp('repetitions of subsample reached');
% i

% Subsample directly, WITH replacement
SubS_direct_wi_rep = datasample(base_sample,sub_samp_size);
% Store relevant statistic (either mean or cohen's d)
if statistic_is_mean
    Stat_SubS_direct_wi_rep = mean(SubS_direct_wi_rep);
else % Cohen's d
    Stat_SubS_direct_wi_rep = mean(SubS_direct_wi_rep)/std(SubS_direct_wi_rep);
end
store_Stat_SubS_direct_wi_rep = cat(1,store_Stat_SubS_direct_wi_rep,[Stat_SubS_direct_wi_rep]);

% Subsample directly, WITHOUT replacement
SubS_direct_wo_rep = datasample(base_sample,sub_samp_size,'Replace',false);
% Store relevant statistic (either mean or cohen's d)
if statistic_is_mean
    Stat_SubS_direct_wo_rep = mean(SubS_direct_wo_rep);
else % Cohen's d
    Stat_SubS_direct_wo_rep = mean(SubS_direct_wo_rep)/std(SubS_direct_wo_rep);
end
store_Stat_SubS_direct_wo_rep = cat(1,store_Stat_SubS_direct_wo_rep,[Stat_SubS_direct_wo_rep]);

% Indirect subsampling, i.e. from what is itself a smaller sample, WITH replacement
SubS_indirect_wi_rep = datasample(Subs_for_indirect,sub_samp_size);
% Store relevant statistic (either mean or cohen's d)
if statistic_is_mean
    Stat_SubS_indirect_wi_rep = mean(SubS_indirect_wi_rep);
else % Cohen's d
    Stat_SubS_indirect_wi_rep = mean(SubS_indirect_wi_rep)/std(SubS_indirect_wi_rep);
end
store_Stat_SubS_indirect_wi_rep = cat(1,store_Stat_SubS_indirect_wi_rep,[Stat_SubS_indirect_wi_rep]);

% Subsample with disjoint splitting of base sample - cannibalising cross validation folding
if k_for_disjoint_split == 1 %% Sub_sample_size is too big to do disjoint splitting
    Stat_SubS_disjoint = NaN;
else  % Sub_sample size small enough that can start splitting
    if i_disj_split == 1  % New splitting when index at 1
        c = cvpartition(base_samp_size,'KFold',k_for_disjoint_split);
    end
    Inds_subs_disjoint = test(c,i_disj_split); % get indices of test set/ split
    Subs_disjoint = base_sample(Inds_subs_disjoint); % get corresponding subset of base_sample
    % Store relevant statistic (either mean or cohen's d)
    if statistic_is_mean
        Stat_SubS_disjoint = mean(Subs_disjoint);
    else % Cohen's d
        Stat_SubS_disjoint = mean(Subs_disjoint)/std(Subs_disjoint);
    end
    if i_disj_split == k_for_disjoint_split % have completed a full cycle
        i_disj_split = 1;                   % of disjoint splits
    else
        i_disj_split = i_disj_split + 1;
    end
end
store_Stat_SubS_disjoint = cat(1,store_Stat_SubS_disjoint,[Stat_SubS_disjoint]);

% Reference estimate of ground truth
Reference_sample = randn(sub_samp_size,1);
% Store relevant statistic (either mean or cohen's d)
if statistic_is_mean
    Stat_Ref_samp = mean(Reference_sample);
else % Cohen's d
    Stat_Ref_samp = mean(Reference_sample)/std(Reference_sample);
end
store_Stat_Ref_samp = cat(1,store_Stat_Ref_samp,[Stat_Ref_samp]);

end

%%%%%%%%%%%%%%%% This commented fragment is not updated to use of two
% statistics - mean and effect size.
% plot histogram of means of base samples
% figure
% hist(store_Mean_base,30);
% xlabel('Mean');
% ylabel('frequency');
% title('Histogram of means of base samples','fontsize', 12);

% Calculate mean across Stats of each subsample, for a particular base
% sample size
Mean_wi_direct = mean(store_Stat_SubS_direct_wi_rep);
Mean_wo_direct = mean(store_Stat_SubS_direct_wo_rep);
Mean_wi_indirect = mean(store_Stat_SubS_indirect_wi_rep);
Mean_disjoint = mean(store_Stat_SubS_disjoint);
Mean_Reference = mean(store_Stat_Ref_samp);

temp_vect_means = [ Mean_wi_direct Mean_wo_direct Mean_wi_indirect Mean_disjoint Mean_Reference ];

% Calculate std across means of each subsample, for a particular base
% sample size
STD_wi_direct = std(store_Stat_SubS_direct_wi_rep);
STD_wo_direct = std(store_Stat_SubS_direct_wo_rep);
STD_wi_indirect = std(store_Stat_SubS_indirect_wi_rep);
STD_disjoint = std(store_Stat_SubS_disjoint);
STD_Reference = std(store_Stat_Ref_samp);

temp_vect_stds = [ STD_wi_direct STD_wo_direct STD_wi_indirect STD_disjoint STD_Reference ];

temp_vects = cat(2,temp_vect_means',temp_vect_stds');

% store everything for this base-sample size (i.e. j)
store_All(:,:,j,r) = temp_vects;

% if j > 8    %%%%% DEBUGGING - NOT UPDATED TO HAVING TWO STATISTICS - MEAN
% & EFFECT SIZE
%
% plot histogram of means of subsamples of base sample taken with
% replacement
% figure
% hist(store_Mean_SubS_direct_wi_rep,30);
% xlabel('Mean');
% ylabel('frequency');
% str = ['Hist of means of (size ',num2str(sub_samp_size),') subsamps taken directly WITH repl'];
% title(str, 'fontsize', 11);
%
% plot histogram of means of subsamples of base sample taken without
% replacement
% figure
% hist(store_Mean_SubS_direct_wo_rep,30);
% xlabel('Mean');
% ylabel('frequency');
% str = ['Hist of means of (size ',num2str(sub_samp_size),') subsamps taken directly WITHOUT repl'];
% title(str, 'fontsize', 11);
%
% plot histogram of means of (indirect) subsamples, i.e. of an already taken subsample, with
% replacement
% figure
% hist(store_Mean_SubS_indirect_wi_rep,30);
% xlabel('Mean');
% ylabel('frequency');
% str = ['Hist of means of (size ',num2str(sub_samp_size),') subsamps taken indirectly WITH repl'];
% title(str, 'fontsize', 11);
%
% plot histogram of means of disjoint subsamples
% figure
% hist(store_Mean_SubS_disjoint,30);
% xlabel('Mean');
% ylabel('frequency');
% str = ['Hist of means of (size ',num2str(sub_samp_size),') subsamps through disjoint splitting'];
% title(str, 'fontsize', 11);
%
% plot histogram of means of referebce sampling
% figure
% hist(store_Mean_Ref_samp,30);
% xlabel('Mean');
% ylabel('frequency');
% str = ['Hist of means of (size ',num2str(sub_samp_size),') samples from unit guassian'];
% title(str, 'fontsize', 11);
%
% end   % DEBUGGING

end

% The following commented out code would show results for a single
% replication of Diego's procedure.

% x = [ 1 : No_subsamp_sizes ];

% plot means for current repetition of Diego's process
% y1 = squeeze(store_All(1,1,:,r));
% y2 = squeeze(store_All(2,1,:,r));
% y3 = squeeze(store_All(3,1,:,r));
% y4 = squeeze(store_All(4,1,:,r));
% y5 = squeeze(store_All(5,1,:,r));

%figure
% plot(x,y1);
% hold on;
% plot(x,y2);
% plot(x,y3);
% plot(x,y4);
% plot(x,y5);
% legend('direct wi rep','direct w/o rep','indirect wi rep', 'disjoint', 'reference');
% ax = gca;
% ax.XTick = 1:11;
% ax.XTickLabel = {'330','300','270','240','210','180','150','120','90','60','30'};
% xlabel('Subsample/sample sizes (big to small)');
% ylabel('Mean');
% title('Estimate of mean by (sub)sample size & (re)sampling method','fontsize', 11);

% 360    30*j
% hold off;

% plot standard deviations for current repetition of Diego's process
% y1 = squeeze(store_All(1,2,:,r));
% y2 = squeeze(store_All(2,2,:,r));
% y3 = squeeze(store_All(3,2,:,r));
% y4 = squeeze(store_All(4,2,:,r));
% y5 = squeeze(store_All(5,2,:,r));

% figure    %%%%%  
% plot(x,y1);
% hold on;
% plot(x,y2);
% plot(x,y3);
% plot(x,y4);
% plot(x,y5);
% legend('direct wi rep','direct w/o rep','indirect wi rep','disjoint','reference','location','northwest');
% ax = gca;
% ax.XTick = 1:11;
% ax.XTickLabel = {'330','300','270','240','210','180','150','120','90','60','30'};
% xlabel('(Sub)sample sizes (big to small)');
% ylabel('Std');
% title('Estimate of standard dev. by (sub)sample size & (re)sampling method','fontsize', 11);

% Trash: Store everything for this replication of Diego's procedure
% store_All = cat(4,store_All,store_this_rep)

end

% Take mean across the last dimension of store_All, which will give us a
% mean estimate across the replications of Diego's procedure that we have run
store_All_means = mean(store_All,4);
% Same but for std across replications
store_All_std = std(store_All,0,4);

x = [ 1 : No_subsamp_sizes ];

% plot estimation of means across all repetitions of Diego's process
y1 = squeeze(store_All_means(1,1,:));
y2 = squeeze(store_All_means(2,1,:));
y3 = squeeze(store_All_means(3,1,:));
y4 = squeeze(store_All_means(4,1,:));
y5 = squeeze(store_All_means(5,1,:));

figure
plot(x,y1);
hold on;
plot(x,y2);
plot(x,y3);
plot(x,y4);
plot(x,y5);
legend('direct wi rep','direct w/o rep','indirect wi rep', 'disjoint', 'reference');
ax = gca;
ax.XTick = 1:11;
ax.XTickLabel = {'330','300','270','240','210','180','150','120','90','60','30'};
xlabel('Subsample/sample sizes (big to small)');
ylabel('Mean');
if statistic_is_mean
    title('Estimate of mean by (sub)sample size','\bf & (re)sampling method, with statistic as Mean','fontsize', 11);
else
    title('Estimate of mean by (sub)sample size','\bf & (re)sampling method, with statistic as Cohen''s d','fontsize', 11);
end
hold off;

% plot mean estimate of standard deviations across all repetitions of Diego's process
y1 = squeeze(store_All_means(1,2,:));
y2 = squeeze(store_All_means(2,2,:));
y3 = squeeze(store_All_means(3,2,:));
y4 = squeeze(store_All_means(4,2,:));
y5 = squeeze(store_All_means(5,2,:));

figure    %%%%%   figure ZZZZ
plot(x,y1);
hold on;
plot(x,y2);
plot(x,y3);
plot(x,y4);
plot(x,y5);
legend('direct wi rep','direct w/o rep','indirect wi rep','disjoint','reference','location','northwest');
ax = gca;
ax.XTick = 1:11;
ax.XTickLabel = {'330','300','270','240','210','180','150','120','90','60','30'};
xlabel('(Sub)sample sizes (big to small)');
ylabel('Std');
if statistic_is_mean
    title('Estimate of standard dev. by (sub)sample size','\bf & (re)sampling method, with statistic as Mean','fontsize', 11);
else
    title('Estimate of standard dev. by (sub)sample size','\bf & (re)sampling method, with statistic as Cohen''s d','fontsize', 11);
end

hold off;

% plot estimate of std of standard deviations across all repetitions of Diego's process
y1 = squeeze(store_All_std(1,2,:));
y2 = squeeze(store_All_std(2,2,:));
y3 = squeeze(store_All_std(3,2,:));
y4 = squeeze(store_All_std(4,2,:));
y5 = squeeze(store_All_std(5,2,:));

figure
plot(x,y1);
hold on;
plot(x,y2);
plot(x,y3);
plot(x,y4);
plot(x,y5);
legend('direct wi rep','direct w/o rep','indirect wi rep','disjoint','reference','location','northwest');
ax = gca;
ax.XTick = 1:11;
ax.XTickLabel = {'330','300','270','240','210','180','150','120','90','60','30'};
xlabel('(Sub)sample sizes (big to small)');
ylabel('Std (of stds)');
if statistic_is_mean
    title('Estimate of std of stds by (sub)sample size','\bf & (re)sampling method, with statistic as Mean','fontsize', 11);
else
    title('Estimate of std of stds by (sub)sample size','\bf & (re)sampling method, with statistic as Cohen''s d','fontsize', 11);
end


disp("at end");

end
