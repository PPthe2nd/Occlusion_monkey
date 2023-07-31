%%%%% This is the main script of the Occlusion Experiment in monkey area V1
% % % 2021 P Papale fecit w/ code and input from AT Morgan, MW Self & F Wang

clear all
addpath('code_utils')
rng('default');
load('PuBuGn_cmap.mat');

% % % Constants
makedata = 0;
nperms = 1000;
nCond = 24;
nIter = 5;
n100 = 70;
nMaxTrials = 104;

% % % Load and reshape data from the two monkeys
if makedata
    tic
    % % % Select good channels
    monkeyB = 1;
    for monkey = 1:2
        if monkey == 2
            monkeyB = 0;
        end
        if monkeyB
            mydir = 'monkeyB\';
            ftl = [mydir,'MuckliSetBNorm'];
            load(ftl,'SNR','madchan');
            ftl = [mydir,'BarMap_monkeyB_20180104'];
            load(ftl,'RF');
            mySNR = SNR;
            v1rf = [zeros(1,96),ones(1,96)];
            chnorder = 1:192;
            goodRF = ones(1,192);
            good_snr = mySNR > 1;
            goodRF(RF.centrex== 0 | RF.centrey == 0) = 0;
            rfcutoff = 2;
            closech = abs(RF.centrex)<rfcutoff|abs(RF.centrey)<rfcutoff;
            goodRF_v1_B = ~closech&goodRF&v1rf&madchan<0.05 & good_snr;
        else
            mydir = 'monkeyL\';
            ftl = [mydir,'MuckliSetLNorm'];
            load(ftl,'SNR','chnorder','good_chns');
            ftl = [mydir,'RFs_L'];
            load(ftl);
            mySNR = SNR';
            v1rf = area_map == 1;
            goodRF = ismember(chnorder,good_chns);
            good_snr = mySNR > 1;
            good_snr = ismember(chnorder,good_chns(good_snr));
            goodRF_v1_L = goodRF & good_snr;
            chnorder = 1:sum(goodRF_v1_L);
        end
    end
    
    % % % load and normalize MUA responses for all good channels
    signal_smooth = nan([770,nMaxTrials*3*(2*nCond),sum(goodRF_v1_B)+sum(goodRF_v1_L)]);
    monkeyB = 1;
    for monkey = 1:2
        clear MUAWIN muadetails maxsignals
        if monkey == 2
            monkeyB = 0;
        end
        if monkeyB
            mydir = '\\10.41.52.130\vs03-vandc-2\Muckli_reboot\monkeyB\';
            ftl = [mydir,'MuckliSetBNorm'];
            load(ftl,'MUAWIN','muadetails','maxsignals');
            start_chn = 0;
            goodRF_v1 = goodRF_v1_B;
        else
            mydir = '\\10.41.52.130\vs03-vandc-2\Muckli_reboot\monkeyL\';
            ftl = [mydir,'MuckliSetLNorm'];
            load(ftl,'tb','goodT','MUAWIN','muadetails','maxsignals');
            tb = tb(goodT);
            start_chn = length(chns);
            goodRF_v1 = goodRF_v1_L;
        end
        chns = find(goodRF_v1);
        for c = 1:length(chns)
            chn = start_chn+c;
            z = 1;
            for pic = 1:(2*nCond)
                for B = 1:3
                    f = find(muadetails(:,1) == chns(c) & muadetails(:,2) == pic & muadetails(:,3) == B);
                    temp = [];
                    for ff = 1:length(f)
                        temp = cat(1,temp,(MUAWIN(f(ff)).vals)./maxsignals(chns(c)));
                    end
                    for s =1:nMaxTrials
                        if s <= size(temp,1)
                            clear sign_temp
                            if monkeyB
                                sign_temp = resample(temp(s,:),770,839);
                            else
                                sign_temp = temp(s,:);
                            end
                            signal_smooth(:,z,chn) = smooth(sign_temp',18);
                        end
                        if chn == 1
                            tot_stim(z) = pic;
                            tot_B(z) = B;
                        end
                        z=z+1;
                    end
                end
            end
            chn
        end
    end
    toc
    save('signal_smooth','signal_smooth','tot_stim','tot_B','tb','-v7.3');
    clear MUAWIN muadetails maxsignals
else
    load('signal_smooth');
end

% % % plot the average response to each stimulus and condition
for a = 1:(2*nCond)
    mean_class_activity(a,:) = nanmean(nanmean(signal_smooth(:,tot_stim==a,:),2),3);
end
gt = tb>-.1 & tb <.4;
tb_short = tb(gt);

figure;
plot(mean_class_activity(1:nCond,gt)','Color',[0 0.5 1])
hold on
plot(nanmean(mean_class_activity(1:nCond,gt))','Color','b','LineWidth',3)
plot(mean_class_activity(nCond+1:end,gt)','Color',[1 0.5 0])
plot(nanmean(mean_class_activity(nCond+1:end,gt))','Color','r','LineWidth',3)
line([n100 n100],[-.5 2.5],'color','black','LineStyle','--','LineWidth',2)
set(gca,'XTick',1:n100:size(tb_short,2))
set(gca,'XTickLabel',1000*round(tb_short((1:n100:end)),2))
ylabel('Normalized activity')
xlabel({'Time','(ms)'})

% % % compute the latencies of the (positive) differences between Occluded stimuli
occ = mean_class_activity(nCond+1:end,:)/max(nanmean(mean_class_activity,2));
z = 1;
for i = 1:nCond
    msg = ['iteration (of ' num2str(nCond) '): '];
    displayProgress(msg, i, nCond);
    for j = 1:nCond
        if i-j>0
            clear lat rs
            if max(occ(i,gt)-occ(j,gt))>abs(min(occ(i,gt)-occ(j,gt)))
                [lat,coeff,rs] = latencyfit4AM(smooth(occ(i,gt)-occ(j,gt),1),tb(gt),1,0);
                lats_occ(z) = lat;
                fits_occ(z) = rs;
                avg_occ(z,:) = smooth(occ(i,gt)-occ(j,gt),20);
            else
                lats_occ(z) = nan;
                fits_occ(z) = nan;
                avg_occ(z,:) = nan(size(smooth(occ(i,gt)-occ(j,gt),1)));
            end
            z=z+1;
        end
    end
end
lat_occ = latencyfit4AM(nanmean(avg_occ)'-nanmean(nanmean(avg_occ(:,1:n100)))',tb(gt),1,1);

% % % compute the latencies of the (positive) differences between NON-Occluded stimuli
nonocc = mean_class_activity(1:nCond,:)/max(nanmean(mean_class_activity,2));
z = 1;
for i = 1:nCond
    msg = ['iteration (of ' num2str(nCond) '): '];
    displayProgress(msg, i, nCond);
    for j = 1:nCond
        if i-j>0
            clear lat rs
            if max(nonocc(i,gt)-nonocc(j,gt))>abs(min(nonocc(i,gt)-nonocc(j,gt)))
                [lat,coeff,rs] = latencyfit4AM(smooth(nonocc(i,gt)-nonocc(j,gt),20),tb(gt),1,0);
                lats_nonocc(z) = lat;
                fits_nonocc(z) = rs;
                avg_nonocc(z,:) = smooth(nonocc(i,gt)-nonocc(j,gt),20);
            else
                lats_nonocc(z) = nan;
                fits_nonocc(z) = nan;
                avg_nonocc(z,:) = nan(size(smooth(nonocc(i,gt)-nonocc(j,gt),20)));
            end
            z=z+1;
        end
    end
end
lat_nonocc = latencyfit4AM(nanmean(avg_nonocc)'-nanmean(nanmean(avg_nonocc(:,1:n100)))',tb(gt),1,1);

% % % plot latencies
figure
th_ = fits_occ>.5;
lats_all_vect = [lats_nonocc(th_) lats_occ(th_)];
mask = [ones([1 length(lats_nonocc(th_))]) ones([1 length(lats_occ(th_))])*2];
boxplot(lats_all_vect,mask,'plotstyle','compact','symbol','','whisker',0.7193,'colorgroup',mask,'color',[[0 0.5 1];[1 0.5 0]])
labels = [{'NON-occluded'},{'Occluded'}];
ylim([0 .15])
ylabel([{'Time'},{'(ms)'}])
set(gca,'XTick',1:2)
set(gca,'XTickLabel',labels)
set(gca,'YTick',[0,0.05,0.1,0.15])
set(gca,'YTickLabel',[0,0.05,0.1,0.15]*1000)

% % % compute the Wilcoxon signed rank test between the two latency distros
pval_latency_diff = signrank(lats_occ(th_),lats_nonocc(th_));
save(['temp-',date],'*','-v7.3')

% % % Decoding: occluded
nTime = size(signal_smooth,1);

occ_stimMsk = tot_stim > nCond;
fbgood = tot_B>1;
tmpStim = tot_stim(occ_stimMsk & fbgood);
tmpSignal = signal_smooth(:,occ_stimMsk & fbgood,:);
des = dummyvar(tmpStim - nCond);

accuracy_occ = nan([nTime nIter]);
for pp = 1:nIter
    clear trInds
    msg = ['iteration (of ' num2str(nIter) '): '];
    displayProgress(msg, pp, nIter);

    trInds = false(length(des),1);
    trInds(randsample(1:length(des),round(4*(length(des)/5)))) = true;
    for cc = 1:nTime
        clear time_tmpSignal y yA yB xA xB
        time_tmpSignal = squeeze(tmpSignal(cc,:,:));
        time_tmpSignal(isnan(time_tmpSignal)) = 0;
        y = zscore(squeeze(time_tmpSignal));
        
        yA = y(trInds,:); yB = y(~trInds,:);
        xA = des(trInds,:); xB = des(~trInds,:);

        clear stim_temp LDAModel pred
        [stim_temp, ~] = find(xA' == max(xA'));
        LDAModel = fitcdiscr(yA,stim_temp);
        pred = predict(LDAModel,yB);
        clear stim_temp
        [stim_temp, ~] = find(xB' == max(xB'));
        accuracy_occ(cc,pp) = mean(pred==stim_temp);
    end
end

% % % Permutations Occludeds
tmpSignal = tmpSignal(gt,:,:);
tmpSignal_occ = tmpSignal;
nTime = size(tmpSignal,1);
save(['temp-',date],'*')

accuracy_occ_perm = nan([nTime nperms]);
for pp = 1:nperms
    clear trInds
    msg = ['iteration (of ' num2str(nperms) '): '];
    displayProgress(msg, pp, nperms);
    
    % choose a random split
    trInds = false(length(des),1);
    trInds(randsample(1:length(des),round(4*(length(des)/5)))) = true;
    for cc = 1:nTime
        clear time_tmpSignal y yA yB xA xB
        time_tmpSignal = squeeze(tmpSignal(cc,:,:));
        time_tmpSignal(isnan(time_tmpSignal)) = 0;
        y = zscore(squeeze(time_tmpSignal));
        yA = y(trInds,:); yB = y(~trInds,:);
        xA = des(trInds,:); xB = des(~trInds,:);
        clear stim_temp perm_pi LDAModel pred
        [stim_temp, ~] = find(xA' == max(xA'));
        perm_pi = stim_temp(randperm(length(stim_temp),length(stim_temp)));
        LDAModel = fitcdiscr(yA,perm_pi);
        pred = predict(LDAModel,yB);
        clear stim_temp
        [stim_temp, ~] = find(xB' == max(xB'));
        accuracy_occ_perm(cc,pp) = mean(pred==stim_temp);
    end
end
save(['temp-',date],'*','-v7.3')

% % % Time-generalization: occluded
timegen_occ_all = nan([nIter nTime nTime]);
diag_acc = nanmean(accuracy_occ(gt,:),2);
[stim_temp, ~] = find(des' == max(des'));

for pp = 1:nIter
    timegen_occ = nan([nTime nTime]);
    for cc = 1:nTime
        msg = ['iteration (of ' num2str(nTime) '): '];
        displayProgress(msg, cc, nTime);
        clear time_tmpSignal y yA yB xA xB
        time_tmpSignal = squeeze(tmpSignal(cc,:,:));
        time_tmpSignal(isnan(time_tmpSignal)) = 0;
        y = zscore(squeeze(time_tmpSignal));
        clear LDAModel
        % choose a random split
        trInds = false(length(des),1);
        trInds(randsample(1:length(des),round(4*(length(des)/5)))) = true;
        LDAModel = fitcdiscr(y(trInds,:),stim_temp(trInds));
        temp = {};
        parfor dd = 1:nTime
            %         clear time_tmpSignal_test ytest pred
            if dd == cc
                temp_{dd} = diag_acc(cc,1);
            else
                time_tmpSignal_test = squeeze(tmpSignal(dd,:,:));
                time_tmpSignal_test(isnan(time_tmpSignal_test)) = 0;
                ytest = zscore(squeeze(time_tmpSignal_test));
                pred = predict(LDAModel,ytest(~trInds,:));
                temp_{dd} = mean(pred==stim_temp(~trInds));
            end
        end
        timegen_occ(cc,:) = cell2mat(temp_);
    end
    timegen_occ_all(pp,:,:) = timegen_occ;
end

% % % compute p-vals and plot
pfit = paretotails(accuracy_occ_perm(:),0,0.7);
x = 0.0001:0.0001:1;
ecdf = cdf(pfit,x);
x = round(x,4);
mean_tg_occ = squeeze(mean(timegen_occ_all));
for cc = 1:nTime
    for dd = 1:nTime
        pvals_tg_occ(cc,dd) = 1-ecdf(x==round(mean_tg_occ(cc,dd),4))+1/length(ecdf);
    end
end
corr_ps_tg_occ = fdr_bh(pvals_tg_occ,0.05,'dep');
bws = bwboundaries(corr_ps_tg_occ);

figure
imagesc(squeeze(mean(timegen_occ_all))')
set(gca,'XTick',1:n100:size(tb_short,2))
set(gca,'XTickLabel',1000*round(tb_short((1:n100:end)),2))
set(gca,'YTick',1:n100:size(tb_short,2))
set(gca,'YTickLabel',1000*round(tb_short((1:n100:end)),2))
xlabel({'Training time','(ms)'})
ylabel({'Generalization time','(ms)'})
colormap(PuBuGn)
caxis([1/nCond .25])
for i = 1:length(bws)
    hold on
    plot(bws{i}(:,1),bws{i}(:,2),'color','black','LineWidth',1)
end
line([n100 n100],[0 nTime],'color','black','LineStyle','--','LineWidth',3)
line([0 nTime],[n100 n100],'color','black','LineStyle','--','LineWidth',3)
axis square
drawnow
save(['decoding_occ_all-',date],'accuracy_occ','timegen_occ_all','accuracy_occ_perm')
save(['temp-',date],'*','-v7.3')

% % % Decoding: NON-occluded
nonocc_stimMsk = tot_stim <= nCond;
fbgood = tot_B~=2;
tmpStim = tot_stim(nonocc_stimMsk & fbgood);
tmpSignal = signal_smooth(:,nonocc_stimMsk & fbgood,:);
des = dummyvar(tmpStim);
nTime = size(tmpSignal,1);
accuracy_nonocc = nan([nTime nIter]);

for pp = 1:nIter
    clear trInds
    msg = ['iteration (of ' num2str(nIter) '): '];
    displayProgress(msg, pp, nIter);
    
    % choose a random split
    trInds = false(length(des),1);
    trInds(randsample(1:length(des),round(4*(length(des)/5)))) = true;
    
    for cc = 1:nTime
        clear time_tmpSignal y yA yB xA xB
        time_tmpSignal = squeeze(tmpSignal(cc,:,:));
        time_tmpSignal(isnan(time_tmpSignal)) = 0;
        y = zscore(squeeze(time_tmpSignal));
        
        yA = y(trInds,:); yB = y(~trInds,:);
        xA = des(trInds,:); xB = des(~trInds,:);
        
        clear stim_temp LDAModel pred
        [stim_temp, ~] = find(xA' == max(xA'));
        LDAModel = fitcdiscr(yA,stim_temp);%change box into box (g) if want to checl best boxv
        pred = predict(LDAModel,yB);
        clear stim_temp
        [stim_temp, ~] = find(xB' == max(xB'));
        accuracy_nonocc(cc,pp) = mean(pred==stim_temp);
    end
end
save(['temp-',date],'*')

% % % Permutations NON-occluded
tmpSignal = tmpSignal(gt,:,:);
tmpSignal_nonocc = tmpSignal;
nTime = size(tmpSignal,1);

accuracy_nonocc_perm = nan([nTime nperms]);
for pp = 1:nperms
    clear trInds
    msg = ['iteration (of ' num2str(nperms) '): '];
    displayProgress(msg, pp, nperms);
    
    % choose a random split
    trInds = false(length(des),1);
    trInds(randsample(1:length(des),round(4*(length(des)/5)))) = true;
    for cc = 1:nTime
        clear time_tmpSignal y yA yB xA xB
        time_tmpSignal = squeeze(tmpSignal(cc,:,:));
        time_tmpSignal(isnan(time_tmpSignal)) = 0;
        y = zscore(squeeze(time_tmpSignal));
        yA = y(trInds,:); yB = y(~trInds,:);
        xA = des(trInds,:); xB = des(~trInds,:);
        clear stim_temp perm_pi LDAModel pred
        [stim_temp, ~] = find(xA' == max(xA'));
        perm_pi = stim_temp(randperm(length(stim_temp),length(stim_temp)));
        LDAModel = fitcdiscr(yA,perm_pi);
        pred = predict(LDAModel,yB);
        clear stim_temp
        [stim_temp, ~] = find(xB' == max(xB'));
        accuracy_nonocc_perm(cc,pp) = mean(pred==stim_temp);
    end
end
save(['temp-',date],'*','-v7.3')

% % % Time-generalization: NON-occluded
timegen_nonocc_all = nan([nIter nTime nTime]);
diag_acc = nanmean(accuracy_nonocc(gt,:),2);
[stim_temp, ~] = find(des' == max(des'));
for pp = 1:nIter
    timegen_nonocc = nan([nTime nTime]);
    for cc = 1:nTime
        msg = ['iteration (of ' num2str(nTime) '): '];
        displayProgress(msg, cc, nTime);
        clear time_tmpSignal y yA yB xA xB
        time_tmpSignal = squeeze(tmpSignal(cc,:,:));
        time_tmpSignal(isnan(time_tmpSignal)) = 0;
        y = zscore(squeeze(time_tmpSignal));
        clear LDAModel trInds
        % choose a random split
        trInds = false(length(des),1);
        trInds(randsample(1:length(des),round(4*(length(des)/5)))) = true;
        LDAModel = fitcdiscr(y(trInds,:),stim_temp(trInds));
        temp = {};
        parfor dd = 1:nTime
            %         clear time_tmpSignal_test ytest pred
            if dd == cc
                temp_{dd} = diag_acc(cc,1);
            else
                time_tmpSignal_test = squeeze(tmpSignal(dd,:,:));
                time_tmpSignal_test(isnan(time_tmpSignal_test)) = 0;
                ytest = zscore(squeeze(time_tmpSignal_test));
                pred = predict(LDAModel,ytest(~trInds,:));
                temp_{dd} = mean(pred==stim_temp(~trInds));
            end
        end
        timegen_nonocc(cc,:) = cell2mat(temp_);
    end
    timegen_nonocc_all(pp,:,:) = timegen_nonocc;
end

% % % compute p-vals and plot
pfit = paretotails(accuracy_nonocc_perm(:),0,0.7);
x = 0.0001:0.0001:1;
ecdf = cdf(pfit,x);
x = round(x,4);
mean_tg_nonocc = squeeze(mean(timegen_nonocc_all));
for cc = 1:nTime
    for dd = 1:nTime
        pvals_tg_nonocc(cc,dd) = 1-ecdf(x==round(mean_tg_nonocc(cc,dd),4))+1/length(ecdf);
    end
end
pvals_tg_nonocc(pvals_tg_nonocc>1) = 1;
corr_ps_tg_nonocc = fdr_bh(pvals_tg_nonocc,0.05,'dep');
bws = bwboundaries(corr_ps_tg_nonocc);

figure
imagesc(squeeze(mean(timegen_nonocc_all))')
set(gca,'XTick',1:n100:size(tb_short,2))
set(gca,'XTickLabel',1000*round(tb_short((1:n100:end)),2))
set(gca,'YTick',1:n100:size(tb_short,2))
set(gca,'YTickLabel',1000*round(tb_short((1:n100:end)),2))
xlabel({'Training time','(ms)'})
ylabel({'Generalization time','(ms)'})
colormap(PuBuGn)
caxis([1/nCond 1])
for i = 1:length(bws)
    hold on
    plot(bws{i}(:,1),bws{i}(:,2),'color','black','LineWidth',1)
end
line([n100 n100],[0 nTime],'color','black','LineStyle','--','LineWidth',3)
line([0 nTime],[n100 n100],'color','black','LineStyle','--','LineWidth',3)
axis square
drawnow
save(['decoding_nonocc_all-',date],'accuracy_nonocc','timegen_nonocc_all','accuracy_nonocc_perm')
save(['temp-',date],'*','-v7.3')

% % % Plot accuracies
mean_acc_occ = squeeze(mean(accuracy_occ(gt,:),2));
pfit = paretotails(accuracy_occ_perm(:),0,0.7);
x = 0.0001:0.0001:1;
ecdf = cdf(pfit,x);
x = round(x,4);
for cc = 1:nTime
    pvals_occ(cc) = 1-ecdf(x==round(mean_acc_occ(cc),4))+1/length(ecdf);
end
corr_ps_occ = fdr_bh(pvals_occ,0.05,'dep');

mean_acc_nonocc = squeeze(mean(accuracy_nonocc(gt,:),2));
pfit = paretotails(accuracy_nonocc_perm(:),0,0.7);
x = 0.0001:0.0001:1;
ecdf = cdf(pfit,x);
x = round(x,4);
for cc = 1:nTime
    pvals_nonocc(cc) = 1-ecdf(x==round(mean_acc_nonocc(cc),4))+1/length(ecdf);
end
corr_ps_nonocc = fdr_bh(pvals_nonocc,0.05,'dep');

figure;
plot(max(accuracy_occ(gt,:),[],2),'Color',[1 0.5 0]);hold on;plot(min(accuracy_occ(gt,:),[],2),'Color',[1 0.5 0]);
plot(mean(accuracy_occ(gt,:),2),'Color','r');
line([find(corr_ps_occ,1) find(corr_ps_occ,1,'last')],[.99 .99],'color','r','LineWidth',2)
plot(max(accuracy_nonocc(gt,:),[],2),'Color',[0 0.5 1]);plot(min(accuracy_nonocc(gt,:),[],2),'Color',[0 0.5 1])
plot(mean(accuracy_nonocc(gt,:),2),'Color','b');
line([find(corr_ps_nonocc,1) find(corr_ps_nonocc,1,'last')],[1 1],'color','b','LineWidth',2)
line([n100 n100],[0 nTime],'color','black','LineStyle','--','LineWidth',2)
ylim([0 1]);
set(gca,'XTick',1:n100:size(tb_short,2))
set(gca,'XTickLabel',1000*round(tb_short((1:n100:end)),2))
ylabel('Decoding accuracy')
xlabel({'Time','(ms)'})
drawnow

% % % Time-generalization: cross-decoding
% train on OCC and test on NONOCC
timegen_cross_occ_to_nonocc_all = nan([nIter nTime nTime]);
for pp = 1:nIter
    timegen_cross_occ_to_nonocc = nan([nTime nTime]);
    for cc = 1:nTime
        msg = ['iteration (of ' num2str(nTime) '): '];
        displayProgress(msg, cc, nTime);
        clear time_tmpSignal y yA yB xA xB
        time_tmpSignal = squeeze(tmpSignal_occ(cc,:,:));
        time_tmpSignal(isnan(time_tmpSignal)) = 0;
        y = zscore(squeeze(time_tmpSignal));
        clear LDAModel trInds
        % choose a random split
        trInds = false(length(des),1);
        trInds(randsample(1:length(des),round(4*(length(des)/5)))) = true;
        LDAModel = fitcdiscr(y(trInds,:),stim_temp(trInds));
        temp = {};
        test_labels = stim_temp(~trInds);
        parfor dd = 1:nTime
            time_tmpSignal_test = squeeze(tmpSignal_nonocc(dd,:,:));
            time_tmpSignal_test(isnan(time_tmpSignal_test)) = 0;
            ytest = zscore(squeeze(time_tmpSignal_test));
            pred = predict(LDAModel,ytest(~trInds,:));
            temp_{dd} = mean(pred==test_labels);
        end
        timegen_cross_occ_to_nonocc(cc,:) = cell2mat(temp_);
    end
    timegen_cross_occ_to_nonocc_all(pp,:,:) = timegen_cross_occ_to_nonocc;
end
save(['temp-',date],'*','-v7.3')

% % % Permutation OCC --> NONOCC
accuracy_cross_occ_to_nonocc_perm = nan([nTime nperms]);
for pp = 1:nperms
    clear trInds
    msg = ['iteration (of ' num2str(nperms) '): '];
    displayProgress(msg, pp, nperms);
    
    % choose a random split
    trInds = false(length(des),1);
    trInds(randsample(1:length(des),round(4*(length(des)/5)))) = true;
    for cc = 1:nTime
        clear time_tmpSignal y yA yB xA xB time_tmpSignal_test ytest
        time_tmpSignal = squeeze(tmpSignal_occ(cc,:,:));
        time_tmpSignal(isnan(time_tmpSignal)) = 0;
        y = zscore(squeeze(time_tmpSignal));
        
        time_tmpSignal_test = squeeze(tmpSignal_nonocc(cc,:,:));
        time_tmpSignal_test(isnan(time_tmpSignal_test)) = 0;
        ytest = zscore(squeeze(time_tmpSignal_test));
        
        yA = y(trInds,:); yB = ytest(~trInds,:);
        xA = des(trInds,:); xB = des(~trInds,:);
        clear stim_temp perm_pi LDAModel pred
        [stim_temp, ~] = find(xA' == max(xA'));
        perm_pi = stim_temp(randperm(length(stim_temp),length(stim_temp)));
        LDAModel = fitcdiscr(yA,perm_pi);
        pred = predict(LDAModel,yB);
        clear stim_temp
        [stim_temp, ~] = find(xB' == max(xB'));
        accuracy_cross_occ_to_nonocc_perm(cc,pp) = mean(pred==stim_temp);
    end
end

% % % compute p-vals and plot
pfit = paretotails(accuracy_cross_occ_to_nonocc_perm(:),0,0.7);
x = 0.0001:0.0001:1;
ecdf = cdf(pfit,x);
x = round(x,4);
mean_tg_n_cross_occ_nonocc = squeeze(mean(timegen_cross_occ_to_nonocc_all));
for cc = 1:nTime
    for dd = 1:nTime
        pvals_tg_cross_occ_nonocc(cc,dd) = 1-ecdf(x==round(mean_tg_n_cross_occ_nonocc(cc,dd),4))+1/length(ecdf);
    end
end
corr_ps_tg_cross_occ_nonocc = fdr_bh(pvals_tg_cross_occ_nonocc,0.05,'dep');
bws = bwboundaries(corr_ps_tg_cross_occ_nonocc);

figure
imagesc(squeeze(mean(timegen_cross_occ_to_nonocc_all))')
set(gca,'XTick',1:n100:size(tb_short,2))
set(gca,'XTickLabel',1000*round(tb_short((1:n100:end)),2))
set(gca,'YTick',1:n100:size(tb_short,2))
set(gca,'YTickLabel',1000*round(tb_short((1:n100:end)),2))
xlabel({'Training time','(ms)'})
ylabel({'Generalization time','(ms)'})
colormap(PuBuGn)
caxis([1/nCond .15])
for i = 1:length(bws)
    hold on
    plot(bws{i}(:,1),bws{i}(:,2),'color','black','LineWidth',1)
end
line([n100 n100],[0 nTime],'color','black','LineStyle','--','LineWidth',3)
line([0 nTime],[n100 n100],'color','black','LineStyle','--','LineWidth',3)
axis square
drawnow
save(['crossdecoding_occ_to_nonocc_all-',date],'timegen_cross_occ_to_nonocc_all','accuracy_cross_occ_to_nonocc_perm')
save(['temp-',date],'*','-v7.3')

% train on NONOCC and test on OCC
timegen_cross_nonocc_to_occ_all = nan([nIter nTime nTime]);
[stim_temp, ~] = find(des' == max(des'));
for pp = 1:nIter
    timegen_cross_nonocc_to_occ = nan([nTime nTime]);
    for cc = 1:nTime
        msg = ['iteration (of ' num2str(nTime) '): '];
        displayProgress(msg, cc, nTime);
        clear time_tmpSignal y yA yB xA xB
        time_tmpSignal = squeeze(tmpSignal_nonocc(cc,:,:));
        time_tmpSignal(isnan(time_tmpSignal)) = 0;
        y = zscore(squeeze(time_tmpSignal));
        clear LDAModel trInds
        % choose a random split
        trInds = false(length(des),1);
        trInds(randsample(1:length(des),round(4*(length(des)/5)))) = true;
        LDAModel = fitcdiscr(y(trInds,:),stim_temp(trInds));
        temp = {};
        parfor dd = 1:nTime
            time_tmpSignal_test = squeeze(tmpSignal_occ(dd,:,:));
            time_tmpSignal_test(isnan(time_tmpSignal_test)) = 0;
            ytest = zscore(squeeze(time_tmpSignal_test));
            pred = predict(LDAModel,ytest(~trInds,:));
            temp_{dd} = mean(pred==stim_temp(~trInds));
        end
        timegen_cross_nonocc_to_occ(cc,:) = cell2mat(temp_);
    end
    timegen_cross_nonocc_to_occ_all(pp,:,:) = timegen_cross_nonocc_to_occ;
end
save(['temp-',date],'*','-v7.3')

% % % Permutation NONOCC --> OCC
accuracy_cross_nonocc_to_occ_perm = nan([nTime nperms]);
for pp = 1:nperms
    clear trInds
    msg = ['iteration (of ' num2str(nperms) '): '];
    displayProgress(msg, pp, nperms);
    
    % choose a random split
    trInds = false(length(des),1);
    trInds(randsample(1:length(des),round(4*(length(des)/5)))) = true;
    for cc = 1:nTime
        clear time_tmpSignal y yA yB xA xB time_tmpSignal_test ytest
        time_tmpSignal = squeeze(tmpSignal_nonocc(cc,:,:));
        time_tmpSignal(isnan(time_tmpSignal)) = 0;
        y = zscore(squeeze(time_tmpSignal));
        
        time_tmpSignal_test = squeeze(tmpSignal_occ(cc,:,:));
        time_tmpSignal_test(isnan(time_tmpSignal_test)) = 0;
        ytest = zscore(squeeze(time_tmpSignal_test));
        
        yA = y(trInds,:); yB = ytest(~trInds,:);
        xA = des(trInds,:); xB = des(~trInds,:);
        clear stim_temp perm_pi LDAModel pred
        [stim_temp, ~] = find(xA' == max(xA'));
        perm_pi = stim_temp(randperm(length(stim_temp),length(stim_temp)));
        LDAModel = fitcdiscr(yA,perm_pi);
        pred = predict(LDAModel,yB);
        clear stim_temp
        [stim_temp, ~] = find(xB' == max(xB'));
        accuracy_cross_nonocc_to_occ_perm(cc,pp) = mean(pred==stim_temp);
    end
end
save(['temp-',date],'*','-v7.3')

% % % compute p-vals and plot
pfit = paretotails(accuracy_cross_nonocc_to_occ_perm(:),0,0.7);
x = 0.0001:0.0001:1;
ecdf = cdf(pfit,x);
x = round(x,4);
mean_tg_n_cross_nonocc_occ = squeeze(mean(timegen_cross_nonocc_to_occ_all));
for cc = 1:nTime
    for dd = 1:nTime
        pvals_tg_cross_nonocc_occ(cc,dd) = 1-ecdf(x==round(mean_tg_n_cross_nonocc_occ(cc,dd),4))+1/length(ecdf);
    end
end
corr_ps_tg_cross_nonocc_occ = fdr_bh(pvals_tg_cross_nonocc_occ,0.05,'dep');
bws = bwboundaries(corr_ps_tg_cross_nonocc_occ);

figure
imagesc(squeeze(mean(timegen_cross_nonocc_to_occ_all))')
set(gca,'XTick',1:n100:size(tb_short,2))
set(gca,'XTickLabel',1000*round(tb_short((1:n100:end)),2))
set(gca,'YTick',1:n100:size(tb_short,2))
set(gca,'YTickLabel',1000*round(tb_short((1:n100:end)),2))
xlabel({'Training time','(ms)'})
ylabel({'Generalization time','(ms)'})
colormap(PuBuGn)
caxis([1/nCond .1])
for i = 1:length(bws)
    hold on
    plot(bws{i}(:,1),bws{i}(:,2),'color','black','LineWidth',1)
end
line([n100 n100],[0 nTime],'color','black','LineStyle','--','LineWidth',3)
line([0 nTime],[n100 n100],'color','black','LineStyle','--','LineWidth',3)
axis square
drawnow
save(['crossdecoding_nonocc_to_occ_all-',date],'timegen_cross_nonocc_to_occ_all','accuracy_cross_nonocc_to_occ_perm')
save(['temp-',date],'*','-v7.3')

% % % Plot cross-decoding accuracies
[xx, yy] = meshgrid(1:size(timegen_cross_nonocc_to_occ,1),1:size(timegen_cross_nonocc_to_occ,1));
for pp = 1:nIter
    clear temp_1 temp_2
    temp_1 = squeeze(timegen_cross_nonocc_to_occ_all(pp,:,:));
    temp_2 = squeeze(timegen_cross_occ_to_nonocc_all(pp,:,:));
    accuracy_cross_nonocc_to_occ(:,pp) = temp_1(logical(xx==yy));
    accuracy_cross_occ_to_nonocc(:,pp) = temp_2(logical(xx==yy));
end

mean_acc_occ_to_nonocc = squeeze(mean(accuracy_cross_occ_to_nonocc,2));
pfit = paretotails(accuracy_cross_occ_to_nonocc_perm(:),0,0.7);
x = 0.0001:0.0001:1;
ecdf = cdf(pfit,x);
x = round(x,4);
for cc = 1:nTime
    pvals_mean_acc_occ_to_nonocc(cc) = 1-ecdf(x==round(mean_acc_occ_to_nonocc(cc),4))+1/length(ecdf);
end
corr_ps_mean_acc_occ_to_nonocc = fdr_bh(pvals_mean_acc_occ_to_nonocc,0.05,'dep');

mean_acc_nonocc_to_occ = squeeze(mean(accuracy_cross_nonocc_to_occ,2));
pfit = paretotails(accuracy_cross_nonocc_to_occ_perm(:),0,0.7);
x = 0.0001:0.0001:1;
ecdf = cdf(pfit,x);
x = round(x,4);
for cc = 1:nTime
    pvals_mean_acc_nonocc_to_occ(cc) = 1-ecdf(x==round(mean_acc_nonocc_to_occ(cc),4))+1/length(ecdf);
end
corr_ps_mean_acc_mean_acc_nonocc_to_occ = fdr_bh(pvals_mean_acc_nonocc_to_occ,0.05,'dep');

figure;
plot(max(accuracy_cross_occ_to_nonocc,[],2),'Color',[1 0.5 0]);hold on;plot(min(accuracy_cross_occ_to_nonocc,[],2),'Color',[1 0.5 0]);
plot(mean(accuracy_cross_occ_to_nonocc,2),'Color','r')
line([find(corr_ps_mean_acc_occ_to_nonocc,1) find(corr_ps_mean_acc_occ_to_nonocc,1,'last')],[.14 .14],'color','r','LineWidth',2)
plot(max(accuracy_cross_nonocc_to_occ,[],2),'Color',[0 0.5 1]);plot(min(accuracy_cross_nonocc_to_occ,[],2),'Color',[0 0.5 1])
plot(mean(accuracy_cross_nonocc_to_occ,2),'Color','b')
line([n100 n100],[-.5 2.5],'color','black','LineStyle','--','LineWidth',2)
if sum(corr_ps_mean_acc_mean_acc_nonocc_to_occ)>0
    line([find(corr_ps_mean_acc_mean_acc_nonocc_to_occ,1) find(corr_ps_mean_acc_mean_acc_nonocc_to_occ,1,'last')],[.15 .15],'color','r','LineWidth',2)
end
ylim([0 .15]);
set(gca,'XTick',1:n100:size(tb_short,2))
set(gca,'XTickLabel',1000*round(tb_short((1:n100:end)),2))
ylabel('Decoding accuracy')
xlabel({'Time','(ms)'})
drawnow

% % % compute RSA between the monkey V1 MUA and human V1 BOLD
rng('default');
nTime = size(signal_smooth,1);
occ_stimMsk = tot_stim > nCond;
fbgood = tot_B>1;
tmpStim = tot_stim(occ_stimMsk & fbgood);
tmpSignal = signal_smooth(:,occ_stimMsk & fbgood,:);
des = dummyvar(tmpStim - nCond);

ldcV_occ = zeros(276,nTime,nIter);
for pp = 1:nIter
    clear trInds
    msg = ['iteration (of ' num2str(nIter) '): '];
    displayProgress(msg, pp, nIter);
    
    % choose a random split
    trInds = false(length(des),1);
    trInds(randsample(1:length(des),round(4*(length(des)/5)))) = true;
    for cc = 1:nTime
        clear time_tmpSignal y yA yB xA xB
        time_tmpSignal = squeeze(tmpSignal(cc,:,:));
        time_tmpSignal(isnan(time_tmpSignal)) = 0;
        y = squeeze(time_tmpSignal);
        
        yA = y(trInds,:); yB = y(~trInds,:);
        xA = des(trInds,:); xB = des(~trInds,:);
        ldcV_occ(:,cc,pp) = abs(fisherDiscrContrast(xA,yA,xB,yB,nCond));
    end
end
ldcMean_occ = mean(ldcV_occ,3);

load('gist3_ldc_V1.mat', 'ldcV')
human_subs = squeeze(nanmean(ldcV(:,:,3,:),2));
human = nanmean(human_subs,2);

for i = 1:size(ldcMean_occ,2)
    ch(i) = corr(ldcMean_occ(:,i),human,'type','Spearman');
end

% % % noise-ceiling
for i = 1:size(human_subs,2)
    upper_bounds(i) = corr(human_subs(:,i),human,'type','Spearman');
end
upper_bound = mean(upper_bounds);
for i = 1:size(human_subs,2)
    lower_bounds(i) = corr(human_subs(:,i),nanmean(human_subs(:,setdiff(1:size(human_subs,2),i)),2),'type','Spearman');
end
lower_bound = nanmean(lower_bounds);

% % Permutations corr
for i = 1:size(ldcMean_occ,2)
    for j = 1:nperms*10
        chp(j,i) = corr(ldcMean_occ(randperm(size(ldcMean_occ,1),size(ldcMean_occ,1)),i),human,'type','Spearman');
    end
end

for i = 1:size(ldcMean_occ,2)
    msg = ['iteration (of ' num2str(size(ldcMean_occ,2)) '): '];
    displayProgress(msg, i, size(ldcMean_occ,2));
    funcia = @(idx) (corr(ldcMean_occ(idx,i),human(idx),'type','Spearman'));
    boot_corrs(i,:) = bootstrp(1000,funcia,1:length(human));
end

% % % plot
ch_gt = ch(gt);
chp_gt = chp(:,gt)';
boot_corrs_gt = boot_corrs(gt,:);
for i = 1:length(ch_gt)
    pvals_corr(i) = sum(chp_gt(i,:)>ch_gt(i))/size(chp_gt,2)+1/size(chp_gt,2);
end
pvals_corr(pvals_corr>1) = 1;
corr_ps_human = fdr_bh(pvals_corr,0.05,'dep');
bws = bwboundaries(corr_ps_human);

figure
plot(ch_gt,'r')
hold on
plot(ch_gt-std(boot_corrs_gt'),'Color',[1 0.5 0])
plot(ch_gt+std(boot_corrs_gt'),'Color',[1 0.5 0])
ylim([-.25 .5])
set(gca,'XTick',1:n100:size(tb_short,2))
set(gca,'XTickLabel',1000*round(tb_short((1:n100:end)),2))
ylabel({'Correlation','(Spearman''s rho)'})
xlabel({'Time','(ms)'})
line([1 size(tb_short,2)],[lower_bound lower_bound],'color','black')
line([1 size(tb_short,2)],[upper_bound upper_bound],'color','black')
for i = 1:length(bws)
    hold on
    line([min(bws{i}(:,2)) max(bws{i}(:,2))],[.48 .48],'color','black','LineWidth',2)
end
drawnow

% % % Compute and print the test-results to be reported
pval_decoding = signrank(mean(accuracy_occ,2),mean(accuracy_nonocc,2));
pval_crossdecoding = signrank(mean(accuracy_cross_occ_to_nonocc,2),mean(accuracy_cross_nonocc_to_occ,2));

disp(['%%%%%%%%% ',' Test results! Pvals are rounded to 4th decimal:'])
disp(['Latencies: occ > nonocc! pval = ',num2str(round(pval_latency_diff,4)),'; Wilcoxon signed-rank test'])
disp(['Dedcoding: nonocc > occ! pval = ',num2str(round(pval_decoding,4)),'; Wilcoxon signed-rank test'])
disp(['Cross-dedcoding: occ_nonocc > nonocc_occ! pval = ',num2str(round(pval_crossdecoding,4)),'; Wilcoxon signed-rank test'])

% % % Save everything
save(['mainV1-',date],'*','-v7.3')

