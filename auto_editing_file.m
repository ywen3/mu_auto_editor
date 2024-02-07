function [fileNameNew, editInfo] = auto_editing_file(fileName, para, plotFlagTmp)
% auto_eiditing_file: automatically editing all motor unit spike trains
% input:
%       fileName: the full file path specify the file containing fastICA
%       output: the edited file saved as fileNameNew
%
% Author: Yue Wen at University of Central Florida
%

if nargin == 0
    % get folder and fileName
    [fileName, path] = uigetfile('*.mat');
    fileName = [path, fileName];
    para = [0.331, 0.339, 0.539, 0.508];    % optimized using majority vote data
    plotFlagTmp = 1;
elseif nargin == 1
    para = [0.331, 0.339, 0.539, 0.508];    % optimized using majority vote data
    plotFlagTmp = 0;
elseif nargin == 2
    plotFlagTmp = 0;
end

% construct file name and path
[filePath, fileTag, fileExt] = fileparts(fileName);
rootFolder = pwd;
if ~isempty(filePath)
    cd(filePath);
end
fileName = [fileTag, fileExt];

% enable plot
global test_mu_id;
global plotFlag;
plotFlag = plotFlagTmp;
global muTag;
if 0
    global plotFlag
    plotFlag = 0;
end

% set parameters
frHighThreshold = 50;  % firing rate threshold
maxCycles = 5;

% load data
data = load(fileName);
MUPulses = data.MUPulses;
IPTs = data.IPTs;
SIG = data.SIG;
MUIDs = data.MUIDs;
PNR = data.PNR;
discardChannelsVec = data.discardChannelsVec;
fsamp = data.fsamp;
refSignal = data.ref_signal;
if isfield(data, 'rule_based_editting')
    processFlag = [data.rule_based_editting.flag];
else
    processFlag = ones(size(PNR));
end

% process reference signal (torque or force)
refSignal(isnan(refSignal)) = 0;
refSignalSmoothed = smooth(refSignal, 2048);
refSignal = refSignal - min(refSignalSmoothed);
refSignalSmoothed = refSignalSmoothed - min(refSignalSmoothed);
if max(refSignalSmoothed)
    refSignal = refSignal./max(refSignalSmoothed);
end
signalLength = length(refSignal);

% create folder to store intermediate results
nameItem = split(fileName, '.');
folder = [nameItem{1}, nameItem{2}];
if ~exist(folder, 'dir')
    mkdir(folder)
end
cd(folder);

%% get the same segments as each operator for online data validation
MUCount = length(MUPulses(:));
segments = zeros(MUCount, 2);
spikeRange = zeros(MUCount, 2);
for mi = 1:MUCount
    if isempty(MUPulses{mi})
        continue;
    end
    spikeRange(mi, 1) = max([1, MUPulses{mi}(1)-2048]);
    spikeRange(mi, 2) = min([signalLength, MUPulses{mi}(end)+2048]);

    if max(IPTs(mi, :))>0
        segments(mi, 1) = find(IPTs(mi, :)>0, 1, 'first');
        segments(mi, 2) = find(IPTs(mi, :)>0, 1, 'last');
    else
        segments(mi, 1) = spikeRange(mi, 1);
        segments(mi, 2) = spikeRange(mi, 2);
    end
end
segments(segments(:,1)==0, :) = [];
startIndex = max(segments(:, 1));
endIndex = min(segments(:, 2));

% extend range based on spike range
if startIndex > min(spikeRange(:, 1))
    startIndex = min(spikeRange(:, 1));
end
if endIndex < max(spikeRange(:, 2))
    endIndex = max(spikeRange(:, 2));
end

% extend range based on reference signal
if startIndex > 1 && refSignal(startIndex)>0.1
    while (startIndex > 1 && refSignal(startIndex)>0.1)
        startIndex = startIndex - 1;
    end
    startIndex = startIndex - 2048*2;
end
if refSignal(endIndex)>0.1
    while (endIndex < signalLength && refSignal(endIndex)>0.1)
        endIndex = endIndex + 1;
    end
    endIndex = endIndex + 2048*2;
end

% contrain starti and endi with signal length
startIndex = max([1, startIndex]);
endIndex = min([endIndex, size(IPTs, 2), signalLength, length(SIG{2})]);

%% update SIG and IPTs as well as MUPulses
for k = 1:length(SIG(:))
    if ~isempty(SIG{k})
        SIG{k} = SIG{k}(startIndex:endIndex);
    end
end

IPTs = IPTs(:,startIndex:endIndex);
for mi = 1:MUCount
    MUPulses{mi} = MUPulses{mi}-startIndex+1;
end

signalLength = length(IPTs(1,:));
[~, ~, wSIG_m, ~] = extend_and_whiten(SIG, discardChannelsVec);

%% update IPTs and spikes for each motor unit
autoEd = auto_editor();
autoEd.disp_flag = 0;
autoEd = autoEd.set_para(para(1), para(2), para(3), para(4));
if plotFlag
    autoEd = autoEd.enable_plot();
end

if isfield(data, 'rule_based_editting')
    editInfo = data.rule_based_editting;
else
    editInfo = struct();
    editInfo.added_spikes = [];
    editInfo.removed_spikes = [];
end
shiftFlag = 1;
for muIndex = 1:MUCount
    if processFlag(muIndex) == 0 && isempty(test_mu_id)
        autoEd.disp_info(['Skip MU (process_flag): ', num2str(muIndex), '/', num2str(MUCount)]);
        continue;
    elseif ~isempty(test_mu_id) && sum(muIndex == test_mu_id)==0
        disp(['Skip MU:', num2str(muIndex)]);
        % editInfo(mu_index).flag = 0;
        continue;
    end
    autoEd.disp_info(['Processing MU ', num2str(muIndex), '/', num2str(MUCount)]);
    muTag = [nameItem{1}, '_MU', num2str(muIndex)];
    timeOutFlag = 0;
    % get ipt and spikes of
    iptTmp = IPTs(muIndex, :);
    spikes = MUPulses{muIndex};
    spikes_old = spikes;
    spikes_o = spikes;
    ipt_o = iptTmp;
    [autoEd, iptTmp, ~] = autoEd.update_ipts(spikes, wSIG_m);
    PNR_o = calculate_PNR(spikes, iptTmp, 2048);
    % reset for each motor unit
    autoEd = autoEd.reset();
    autoEd.adjust_enable = 0;
    autoEd.fr_std_limit = 2;
    autoEd.fr_min = 2;

    if length(spikes)<10
        editInfo(muIndex).added_spikes = 0;
        editInfo(muIndex).removed_spikes = 0;
        editInfo(muIndex).tuned_spikes = 0;
        editInfo(muIndex).cycles = 0;
        editInfo(muIndex).flag = 0;
        autoEd.disp_info(['Skip MU (few spikes): ', num2str(muIndex), '/', num2str(MUCount)]);
        % skip if number of spikes is less than 10
        continue;
    elseif PNR_o <= 10
        % skip MU if PNR is lower than 20
        editInfo(muIndex).added_spikes = 0;
        editInfo(muIndex).removed_spikes = 0;
        editInfo(muIndex).tuned_spikes = 0;
        editInfo(muIndex).cycles = 0;
        editInfo(muIndex).flag = 0;
        [autoEd, fig] = autoEd.plot_spikes_ipt_fr(iptTmp, spikes);
        sgtitle(['MU', num2str(muIndex), '; PNR: ', num2str(PNR_o)]);
        saveas(fig, [fileTag, '_MU', num2str(muIndex), '_skip.fig']);
        autoEd.disp_info(['Skip MU (low PNR): ', num2str(muIndex), '/', num2str(MUCount)]);
        continue;
    end

    % update spikes to maximize the pnr
    if shiftFlag
        pnrList = zeros(20, 4);
        for ti = 1:20
            shift = 10 - ti;
            spikes_n = spikes + shift;
            spikes_n(spikes_n<1) = [];
            spikes_n(spikes_n>signalLength) = [];
            [autoEd, ipt_n, ~] = autoEd.update_ipts(spikes, wSIG_m);
            spikes_ipt = ipt_n(spikes_n);
            pnrList(ti, 1) = shift;
            pnrList(ti, 2) = calculate_PNR(spikes_n, ipt_n, 2048);
            pnrList(ti, 3) = mean(spikes_ipt);
            pnrList(ti, 4) = std(spikes_ipt);
        end
        [PNR_n, i] = max(pnrList(:, 2));
        [ip_mean, i2] = max(pnrList(:, 3));
        shift = pnrList(i, 1);
        shift2 = pnrList(i2, 1);
        PNR_n2 = pnrList(i2, 2);
        %             if PNR_n2 > 30 || abs(shift2) >= 3
        if PNR_n2 > 35 && abs(shift2) >= 3
            %               shift2 = features{mu_index, 4};
            spikes = spikes + shift2;
            spikes(spikes<1) = [];
            spikes(spikes>signalLength) = [];
            [autoEd, ipt_n, ~] = autoEd.update_ipts(spikes, wSIG_m);
            spikes_old = spikes;
            MUPulses{muIndex} = spikes;
        elseif abs(shift) >= 3
            spikes = spikes + shift;
            spikes(spikes<1) = [];
            spikes(spikes>signalLength) = [];
            [autoEd, ipt_n, ~] = autoEd.update_ipts(spikes, wSIG_m);
            ipt_o = ipt_n;
            spikes_old = spikes;
            MUPulses{muIndex} = spikes;
        end
    end

    % iteratively update the spikes and ipts
    adjusted_cnt = 10;
    cycleIndex = 0;
    % max_cycles = 5;
    adjusted_cnt_hist = zeros(maxCycles, 1);
    [autoEd, iptTmp, ~] = autoEd.update_ipts(spikes, wSIG_m);
    spike_updated = cell(1,1);
    while cycleIndex <= maxCycles && adjusted_cnt > 0
        cycleIndex = cycleIndex + 1;
        autoEd.disp_info(['Cycles: ', num2str(cycleIndex), '/', num2str(maxCycles)]);

        % remove spikes with high firing rate and update the IPTs
        [autoEd, spikes, removed_spikes_fr] = autoEd.remove_high_fr_pre(iptTmp, spikes, frHighThreshold);
        while ~isempty(removed_spikes_fr)
            [autoEd, iptTmp, ~] = autoEd.update_ipts(spikes, wSIG_m);
            [autoEd, spikes, removed_spikes_fr] = autoEd.remove_high_fr_pre(iptTmp, spikes, frHighThreshold);
        end

        % remove spikes with low ipt and update the IPTs
        [autoEd, spikes, removed_spikes_ipt]= autoEd.remove_low_ipt(iptTmp, spikes, frHighThreshold);
        while ~isempty(removed_spikes_ipt)
            [autoEd, iptTmp, ~] = autoEd.update_ipts(spikes, wSIG_m);
            [autoEd, spikes, removed_spikes_ipt]= autoEd.remove_low_ipt(iptTmp, spikes, frHighThreshold);
        end

        % adjust spikes with surrounding peaks
        [autoEd, spikes, removed_spikes] = autoEd.adjust_spikes(iptTmp, spikes, frHighThreshold);
        while ~isempty(removed_spikes)
            [autoEd, iptTmp, ~] = autoEd.update_ipts(spikes, wSIG_m);
            [autoEd, spikes, removed_spikes] = autoEd.adjust_spikes(iptTmp, spikes, frHighThreshold);
        end

        if cycleIndex==1
            % refine spikes with a large margin for the first round
            [autoEd, spikes, updated_cnt] = autoEd.refine_spikes(iptTmp, spikes, 100);
        else
            [autoEd, spikes, updated_cnt] = autoEd.refine_spikes(iptTmp, spikes, 20);
        end

        [autoEd, iptTmp, ~] = autoEd.update_ipts(spikes, wSIG_m);
        % go through spikes and ipts forward and backward
        for di = 1:2
            [iptTmp, spikes] = flip_ipts_spikes(iptTmp, spikes);
            if di == 2
                [autoEd, iptTmp, ~] = autoEd.update_ipts(spikes, wSIG_m);
            end
            % add spikes with high ipt
            [autoEd, spikes, added_spikes_ipt] = autoEd.add_high_ipt(iptTmp, spikes, frHighThreshold);
            while ~isempty(added_spikes_ipt)
                [autoEd, spikes, added_spikes_ipt] = autoEd.add_high_ipt(iptTmp, spikes, frHighThreshold);
            end
            % remove spikes with hihg fr or spikes with low ipts
            [autoEd, spikes, removed_spikes] = autoEd.adjust_spikes(iptTmp, spikes, frHighThreshold);
            [autoEd, spikes, removed_spikes_fr] = autoEd.remove_high_fr_pre(iptTmp, spikes, frHighThreshold);
            [autoEd, spikes, removed_spikes_ipt]= autoEd.remove_low_ipt(iptTmp, spikes, frHighThreshold);
            [autoEd, spikes, added_spikes_ipt] = autoEd.add_first_spike(iptTmp, spikes, frHighThreshold);
        end

        % add potential spike if firing rate is low
        [autoEd, iptTmp, ~] = autoEd.update_ipts(spikes, wSIG_m);
        [autoEd, spikes, added_spikes_fr] = autoEd.add_low_fr(iptTmp, spikes, frHighThreshold);
        while ~isempty(added_spikes_fr)
            [autoEd, iptTmp, ~] = autoEd.update_ipts(spikes, wSIG_m);
            [autoEd, spikes, added_spikes_fr] = autoEd.add_low_fr(iptTmp, spikes, frHighThreshold);
        end

        % remove low ipt
        [autoEd, spikes, removed_spikes_ipt]= autoEd.remove_low_ipt(iptTmp, spikes, frHighThreshold);
        [autoEd, spikes, removed_spikes] = autoEd.remove_first_low_fr(iptTmp, spikes, frHighThreshold);
        [autoEd, spikes, removed_spikes] = autoEd.adjust_spikes(iptTmp, spikes, frHighThreshold);
        [autoEd, spikes, removed_spikes_fr] = autoEd.remove_high_fr_pre(iptTmp, spikes, frHighThreshold);
        [autoEd, spikes, added_spikes_fr] = autoEd.add_low_fr(iptTmp, spikes, frHighThreshold);

        % refine spikes
        interval = diff(spikes);
        interval_min = min(interval);
        margin_s = round(interval_min/3);
        updated_cnt = 10;
        cycles_u = 20;
        while cycles_u > 0 && updated_cnt > 0
            [autoEd, iptTmp, ~] = autoEd.update_ipts(spikes, wSIG_m);
            [autoEd, spikes, updated_cnt] = autoEd.refine_spikes(iptTmp, spikes, margin_s);
            cycles_u = cycles_u - 1;
        end

        % check number of adjusted spikes
        [autoEd, iptTmp, ~] = autoEd.update_ipts(spikes, wSIG_m);
        [added_spikes, removed_spikes, fine_tuned] = compare_mupulses(spikes_old, spikes);
        spikes_old = spikes;
        adjusted_cnt = length(added_spikes)+length(removed_spikes);
        adjusted_cnt_hist(cycleIndex) = adjusted_cnt;
        spike_updated{cycleIndex} = [removed_spikes, added_spikes];
        if length(spikes) < 10
            break;
        end

        % update parameters based on signal quality
        % these functions are for noisy signals
        autoEd.adjust_enable = 1;
        if PNR_n >= 35
            autoEd.fr_std_limit = 5;
            autoEd.fr_min = 0.5;
        elseif PNR_n > 30
            autoEd.fr_std_limit = 3;
            autoEd.fr_min = 1;
        elseif PNR_n > 25
            autoEd.fr_std_limit = 2;
            autoEd.fr_min = 2;
        else
            autoEd.fr_std_limit = 1;
            autoEd.fr_min = 5;
        end
    end

    if adjusted_cnt > 0 && cycleIndex > maxCycles
        autoEd.disp_info('Time out!');
        timeOutFlag = 1;
    end

    %% compare differences and save the difference
    [added_spikes, removed_spikes, fine_tuned] = compare_mupulses(MUPulses{muIndex}, spikes);
    % save removed spikes
    if plotFlag && ~isempty(removed_spikes)
        [~, fig] = autoEd.plot_spike_firing(ipt_o, MUPulses{muIndex}, removed_spikes);
        sgtitle(['MU', num2str(muIndex), ': removed']);
        % saveas(fig, [file_tag, '_MU', num2str(mu_index), '_removed.png']);
    end
    % save added spikes
    if plotFlag && ~isempty(added_spikes)
        [~, fig] = autoEd.plot_spike_firing(iptTmp, spikes, added_spikes);
        sgtitle(['MU', num2str(muIndex), ': added']);
        % saveas(fig, [file_tag, '_MU', num2str(mu_index), '_added.png']);
    end

    % save processed spike train
    iptTmp = iptTmp./max(iptTmp);
    if plotFlag
        [autoEd, iptTmp, ~] = autoEd.update_ipts(spikes, wSIG_m);
        [~, fig] = autoEd.plot_spikes_ipt_fr(ipt_o, spikes_o);
        [~, fig] = autoEd.plot_spikes_ipt_fr(iptTmp, spikes);
        sgtitle(['MU', num2str(muIndex), ' - ', MUIDs{muIndex}]);
        saveas(fig, [fileTag, '_MU', num2str(muIndex), '_2post.png']);
    end

    % summarize information
    editInfo(muIndex).added_spikes = length(added_spikes);
    editInfo(muIndex).removed_spikes = length(removed_spikes);
    editInfo(muIndex).tuned_spikes = length(fine_tuned);
    editInfo(muIndex).cycles = cycleIndex;
    editInfo(muIndex).flag = timeOutFlag;

    % display if the PNR is higher than 35
    PNR_n = calculate_PNR(spikes, iptTmp, fsamp);
    if PNR_n >= 35
        [~, fig] = autoEd.plot_spikes_ipt_fr(iptTmp, spikes);
        sgtitle(['MU', num2str(muIndex), '; PNR: ', num2str(PNR_n)]);
        % saveas(fig, [file_tag, '_MU', num2str(mu_index), '_2post.fig']);
    end

    % update spike train and IPT
    MUPulses{muIndex} = spikes;
    IPTs(muIndex, 1:length(iptTmp)) = iptTmp;
    PNR(muIndex) = PNR_n;
    clc;
    close all;
end

%% post process
if ~isempty(filePath)
    cd(filePath);
else
    cd('..')
end

% update spikes (MUPulses) and ipts (IPTs)
for mi = 1:length(MUPulses)
    MUPulses{mi} = MUPulses{mi}+startIndex-1;
end
IPTs_t = data.IPTs;
if size(IPTs, 2) < size(IPTs_t, 2)
    IPTs_t(:, startIndex:endIndex) = IPTs(:, 1:endIndex-startIndex+1);
else
    IPTs_t = IPTs(:, 1:signalLength);
end

% save results
data.MUPulses = MUPulses;
data.IPTs(:, 1:size(IPTs_t, 2)) = IPTs_t;
data.rule_based_editting = editInfo;
data.PNR = PNR;
if fileName(1) == 'a'
    fileNameNew = fileName;
elseif fileName(1) == 'r'
    fileNameNew = ['a', fileName];
else
    fileNameNew = ['ar', fileName];
    %     fileName_new = fileName;
end
save(fileNameNew, '-struct', 'data');
if ~isempty(filePath)
    fileNameNew = [filePath, '\', fileNameNew];
end
cd(rootFolder);
end