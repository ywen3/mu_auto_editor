classdef auto_editor
%% editing rules for motor unit refinement using height of IPT and FR
% Rule 1: remove spikes with high firing rate
% Rule 2: remove spikes with low ipt
% Rule 3: add spikes with high ipt
% Rule 4: add spikes filling low firing rate
% 
% Author: Yue Wen at Shirley Ryan AbilityLab, 2022.
%
    properties
        removed_spikes;
        added_spikes;
        adjusted_spikes;
        low_fr_ignored;
        low_fr_processed;
        high_fr_ignored;
        high_fr_processed;
        debug_flag;
        abnormal_identified;
        lfr_spikes;
        dlfr_spikes;
        hfr_spikes;
        dhfr_spikes;
        hlfr_spikes;
        tuned_spikes;
        flag_hist;
        added_count;
        removed_count;

        % video record variables
        video_handle;
        edit_type;
        record_type;
        video_enable;

        save_ignored_flag;
        save_processed_flag;

        % ratio
        ipt_ratio_r_g;
        ipt_ratio_r_l;
        ipt_ratio_a_g;
        ipt_ratio_a_l;
        fr_ratio_limit;
        fr_std_limit;
        fr_min;
        ipt_threshold_stats;
        non_spike_ipt_mean;
        spike_ipt_mean;
        non_spike_ipt_std;
        spike_ipt_std;
        
        %
        ipt_check;
        fig_format;
        adjust_enable;
        save_flag;

        %
        disp_flag;
    end

    methods

        function obj = auto_editor()
            obj.low_fr_ignored = 1;
            obj.low_fr_processed = 1;
            obj.high_fr_ignored = 1;
            obj.high_fr_processed = 1;
            obj.abnormal_identified = 1;
            obj.debug_flag = 0;
            obj.added_count = 0;
            obj.removed_count = 0;
            obj.edit_type = -2;      % 1 for high ipt; 2 for low fr; -1 for low ipt; -2 for high fr
            obj.record_type = 0;
            obj.adjust_enable = 0;

            obj.save_processed_flag = 1;
            obj.save_ignored_flag = 0;
            obj.video_enable = 0;
            
            obj.ipt_threshold_stats = 0;
            obj.non_spike_ipt_mean = 0.35;
            obj.spike_ipt_mean = 0.75;
            obj.non_spike_ipt_std = 0.1;
            obj.spike_ipt_std = 0.1;

            obj.ipt_ratio_r_g = 0.2;
            obj.ipt_ratio_r_l = 0.3;

            obj.ipt_ratio_a_g = 0.4;
            obj.ipt_ratio_a_l = 0.4;

            obj.fr_ratio_limit = 2.5;
            obj.fr_min = 1;  % 1 Hz or 2 Hz
            obj.fr_std_limit = 2;

            obj.ipt_check = 0;
            obj.save_flag = 0;
            obj.fig_format = '.fig';
            obj.disp_flag = 1;
        end

        function autoEd = reset(autoEd)
            autoEd.ipt_threshold_stats = 0;
            autoEd.non_spike_ipt_mean = 0.35;
            autoEd.spike_ipt_mean = 0.75;
            autoEd.non_spike_ipt_std = 0.1;
            autoEd.spike_ipt_std = 0.1;
            autoEd.fr_min = 2;          % 1 Hz or 2 Hz
            autoEd.fr_ratio_limit = 1.5;
            autoEd.fr_std_limit = 2;
            autoEd.adjust_enable = 0;
            autoEd.ipt_check = 0;
        end

        function autoEd = set_para(autoEd, ratio_r_g, ratio_r_l, ratio_a_g, ratio_a_l)
            % remove spike if lower than both global and local threshold
            autoEd.ipt_ratio_r_g = ratio_r_g;
            autoEd.ipt_ratio_r_l = ratio_r_l;

            % add spike if higher than both global and local threshold
            autoEd.ipt_ratio_a_g = ratio_a_g;
            autoEd.ipt_ratio_a_l = ratio_a_l;
        end

        function autoEd = video_record_init(autoEd, fileName)
            vidobj = VideoMaker(fileName);
            autoEd.video_handle = vidobj;
            autoEd.video_enable = 1;
            autoEd.debug_flag = 1;
            autoEd.record_type = 0;
            autoEd.save_flag = 1;
        end

        function autoEd = video_record_write(autoEd, fig, type)
            if autoEd.video_enable == 1
                if type == autoEd.record_type
                    autoEd.video_handle = autoEd.video_handle.writeFrame(fig, 10);
                elseif autoEd.record_type == 0
                    autoEd.video_handle = autoEd.video_handle.writeFrame(fig, 10);
                end
            end
        end

        function autoEd = video_record_close(autoEd)
            autoEd.video_handle.closeVideo();
        end

        function autoEd = enable_plot(autoEd)
            autoEd.debug_flag = 1;
        end

        %% Rule 1: remove spikes if ipt is low
        % remove spikes with low ipt
        % remove spikes if two consecutive spikes has low ipts
        % if ipt is lower than 30% of mean ipts of surrounding spikes
        function [autoEd, spikes, removed_spikes]= remove_low_ipt(autoEd, ipt_t, spikes, fr_high_threshold)
            % skip if number of spikes is less than 10
            if length(spikes) < 10
                removed_spikes = [];
                return;
            end

            % S1. get firing rate
            fsamp = 2048;
            fr = fsamp./diff(spikes);
            spikes(fr==inf) = [];
            fr = fsamp./diff(spikes);
            fr = [fr, fr(end)];

            % S2. get ipt threshold
            ipt_spikes = ipt_t(spikes);
            ipt_mean = mean(ipt_spikes);

            % initialize parameters
            spike_cnt = length(spikes);
            remove_flag = spikes*0;
            editing_hist = struct();
            removed_cnt = 0;

            %% remove consecutive spikes if both have low ipt
            for k = 1:spike_cnt-2
                % get spike index, firing rate, and ipt
                spike_cur = spikes(k);
                spike_post = spikes(k+1);
                ipt_cur = ipt_t(spike_cur);
                ipt_post = ipt_t(spike_post);

                % R1: remove spikes if two consecutive spikes has ipt lower than the global threshold
                ipt_threshold_g = ipt_mean*autoEd.ipt_ratio_r_g;
                flag_1 = ipt_cur < ipt_threshold_g;
                flag_2 = ipt_post < ipt_threshold_g;
                flag_sum = flag_1 + flag_2;

                if flag_sum == 2
                    potential_spike_display(ipt_t, spikes, [spike_cur, spike_post], 2, 'Low ipt');
                    remove_flag(k) = 3;
                    remove_flag(k+1) = 3;
                    removed_cnt = removed_cnt + 2;
                    autoEd.write_csv(['Removing double low ipt - index: ', num2str(spike_cur), '; ipt: ', num2str(ipt_cur)]);
                end
            end

            %% remove a spike if its ipt is low
            spikes_t = [spikes(1), spikes, spikes(end)];
            for k = 2:spike_cnt-1
                if remove_flag(k) > 0
                    continue;
                end

                % get spike index and firing rate
                spike_cur = spikes(k);
                ipt_cur = ipt_t(spike_cur);

                % remove low ipts
                nearby_index = [k-2, k-1, k+1, k+2];
                nearby_ipt_mean = mean(ipt_t(spikes_t(nearby_index+1)));
                ipt_threshold_g = ipt_mean*autoEd.ipt_ratio_r_g;
                ipt_threshold_l = nearby_ipt_mean*autoEd.ipt_ratio_r_l;
                ipt_threshold = min([ipt_threshold_g, ipt_threshold_l]);
                ipt_low_flag = ipt_cur < ipt_threshold;

                % display removed spikes
                if ipt_low_flag
                    remove_flag(k) = 1;
                    removed_cnt = removed_cnt + 1;
                    % fig = potential_spike_display(ipt_t, spikes, spike_cur, 2, 'Low ipt');
                    % autoEd.write_csv(['Removing low ipt - index: ', num2str(spike_cur), '; ipt: ', num2str(ipt_cur)]);
                    % if autoEd.save_flag && autoEd.debug_flag
                    %     saveas(fig, ['removed_low_ipt_', num2str(autoEd.removed_count), autoEd.fig_format]);
                    %     autoEd.removed_count = autoEd.removed_count + 1;
                    %     autoEd = autoEd.video_record_write(fig, -1);
                    % end
                end

                % save data for debug purpose
                editing_hist(k).spike = spikes(k);
                editing_hist(k).fr = fr(k);
                editing_hist(k).ipt = ipt_cur;
                editing_hist(k).ipt_threshold = ipt_threshold;
                editing_hist(k).flag_1 = ipt_low_flag;
            end

            % update results
            removed_spikes = spikes(remove_flag>0);
            autoEd.removed_spikes = [autoEd.removed_spikes, removed_spikes];
            autoEd.flag_hist = editing_hist;
            spikes(remove_flag>0) = [];
            if removed_cnt
                autoEd.disp_info(['Spikes removed (low ipt):', num2str(removed_cnt)]);
            end
        end


        %% Rule 2: Remove spikes with high firing rate
        function [autoEd, spikes, removed_spikes] = remove_high_fr_pre(autoEd, ipt_t, spikes, fr_high_threshold)
            % S1. get firing rate threshold
            fsamp = 2048;
            fr = fsamp./diff(spikes);
            spikes(fr==inf) = [];               % remove repeat spikes
            fr = fsamp./diff(spikes);
            fr = [fr, fr(end)];
            % if mean(fr) > 5 % if fr_mean greater than 5, remove fr lower than 5 Hz
            %     fr_g5 = fr(fr>=5);
            % else
            %     fr_g5 = fr(fr>=2);
            % end
            % fr_mean = mean(fr_g5);

            % check fr standard deviation of 10-consecutive spikes
            % si_index = 1:5:length(fr)-10;
            % fr_std_hist = si_index*0;
            % for k = 1:length(si_index)
            %     fr_std_hist(k) = std(fr(si_index(k):si_index(k)+10));
            % end
            % fr_std = median(fr_std_hist);
            fr_hard_threshold = fr_high_threshold;

            % S2. initialize variables
            remove_flag = spikes*0;
            spike_cnt = length(spikes);
            removed_cnt = 0;
            ipt_mean = mean(ipt_t(spikes));
            ipt_std = std(ipt_t(spikes));
            ipt_threshold = ipt_mean - 2*ipt_std;
            % to avoid negtive ipt_threshold
            ipt_threshold = max([ipt_threshold, ipt_mean*0.3]);

            editing_hist = struct();
            spike_index = 1;
            fr_std_hist = spikes*0;
            % S3. go through each spike
            seg_len = 10;
            seg_len_half = seg_len/2;
            % set to 0 if no constraints
            std_flag_threshold = 0.25;
            for k = 2:spike_cnt-1
                seg_si = max([1, k-seg_len_half]);
                seg_ei = min([k+seg_len_half, spike_cnt-1]);
                fr_seg = fr(seg_si:seg_ei);
                fr_seg(fr_seg<3) = [];
                fr_seg(fr_seg>50) = [];
                fr_seg_median = median(fr_seg);
                fr_seg_std = std(fr_seg);
                fr_std_hist(k) = fr_seg_std;
                
                % if the variance is too high, reduce the fr boundary
                if fr_seg_std >= fr_seg_median*0.3
                    fr_seg_std = fr_seg_median*0.3;
                    fr_std_limit_l = autoEd.fr_std_limit-1;
                else
                    fr_std_limit_l = autoEd.fr_std_limit;
                end

                fr_threshold = fr_seg_median + fr_seg_std*fr_std_limit_l;
                fr_threshold = max([fr_threshold, fr_seg_median*1.5]);

                % initialize
                spike_cur = spikes(k);
                spike_post = spikes(k+1);
                ipt_cur = ipt_t(spike_cur);
                ipt_post = ipt_t(spike_post);
                fr_cur = fr(k);

                % fr_pre_post_flag = 0;
                fr_pre = fr(k-1);
                fr_post = fr(k+1);
                % fr_pre_post_ratio = fr_pre/fr_post;
                % fr_pre_post_flag = fr_pre_post_ratio > 0.75 && fr_pre_post_ratio < 1.25;

                % check high firing rate
                remove_high_fr_flag = 0;
                ipt_flag = 1;
                % high firing rate > 
                if fr_cur >= fr_hard_threshold || fr_cur >= fr_seg_median*2
                    if ipt_cur > ipt_post
                        remove_high_fr_flag = 1;
                    elseif ipt_post > ipt_cur
                        remove_high_fr_flag = 2;
                    else
                        remove_high_fr_flag = -1;
                    end
                    % if remove_high_fr_flag >= 1 && fr_pre_post_flag == 1
                    %     potential_spike_display(ipt_t, spikes, spike_cur, 10, 'High firing rate:: potential');
                    %     autoEd.write_csv(['False high fr - high fr threshold - normal ratio: ', num2str(fr_pre_post_flag)]);
                    % end
                elseif fr_cur >= fr_threshold
                    if ipt_cur >= ipt_post*1.5
                        remove_high_fr_flag = 1;
                        ipt_flag = ipt_post <= ipt_threshold;
                    elseif ipt_post >= ipt_cur*1.5
                        remove_high_fr_flag = 2;
                        ipt_flag = ipt_cur <= ipt_threshold;
                    else
                        remove_high_fr_flag = -2;
                    end
                    % if remove_high_fr_flag >= 1 && ipt_flag && fr_pre_post_flag == 1
                    %     potential_spike_display(ipt_t, spikes, spike_cur, 10, 'High firing rate (current)');
                    %     autoEd.write_csv(['False high fr - low fr threshold -  ratio: ', num2str(fr_pre_post_flag)]);
                    % end
                end

                % display
                if remove_high_fr_flag == 1 && ipt_flag
                    std_flag = std([fr_pre, fr_cur, fr_post, fr_seg_median])/fr_seg_median;
                    if std_flag > std_flag_threshold
                        remove_flag(k+1) = 1;
                        fr(k+1) = 0;
                        removed_cnt = removed_cnt + 1;
                        % autoEd.write_csv(['Remove high fr - index: ',num2str(spike_post), '; ipt: ', num2str(ipt_post)]);
                    else
                        % autoEd.write_csv(['False high fr - low std 1 -', num2str(std_flag)]);
                    end
                    % fig = potential_spike_display(ipt_t, spikes, [spike_cur, spike_post], 2, 'High firing rate (remove 2nd)');
                    % if autoEd.save_flag && autoEd.debug_flag
                    %     saveas(fig, ['removed_hfrp_', num2str(autoEd.removed_count), autoEd.fig_format]);
                    %     autoEd.removed_count = autoEd.removed_count + 1;
                    %     autoEd = autoEd.video_record_write(fig, -2);
                    % end
                elseif remove_high_fr_flag == 2 && ipt_flag
                    std_flag = std([fr_pre, fr_cur, fr_post, fr_seg_median])/fr_seg_median;
                    if std_flag > std_flag_threshold
                        remove_flag(k) = 1;
                        fr(k) = 0;
                        removed_cnt = removed_cnt + 1;
                        % autoEd.write_csv(['Remove high fr - index: ',num2str(spike_cur), '; ipt: ', num2str(ipt_cur)]);
                    else
                        % autoEd.write_csv(['False high fr - low std 2 -', num2str(std_flag)]);
                    end
                    % fig = potential_spike_display(ipt_t, spikes, [spike_cur, spike_post], 2, 'High firing rate (remove 1st)');
                    % if autoEd.save_flag && autoEd.debug_flag
                    %     saveas(fig, ['removed_hfrp_', num2str(autoEd.removed_count), autoEd.fig_format]);
                    %     autoEd.removed_count = autoEd.removed_count + 1;
                    %     autoEd = autoEd.video_record_write(fig, -2);
                    % end
                elseif remove_high_fr_flag == -1
                    fig = potential_spike_display(ipt_t, spikes, spike_cur, 0, 'High firing rate (ignore)');
                    if autoEd.save_flag && autoEd.debug_flag
                        saveas(fig, ['ignore_hfrp_', num2str(autoEd.removed_count), autoEd.fig_format]);
                        autoEd = autoEd.video_record_write(fig, -2);
                    end
                elseif remove_high_fr_flag == -2
                    fig = potential_spike_display(ipt_t, spikes, spike_cur, 0, 'High firing rate (ignore)');
                    if autoEd.save_flag && autoEd.debug_flag
                        saveas(fig, ['ignore_hfrp_', num2str(autoEd.removed_count), autoEd.fig_format]);
                        autoEd = autoEd.video_record_write(fig, -2);
                    end
                end

                % save data for debug purpose
                editing_hist(spike_index).spike = spike_cur;
                editing_hist(spike_index).ipt_cur = ipt_cur;
                editing_hist(spike_index).ipt_post = ipt_post;
                editing_hist(spike_index).remove_flag = remove_high_fr_flag;
                spike_index = spike_index + 1;
            end

            % update results
            removed_spikes = spikes(remove_flag==1);
            spikes(remove_flag==1) = [];
            autoEd.removed_spikes = [autoEd.removed_spikes, removed_spikes];
            if removed_cnt
                autoEd.disp_info(['Spikes removed (high fr pre):', num2str(removed_cnt)]);
            end
        end


        %% Rule 3: add spike based on ipt and fr threshold
        % the purpose is to add potential spikes based on ipt and firing rate
        function [autoEd, spikes, added_spikes_t] = add_high_ipt(autoEd, ipt_t, spikes, fr_high_threshold)
            % S1. get firing rate threshold
            fsamp = 2048;
            fr = fsamp./diff(spikes);
            spikes(fr==inf) = [];
            fr = fsamp./diff(spikes);
            fr = [fr, fr(end)];
            if mean(fr) > 5
                fr_g5 = fr(fr>=5);
            else
                fr_g5 = fr(fr>=2);
            end
            fr_mean = mean(fr_g5);
            fr_threshold = min([fr_mean*2, fr_high_threshold]);
            % fr_threshold_low = autoEd.fr_min;

            % S2. get ipt threshold
            ipt_spikes = ipt_t(spikes);
            ipt_mean = mean(ipt_spikes);

            % initialize parameters
            window_size = round(fsamp/fr_threshold);
            if isnan(window_size)
                window_size = round(fsamp/50);
            end
            added_cnt = 0;
            added_spikes_t = [];        % store potential spikes
            spike_index = 1;
            editing_hist = struct();
            start_index = find(ipt_t>ipt_mean*0.1, 1, 'first');
            % start_index = start_index - round(rand(1)*window_size/2); % randomize the window segmentation
            start_index = max([start_index, 100]);
            for k = start_index:window_size:length(ipt_t)-window_size
                % use window size
                ipt_seg = ipt_t(k:k+window_size-1);
                % skip if maximum of ipt in a window is low
                if max(ipt_seg) <= ipt_mean*0.1
                    continue;
                end

                % find potential spike
                potential_spikes = find_potential_spike(ipt_t, k-10, k+window_size+9, 5);
                potential_spike = potential_spikes(1);
                potential_spike_ipt = ipt_t(potential_spike);

                % get mean ipt of nearby spikes - compare potential spike with existing spikes
                spikes_t = sort([spikes, added_spikes_t]);
                spike_margins = abs(spikes_t - potential_spike);

                % get nearby spikes and nearby_ipt_mean
                [sv, si] = sort(spike_margins);
                if length(spikes_t) >= 4
                    nearby_spikes = spikes_t(si(1:4));
                else
                    nearby_spikes = spikes_t;
                end
                nearby_ipt_mean = mean(ipt_t(nearby_spikes));

                % ipt threshold
                ipt_threshold_g = ipt_mean*autoEd.ipt_ratio_a_g;
                ipt_threshold_l = nearby_ipt_mean*autoEd.ipt_ratio_a_l;
                ipt_threshold = max([ipt_threshold_g, ipt_threshold_l]);
                ipt_threshold_low = min([ipt_threshold_g, ipt_threshold_l]);

                % potential spike firing rate
                spike_margin = sv(1);
                potential_spike_fr = fsamp/spike_margin;

                % add spikes
                add_flag = 0;
                ipt_accept = 0;
                if spike_margin > 5
                    if potential_spike_ipt > ipt_threshold && potential_spike_fr < fr_high_threshold
                        ipt_accept = 1;
                        add_flag = 1;
                    elseif potential_spike_ipt > ipt_threshold_low && potential_spike_fr < fr_threshold && potential_spike_fr > autoEd.fr_min
                        ipt_accept = 1;
                        add_flag = 1;
                    end
                else
                    add_flag = -1;
                end

                % add potential spikes and display results
                if add_flag > 0
                    added_spikes_t = [added_spikes_t, potential_spike];
                    added_cnt = added_cnt + 1;
                    % autoEd.write_csv(['Adding high IPT - index: ', num2str(potential_spike), '; ipt: ', num2str(potential_spike_ipt)]);
                    % fig = potential_spike_display(ipt_t, spikes_t, potential_spike, 1, 'High ipt');
                    % if autoEd.save_flag && autoEd.debug_flag
                    %     saveas(fig, ['added_hipt_', num2str(autoEd.added_count), autoEd.fig_format]);
                    %     autoEd.added_count = autoEd.added_count + 1;
                    %     autoEd = autoEd.video_record_write(fig, 1);
                    % end
                end

                % record information
                editing_hist(spike_index).potential_spike = potential_spike;
                editing_hist(spike_index).ipt = potential_spike_ipt;
                editing_hist(spike_index).fr = potential_spike_fr;
                editing_hist(spike_index).ipt_accept = ipt_accept;
                editing_hist(spike_index).add_flag = add_flag;
                editing_hist(spike_index).ipt_global = potential_spike_ipt > ipt_threshold_g;
                editing_hist(spike_index).ipt_local = potential_spike_ipt > ipt_threshold_l;
                spike_index = spike_index + 1;
            end

            % add spikes and return
            spikes = sort([spikes, added_spikes_t]);
            autoEd.flag_hist = editing_hist;
            if added_cnt
                autoEd.disp_info(['Spikes added (high ipt):', num2str(added_cnt)]);
            end
        end


        %% Rule 4: add spike if the firing rate is low
        % the purpose is to add potential spikes based firing rate
        function [autoEd, spikes, added_spikes]= add_low_fr(autoEd, ipt_t, spikes, fr_high_threshold)
            % skip if number of spikes is lower than 10
            if length(spikes) < 10
                added_spikes = [];
                return;
            end

            % S1. get firing rate threshold
            added_cnt = 0;
            fsamp = 2048;
            fr = fsamp./diff(spikes);
            spikes(fr==inf) = [];
            fr = fsamp./diff(spikes);
            fr = [fr, fr(end)];
            if mean(fr) > 5
                fr_g5 = fr(fr>=5);
            else
                fr_g5 = fr(fr>=2);
            end
            fr_mean = mean(fr_g5);
            fr_std = std(fr_g5);

            % S2. get ipt threshold
            ipt_spikes = ipt_t(spikes);
            ipt_mean = mean(ipt_spikes);
            ipt_std = std(ipt_spikes);

            % initialize parameters
            editing_hist = struct();
            added_spikes = [];
            spikes_t = spikes;
            spike_cnt = length(spikes_t);
            seg_len = 10;
            seg_len_half = seg_len/2;
            for k = 1:spike_cnt-1
                % get start index and ending index of each segment
                seg_si = max([1, k-seg_len_half]);
                seg_ei = min([k+seg_len_half, spike_cnt-1]);

                % refine firing rate of a segment of spikes
                fr_seg = fr(seg_si:seg_ei);
                fr_seg_len = length(fr_seg);
                if  fr_seg_len <= seg_len
                    if k-seg_len < 1
                        fr_seg_t = ones(1, seg_len+1)*fr_seg(1);
                        fr_seg_t(1:fr_seg_len) = fr_seg;
                        fr_seg = fr_seg_t;
                    elseif k+seg_len > spike_cnt
                        fr_seg_t = ones(1, seg_len+1)*fr_seg(end);
                        fr_seg_t(1:fr_seg_len) = fr_seg;
                        fr_seg = fr_seg_t;
                    end
                end
                fr_seg(fr_seg<3) = [];
                fr_seg(fr_seg>50) = [];
                fr_seg_median = max([mean(fr_seg), median(fr_seg)]);
                fr_seg_std = min([std(fr_seg), 2*fr_std]);

                % get fr_threshold_seg
                fr_threshold_seg = fr_seg_median + fr_seg_std*autoEd.fr_std_limit;
                fr_threshold_seg = max([fr_threshold_seg, fr_seg_median*1.5]);

                % get spike index and firing rate
                spike_cur = spikes(k);
                spike_post = spikes(k+1);
                fr_cur = fr(k);
                fr_post = fr(k+1);

                % check if fr is low
                flag_fr = fr_cur>autoEd.fr_min && fr_post>autoEd.fr_min; %
                fr_ratio_t = fr_cur/fr_seg_median;
                flag_ratio = fr_ratio_t < 0.7;    % corresponding to 30% of variance tolerance
                flag_low = flag_fr + flag_ratio;

                fr_accept = 0;
                ipt_accept = 0;
                std_flag = 1;
                mean_flag = 1;
                fr_ratio_accept = 0;

                % firing rate and ipt requirements
                if flag_low == 2
                    % find potential spike
                    if fr_ratio_t < 0.4
                        potential_spikes = find_potential_spike_loc(ipt_t, spike_cur, spike_post, 0.33);
                    else
                        potential_spikes = find_potential_spike(ipt_t, spike_cur, spike_post, 5);
                    end
                    potential_spike = potential_spikes(1);

                    % Rule 1: check firing rate of the potential spike
                    fr_potential_pre = fsamp/(potential_spike-spike_cur);
                    fr_potential_post = fsamp/(spike_post-potential_spike);
                    fr_ratio = fr_potential_pre/fr_potential_post;
                    fr_accept = fr_potential_pre < fr_threshold_seg & fr_potential_post < fr_threshold_seg;

                    % check ipt of the potential spike
                    ipt_potential = ipt_t(potential_spike);
                    spike_margins = abs(spikes - potential_spike);
                    [sv, si] = sort(spike_margins);
                    if length(si) >= 4
                        nearby_spikes = spikes(si(1:4));
                    else
                        nearby_spikes = spikes_t(si);
                    end
                    % mean ipts of surrounding four spikes
                    nearby_ipt_mean = mean(ipt_t(nearby_spikes));
                    % get ipt threshold
                    ipt_threshold_l = nearby_ipt_mean*autoEd.ipt_ratio_a_l;
                    ipt_threshold_g = ipt_mean*autoEd.ipt_ratio_a_g;
                    ipt_threshold = max([ipt_threshold_l, ipt_threshold_g])*0.5;
                    ipt_accept = ipt_potential >= ipt_threshold;

                    % Rule 2: other assessment
                    fr_ratio_accept = fr_ratio >= 0.7 && fr_ratio <= 1.3;
                    ipt_threshold_low = min([ipt_threshold_l, ipt_threshold_g])*0.5;
                    ipt_accept_low = ipt_potential >= ipt_threshold_low;

                    std_flag = std([fr_post, fr_potential_pre, fr_potential_post]) > std([fr_cur, fr_post]);
                    mean_flag = abs(fr_seg_median - (fr_potential_pre+fr_potential_post)/2) > abs(fr_seg_median - fr_cur);
                end

                % add spike if fullfill firing rate and ipt requirements
                if fr_accept && ipt_accept
                    added_spikes = [added_spikes, potential_spike];
                    added_cnt = added_cnt + 1;
                    spikes_t = sort([spikes, added_spikes]);
                    % if autoEd.save_flag && autoEd.debug_flag
                    %     fig = potential_spike_display(ipt_t, spikes_t, potential_spike, 1, 'Low firing rate');
                    %     saveas(fig, ['added_lfr_', num2str(autoEd.added_count), autoEd.fig_format]);
                    %     autoEd.added_count = autoEd.added_count + 1;
                    %     autoEd = autoEd.video_record_write(fig, 2);
                    % end
                    % autoEd.write_csv(['Adding low fr - index: ', num2str(potential_spike), '; ipt: ', num2str(ipt_potential)]);
                elseif fr_ratio_accept && ipt_accept_low
                    added_spikes = [added_spikes, potential_spike];
                    added_cnt = added_cnt + 1;
                    spikes_t = sort([spikes, added_spikes]);
                    % fig = potential_spike_display(ipt_t, spikes_t, potential_spike, 11, 'Low firing rate');
                    % autoEd.write_csv(['Adding spike with ratio and ipt global - index: ', num2str(potential_spike), '; ipt: ', num2str(ipt_potential)]);
                elseif ipt_accept && (std_flag == 0 || mean_flag == 0)
                    % to delete
                    added_spikes = [added_spikes, potential_spike];
                    added_cnt = added_cnt + 1;
                    fig = potential_spike_display(ipt_t, spikes_t, potential_spike, 1, 'Low firing rate');
                    spikes_t = sort([spikes, added_spikes]);
                    autoEd.write_csv(['Adding spike with std ratio and mean ratio - index: ', num2str(potential_spike), '; ipt: ', num2str(ipt_potential)]);
                end

                % save data for debug purpose
                editing_hist(k).spike = spikes(k);
                editing_hist(k).fr = fr(k);
                editing_hist(k).flag_fr = flag_fr;
                editing_hist(k).flag_ratio  = flag_ratio;
                editing_hist(k).flag_sum  = flag_low;
                editing_hist(k).fr_accept = fr_accept;
                editing_hist(k).ipt_accept = ipt_accept;
                editing_hist(k).fr_ratio_t = fr_ratio_t;
            end

            % add spikes
            spikes = sort([spikes, added_spikes]);
            autoEd.added_spikes = [autoEd.added_spikes, added_spikes];
            if added_cnt
                autoEd.disp_info(['Spikes added (low fr):', num2str(added_cnt)]);
            end
        end

        
        %% Rule: add spike if the firing rate match
        function [autoEd, spikes, added_spikes_t] = add_first_spike(autoEd, ipt_t, spikes, fr_high_threshold)
            %% Add spike at begining if fulfill requirements
            % get firing rate threshold
            added_cnt = 0;
            fsamp = 2048;
            fr = fsamp./diff(spikes);
            spikes(fr==inf) = [];
            fr = fsamp./diff(spikes);
            added_spikes_t = [];

            if length(spikes) < 10
                return;
            end

            % get ipt mean
            ipt_spikes = ipt_t(spikes);
            ipt_mean = mean(ipt_spikes);
            spike_cur = spikes(1);
            fr_cur = fr(1);

            % find potential spike
            margin = round(2048/autoEd.fr_min);  % use the max margin corresponding to 1 Hz
            potential_spikes = find_potential_spike(ipt_t, max([1, spike_cur-margin]), spike_cur, 5);
            if length(potential_spikes) > 2
                margin = round(2*2048/fr_cur);
                potential_spikes = find_potential_spike(ipt_t, max([1, spike_cur-margin]), spike_cur, 5);
            end
            potential_spike = potential_spikes(1);
            fr_potential = fsamp/(spike_cur-potential_spike);
            ipt_potential = ipt_t(potential_spike);

            % get nearby spikes
            if length(spikes) >= 4
                nearby_spikes = spikes(1:4);
            else
                nearby_spikes = spikes;
            end
            nearby_ipt_mean = mean(ipt_t(nearby_spikes));

            % get ipt threshold
            ipt_threshold_g = ipt_mean*autoEd.ipt_ratio_a_g;
            ipt_threshold_l = nearby_ipt_mean*autoEd.ipt_ratio_a_l;
            ipt_threshold = max([ipt_threshold_l, ipt_threshold_g])*0.5;

            % add if fulfill firing rate
            fr_ratio = fr_potential/fr_cur;
            if fr_ratio < 1.5 && fr_ratio > 0.35 && fr_potential >= autoEd.fr_min
                if ipt_potential > ipt_threshold && length(potential_spikes) <= 2
                    added_spikes_t = potential_spike;
                    added_cnt = added_cnt + 1;
                    potential_spike_display(ipt_t, spikes, potential_spike, 10, 'Low firing rate');
                    autoEd.disp_info('Matched firing rate!');
                    autoEd.write_csv(['Adding first spike - ', num2str(potential_spike)]);
                elseif length(potential_spikes) == 1
                    potential_spike_display(ipt_t, spikes, potential_spike, 10, 'Low firing rate');
                end
            elseif ipt_potential > ipt_threshold && length(potential_spikes) <= 2
                if length(fr) > 10
                    fr_seg = fr(1:10);
                    fr_seg_median = max([mean(fr_seg), median(fr_seg)]);
                    fr_seg_std = std(fr_seg);
                    fr_threshold = fr_seg_median + fr_seg_std*autoEd.fr_std_limit;
                    if fr_potential < fr_threshold && fr_potential >= autoEd.fr_min
                        potential_spike_display(ipt_t, spikes, potential_spike, 11, 'First low firing rate');
                        added_spikes_t = potential_spike;
                        added_cnt = added_cnt + 1;
                        autoEd.write_csv(['Adding first spike - Rule 2 -', num2str(potential_spike)]);
                    end
                end
            end

            % add spikes
            spikes = sort([spikes, added_spikes_t]);
            autoEd.added_spikes = [autoEd.added_spikes, added_spikes_t];
            if added_cnt
                autoEd.disp_info(['Spikes added (first spike):', num2str(added_cnt)]);
            end
        end

        %% Rule: remove first spike if fr and ipt are low
        function [autoEd, spikes, removed_spikes]= remove_first_low_fr(autoEd, ipt_t, spikes, fr_high_threshold)
            % remove spikes with high FR
            % get firing rate threshold
            fsamp = 2048;
            fr = fsamp./diff(spikes);
            spikes(fr==inf) = [];
            fr = fsamp./diff(spikes);
            fr = [fr, fr(end)];
            if mean(fr) > 5
                fr_g5 = fr(fr>=5);
            else
                fr_g5 = fr(fr>=2);
            end
            fr_mean = mean(fr_g5);
            fr_std = std(fr_g5);

            % get ipt threshold
            ipt_spikes = ipt_t(spikes);
            ipt_mean = mean(ipt_spikes);
            removed_spikes = [];

            % remove first spike
            if (fr(1)<2)
                if ipt_spikes(1) < ipt_mean*0.5
                    removed_spikes = [spikes(1)];
                    fig = potential_spike_display(ipt_t, spikes, spikes(1), 10, 'First low fr');
                    autoEd.write_csv(['Removing first low - index: ', num2str(spikes(1)), '; ipt: ', num2str(ipt_spikes(1))]);
                    spikes(1) = [];
                end
            end

            % remove first low firing rate and last low firing rate
            % fr = 2048./diff(spikes);
            % if fr(1) < 1.5
            %     spikes(1) = [];
            % end
            % if fr(end) < 1.5
            %     spikes(end) = [];
            % end
        end

        %% Rule: adjust spikes
        function [autoEd, spikes, removed_spikes]= adjust_spikes(autoEd, ipt_t, spikes, fr_high_threshold)
            if autoEd.adjust_enable == 0
                removed_spikes = [];
                return;
            end

            % S1. get firing rate threshold and ipt threshold
            fsamp = 2048;
            fr = fsamp./diff(spikes);
            spikes(fr==inf) = [];
            fr = fsamp./diff(spikes);
            fr = [fr, fr(end)];
            if mean(fr) > 5 % remove fr lower than 5 Hz
                fr_g5 = fr(fr>=5);
            else
                fr_g5 = fr(fr>=2);
            end
            fr_mean = mean(fr_g5);
            fr_std = std(fr_g5);

            fr_threshold = min([fr_mean*1.5, fr_high_threshold]);

            ipt_spikes = ipt_t(spikes);
            ipt_mean = mean(ipt_spikes);

            % S2. initialize parameters and go through each spike
            spike_cnt = length(spikes);
            editing_hist = struct();
            remove_flag = spikes*0;
            removed_spikes = [];
            spike_cnt = length(spikes);
            count = zeros(1,3);
            adjusted_cnt = 0;
            seg_len = 10;
            seg_len_half = seg_len/2;
            for k = 2:spike_cnt-2
                % basic info
                seg_si = max([1, k-seg_len_half]);
                seg_ei = min([k+seg_len_half, spike_cnt-1]);
                fr_seg = fr(seg_si:seg_ei);
                fr_seg(fr_seg<3) = [];
                fr_seg(fr_seg>50) = [];
                fr_seg_median = median(fr_seg);

                spike_pre = spikes(k-1);
                spike_cur = spikes(k);
                spike_post = spikes(k+1);
                spike_post2 = spikes(k+2);

                fr_pre = fr(k-1);
                fr_cur = fr(k);
                fr_post = fr(k+1);
                fr_post2 = fr(k+2);

                if fr_cur > 100 || fr_pre > 100
                    continue;
                end

                % shift
                % %                 hlfr_flag = fr_cur>fr_seg_median*1.2 && fr_post<fr_seg_median*0.8;
                %                 hlfr_flag = fr_cur > fr_post*1.5;
                % %                 lhfr_flag = fr_cur<fr_seg_median*0.8 && fr_post>fr_seg_median*1.2;
                %                 lhfr_flag = fr_cur*1.5 < fr_post;
                %                 hlfr_flag = fr_cur > fr_post*1.4 && fr_post > 3;
                %                 %                 lhfr_flag = fr_cur<fr_seg_median*0.8 && fr_post>fr_seg_median*1.2;
                %                 lhfr_flag = fr_cur*1.4 < fr_post && fr_cur > 3;
                hlfr_flag = fr_cur > fr_post*1.5 && fr_post > 3;
                lhfr_flag = fr_cur*1.5 < fr_post && fr_cur > 3;
                sfr_flag = hlfr_flag + lhfr_flag;
                fr_accept = 0;
                ipt_accept = 0;

                if sfr_flag
                    if lhfr_flag
                        %                         potential_spike = find_potential_spike(ipt_t, spike_cur, spike_post, 1);
                        % find potential spike
                        potential_spike = find_potential_spike_loc(ipt_t, spike_cur, spike_post, 0.75);
                        if abs(spike_post - potential_spike) > 20
                            % calculate potential firing rate
                            fr_potential_pre = fsamp./(potential_spike-spike_cur);
                            fr_potential_post = fsamp./(spike_post2-potential_spike);
                            fr_pre_post_ratio = fr_potential_pre./fr_potential_post;
                            %                             fr_pre_post_ratio_old = fr_cur./fr_post;
                            fr_accept = fr_pre_post_ratio>0.75 & fr_pre_post_ratio<1.33;
                            fr_pre_post_cov = std([fr_potential_pre, fr_potential_post])/mean([fr_potential_pre, fr_potential_post]);
                            fr_pre_post_cov_old = std([fr_cur, fr_post])/mean([fr_cur, fr_post]);
                            if fr_accept == 0 && fr_pre_post_cov_old > fr_pre_post_cov && fr_pre_post_cov < 0.3 % && fr_pre_post_cov_old > 0.3
                                potential_spike_display(ipt_t, spikes, [spike_cur, potential_spike, spike_post], 10, 'remove&add one');
                                fr_accept = 2;
                            end
                        else
                            fr_accept = 0;
                        end
                    else
                        %                         potential_spike = find_potential_spike(ipt_t, spike_post, spike_post2, 1);
                        potential_spike = find_potential_spike_loc(ipt_t, spike_post, spike_post2, 0.25);
                        if abs(spike_post - potential_spike) > 20
                            fr_potential_pre = fsamp./(potential_spike-spike_cur);
                            fr_potential_post = fsamp./(spike_post2-potential_spike);
                            fr_pre_post_ratio = fr_potential_pre./fr_potential_post;
                            %                             fr_pre_post_ratio_old = fr_cur./fr_post;
                            fr_accept = fr_pre_post_ratio>0.75 & fr_pre_post_ratio<1.33;
                            fr_pre_post_cov = std([fr_potential_pre, fr_potential_post])/mean([fr_potential_pre, fr_potential_post]);
                            fr_pre_post_cov_old = std([fr_cur, fr_post])/mean([fr_cur, fr_post]);
                            if fr_accept == 0 && fr_pre_post_cov_old > fr_pre_post_cov && fr_pre_post_cov < 0.3 % && fr_pre_post_cov_old > 0.3
                                potential_spike_display(ipt_t, spikes, potential_spike, 10, 'remove&add one');
                                fr_accept = 2;
                            end
                        else
                            fr_accept = 0;
                        end
                    end

                    ipt_old = ipt_t(spike_post);
                    spike_old = spike_post;
                    ipt_threshold_l = min(ipt_t([spike_cur, spike_post]))*autoEd.ipt_ratio_a_l;
                    ipt_threshold_g = ipt_mean*autoEd.ipt_ratio_a_g;

                    % if ipt_check is disabled, release the ipt requirement
                    %                     if autoEd.ipt_check
                    %                         ipt_threshold = max([ipt_threshold_l, ipt_threshold_g]);
                    %                     else
                    %                         ipt_threshold = min([ipt_threshold_l, ipt_threshold_g]);
                    %                     end
                    ipt_threshold = min([ipt_threshold_l, ipt_threshold_g]);

                    ipt_accept = ipt_t(potential_spike) >= ipt_threshold || ipt_t(potential_spike) >= ipt_old*0.5;

                    if fr_accept && ipt_accept
                        fig = potential_spike_display(ipt_t, spikes, [spike_cur, potential_spike], 1, 'remove&add one');
                        if autoEd.debug_flag
                            saveas(fig, ['hlfr_adjusted_', num2str(autoEd.low_fr_processed), autoEd.fig_format]);
                            autoEd.low_fr_processed = autoEd.low_fr_processed + 1;
                        end
                        spikes(k+1) = potential_spike;
                        fr(k) = fr_potential_pre;
                        fr(k+1) = fr_potential_post;
                        adjusted_cnt = adjusted_cnt + 1;
                        if fr_accept == 2
                            autoEd.disp_info('Adjusted to improve FR!');
                        end
                    elseif ipt_accept
                        fig = potential_spike_display(ipt_t, spikes, potential_spike, 0, 'adjust high fr');
                        if autoEd.debug_flag && autoEd.save_ignored_flag
                            saveas(fig, ['hfr_ignored', num2str(autoEd.low_fr_ignored), autoEd.fig_format]);
                            autoEd.low_fr_ignored = autoEd.low_fr_ignored + 1;
                        end
                    elseif fr_accept
                        fig = potential_spike_display(ipt_t, spikes, potential_spike, 0, 'adjust high fr');
                        if autoEd.debug_flag && autoEd.save_ignored_flag
                            saveas(fig, ['hfr_ignored', num2str(autoEd.low_fr_ignored), autoEd.fig_format]);
                            autoEd.low_fr_ignored = autoEd.low_fr_ignored + 1;
                        end
                    end
                end

                editing_hist(k).spike = spikes(k);
                editing_hist(k).fr = fr(k);
                editing_hist(k).hlfr_flag = hlfr_flag;
                editing_hist(k).lhfr_flag  = lhfr_flag;
                editing_hist(k).sfr_flag  = sfr_flag;
                editing_hist(k).fr_accept  = fr_accept;
                editing_hist(k).ipt_accept = ipt_accept;
            end

            % update spikes
            autoEd.flag_hist = editing_hist;
            if adjusted_cnt
                autoEd.disp_info(['Adjusted to improve FR: ', num2str(adjusted_cnt)]);
            end
            %             autoEd.adjusted_spikes = [autoEd.adjusted_spikes, ]
        end

        %% Rule: refine spikes
        function [autoEd, updated_spikes, updated_cnt] = refine_spikes(autoEd, ipt_t, spikes, margin_s)
            % fine tune the spike within margin_s
            if nargin == 3
                margin_s = 5;
            end
            updated_spikes = spikes;
            tuned_spikes_t = spikes*0;
            data_len = length(ipt_t);
            updated_cnt = 0;

            % go through each spike
            for si = 1:length(updated_spikes)
                spike_cur = updated_spikes(si);
                ipt_cur = ipt_t(spike_cur);
                try
                    if spike_cur <= margin_s || spike_cur>data_len-margin_s
                        continue;
                    end
                catch
                    autoEd.disp_info('e')
                end

                % fine tune with margin_s
                ipt_seg = ipt_t(spike_cur-margin_s:spike_cur+margin_s);
                [spike_ipt, spike_i] = max(ipt_seg);
                potential_spike = spike_i+spike_cur-margin_s-1;

                % update spike if potential spike has higher ipt
                if (spike_ipt>ipt_cur) %&& abs(potential_spike - spike_cur)<40
                    fig = potential_spike_display(ipt_t, updated_spikes, [spike_cur, potential_spike], 3, 'refine spike');
                    updated_spikes(si) = potential_spike;
                    updated_cnt = updated_cnt + 1;
                    tuned_spikes_t(si) = 1;
                elseif spike_ipt>ipt_cur
                    fig = potential_spike_display(ipt_t, updated_spikes, [spike_cur, potential_spike], 3, 'refine spike');
                    updated_spikes(si) = potential_spike;
                    updated_cnt = updated_cnt + 1;
                    tuned_spikes_t(si) = 1;
                end
            end

            % update spikes
            autoEd.tuned_spikes = spikes(tuned_spikes_t==1);
            updated_spikes(diff(updated_spikes)==0) = [];
            if updated_cnt
                autoEd.disp_info(['Spike tuned: ', num2str(updated_cnt)]);
            end
        end

        function [autoEd, updated_spikes, updated_cnt] = refine_spikes_ipt(autoEd, ipt_t, spikes)
            % update spike if there is nearby high IPT
            % the nearby spike should be within 33% of the interval
            fsamp = 2048;
            updated_spikes = spikes;
            tuned_spikes_t = spikes*0;
            updated_cnt = 0;
            %             ipt_hist = zeros(length(spikes), 5);

            for si = 1:length(updated_spikes) - 1
                spike_cur = updated_spikes(si);
                spike_post = updated_spikes(si+1);
                interval_max = (spike_post - spike_cur)*0.33;
                ipt_cur = ipt_t(spike_cur);
                ipt_post = ipt_t(spike_post);

                % fine tune with margin_s
                ipt_seg = ipt_t(spike_cur+10:spike_post-10);
                if isempty(ipt_seg)
                    continue;
                end
                [spike_ipt, spike_i] = max(ipt_seg);
                potential_spike = spike_i+spike_cur+10-1;

                % update spike if potential spike has higher ipt
                if spike_ipt > ipt_cur && abs(potential_spike - spike_cur) < interval_max
                    fig = potential_spike_display(ipt_t, updated_spikes, [potential_spike, spike_cur], 11, 'refine spike');
                    updated_spikes(si) = potential_spike;
                    updated_cnt = updated_cnt + 1;
                    tuned_spikes_t(si) = 1;
                elseif spike_ipt > ipt_post && abs(potential_spike - spike_post) < interval_max
                    fig = potential_spike_display(ipt_t, updated_spikes, [potential_spike, spike_post], 11, 'refine spike');
                    updated_spikes(si+1) = potential_spike;
                    updated_cnt = updated_cnt + 1;
                    tuned_spikes_t(si+1) = 1;
                end
                %                 ipt_hist(si, :) = [spike_cur, spike_ipt, ipt_cur, ipt_post, interval_max, abs(potential_spike - spike_post)];
            end

            autoEd.tuned_spikes = spikes(tuned_spikes_t==1);
            updated_spikes(diff(updated_spikes)==0) = [];
            if updated_cnt
                autoEd.disp_info(['Spike refined (high IPT): ', num2str(updated_cnt)]);
            end
        end

        %% other functions
        function autoEd = disp_info(autoEd, msg)
            if autoEd.disp_flag
                disp(msg);
                autoEd.write_csv(msg);
            end
        end

        function autoEd = write_csv(autoEd, msg)
            % disp(msg);
            % write_csv(msg);
        end

        function [autoEd, ipt_t, Wts] = update_ipts(autoEd, spikes, wSIG_m)
            if length(spikes)<=5
                ipt_t = wSIG_m(1, :)*0;
                Wts = wSIG_m(:, 1)*0;
                return;
            end
            Wt = sum(wSIG_m(:,spikes),2);
            icasig = Wt'*wSIG_m;
            icasig = abs(icasig).*icasig; % using square to increase the difference between spikes and non-spikes

            % get normalization factor
            spikes_icasig = icasig(spikes);
            [spikes_icasig_sorted, index] = sort(spikes_icasig, 'descend');
            if spikes_icasig_sorted(1)/spikes_icasig_sorted(2) > 5
                spike_t = spikes(index(1));
                icasig(spike_t-2:spike_t+2) = 0;
                spikes_icasig_sorted(1) = 0;
            elseif spikes_icasig_sorted(1)/spikes_icasig_sorted(2) > 3
                spike_t = spikes(index(1));
                if spike_t > 2 && spike_t < length(icasig)-2
                    icasig(spike_t-2:spike_t+2) = icasig(spike_t-2:spike_t+2)*spikes_icasig_sorted(2)/spikes_icasig_sorted(1);
                else
                    icasig(spike_t) = icasig(spike_t)*spikes_icasig_sorted(2)/spikes_icasig_sorted(1);
                end
            end
            icasig_max = mean(spikes_icasig_sorted) + std(spikes_icasig_sorted);

            % normalize icasig
            ipt_t = icasig./icasig_max;
            Wts = Wt;
        end

        function [autoEd, fig] = plot_spikes_ipt_fr(autoEd, ipt_t, spikes)
            fsamp = 2048;
            fr = [fsamp./diff(spikes), 0];
            ti = 1:length(ipt_t);
            % ti = ti./fsamp;

            fig = figure();
            % plot ipt and spikes
            subplot(2, 1, 1);
            plot(ti, ipt_t);
            hold on;
            plot(ti(spikes), ipt_t(spikes), 'ro');
            xlim([ti(1), ti(end)+1]);
            %             xlim([ti(spikes(1)- 100), ti(spikes(end)+100)]);

            % plot firing rate
            subplot(2, 1, 2);
            plot(ti(spikes), fr, 'bo');
            xlim([ti(1), ti(end)+1]);
            %             xlim([ti(spikes(1)- 100), ti(spikes(end)+100)]);
        end

        function [autoEd, fig] = plot_spike_firing(autoEd, ipt, spikes, spikes_bad)
            seg_index = find(diff(spikes)>2048*5);
            pause_start = [spikes(seg_index); spikes(seg_index+1)];
            spike_segments = [spikes(1), pause_start(:)', spikes(end)];
            seg_count = length(spike_segments)/2;
            spike_segments = reshape(spike_segments, 2, seg_count)';
            duration = spike_segments(:,2)- spike_segments(:,1);
            spike_segments(duration<1000, :) = [];
            seg_count = size(spike_segments, 1);
            fig = figure();
            for seg_i = 1:seg_count
                starti = spike_segments(seg_i, 1);
                endi = spike_segments(seg_i, 2);
                subplot(seg_count, 1, seg_i);
                plot(ipt);
                hold on; plot(spikes, ipt(spikes), 'bo');
                if nargin == 3
                    if sum(spikes==spikes_bad(1))
                        plot(spikes_bad, ipt(spikes_bad), 'rx');
                        plot(spikes_bad, ipt(spikes_bad), 'ro');
                    else
                        plot(spikes_bad, ipt(spikes_bad), 'ro');
                        plot(spikes_bad, ipt(spikes_bad), 'yd');
                    end
                end
                xlim([starti-1000, endi+1000]);
            end
        end

    end
end