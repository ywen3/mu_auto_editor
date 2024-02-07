function fig = potential_spike_display(ipt, spikes, potential_spike, type_flag, caption)
% potential add (high ipt, low fr, etc), 
% 0 for ignore; 1 for add; 2 for remove; 3 for fine tune
% remove (low ipt, high fr, etc),
% adjustment
% 
if nargin == 4
   caption = ''; 
end
global plot_flag;
global mu_tag;
% 
% 
fsamp = 2048;
interval = diff(spikes);
fr = [fsamp./interval, 0];
xrange = 1.5;
xt = 1:length(ipt);
x_scale = fsamp;
% 
% spike_info.type = type_flag;
% potential_index = find(max(spikes-potential_spike(1)));
% ipt(potential_spike);
% potential_index = find(spikes==potential_spike(1));

% plot_flag = 0;
fig = 0;
if plot_flag %|| type_flag == 11 %|| contains(caption, 'Double') %|| contains(caption, 'Low ipt')
    fig = figure(10);
    movegui('east');
    % plot firing rate
    subplot(2, 1, 1);
    plot(spikes/x_scale, fr, 'bo');
    hold on;
    potential_index = find(spikes==potential_spike(1));
    if ~isempty(potential_index)
        plot(potential_spike/x_scale, fr(potential_index), 'ro');
    else
        [~, ti]= min(abs(spikes-potential_spike(1)));
        potential_spike_c = spikes(ti);
        plot(potential_spike_c/x_scale, fr(ti), 'ro');
    end
    xlim([mean(potential_spike)-2048*xrange, mean(potential_spike)+2048*xrange]/x_scale);
    ylabel('Firing rate');
    hold off;
    % plot ipt
    subplot(2, 1, 2);
    plot(xt/x_scale, ipt);
    hold on;
    plot(spikes/x_scale, ipt(spikes), 'bo')
    xlim([mean(potential_spike)-1000, mean(potential_spike)+1000]/x_scale);
    if type_flag == 0 || type_flag == 10
        plot(potential_spike/x_scale, ipt(potential_spike), 'rx');
        title(['ignore: ', caption]);
    elseif type_flag == 1
        plot(potential_spike/x_scale, ipt(potential_spike), 'rd');
        title(['add: ', caption]);
    elseif type_flag == 2
        plot(potential_spike/x_scale, ipt(potential_spike), 'rx');
        title(['remove: ', caption]);
        if ipt(potential_spike)>0.4
            disp('high ipt');
        end
    elseif type_flag == 3
        plot(potential_spike(1)/x_scale, ipt(potential_spike(1)), 'rd');
        plot(potential_spike(2)/x_scale, ipt(potential_spike(2)), 'rx');
        title(['fine tune: ', caption]);
    elseif type_flag == 5
        plot(potential_spike/x_scale, ipt(potential_spike), 'ro');
        title(['debug: ', caption]);
    elseif type_flag == 11 || type_flag == 12
        plot(potential_spike/x_scale, ipt(potential_spike), 'ro');
        title(['debug: ', caption]);
    end
%     xlim([mean(potential_spike)-2000, mean(potential_spike)+2000]);
    ylabel('IPT height');
    hold off;
    xlim([mean(potential_spike)-2048*xrange, mean(potential_spike)+2048*xrange]/x_scale);
    
%     if type_flag == 3
%         pause(0.001);
%     end

    if type_flag == 11
        folder = 'C:\Users\ywen\Documents\Auto editing\auto edit data\debug\';
%         folder = 'F:\Auto editting\online data\debug\';
        if exist(folder, 'dir')
           saveas(fig, [folder, mu_tag, '_debug_', caption, '_', num2str(round(potential_spike(1)/10)), '.png']); % datestr(now,'HH_MM_SS.FFF'),
        end
    end

    if type_flag == 1 || type_flag == 2
        % folder = 'C:\Users\ywen\Documents\Auto_editing_pp\vi\';
        root_folder = go_root_folder(1);
        folder = [root_folder, 'vi\'];
        if exist(folder, 'dir')
            if type_flag == 1
                saveas(fig, [folder, mu_tag, '_add_', caption, '_', num2str(round(potential_spike/10)), '.fig']); % datestr(now,'HH_MM_SS.FFF'),
            else
                saveas(fig, [folder, mu_tag, '_del_', caption, '_', num2str(round(potential_spike/10)), '.fig']); % datestr(now,'HH_MM_SS.FFF'),
            end
        end
    end
end
end