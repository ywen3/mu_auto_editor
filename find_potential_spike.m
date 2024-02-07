function potential_spikes = find_potential_spike(ipt, spike_cur, spike_post, tops)
if nargin == 1
    spike_cur = 1;
    spike_post = length(ipt);
    tops = 1;
elseif nargin == 2
    spike_cur = 1;
    spike_post = length(ipt);
    tops = 1;
elseif nargin == 3
	tops = 1; 
end

if length(ipt) < 20
    potential_spikes = spike_cur;
    return;
end

ipt_seg = ipt(spike_cur+10:spike_post-10);
% sort the ipt by magnitude
[v, potential_spikes] = sort(ipt_seg, 'descend');
if isempty(potential_spikes)
    disp('Bad!');
    potential_spikes = spike_cur;
    return;
end
if length(potential_spikes)>10
    potential_spikes = potential_spikes(1:10);
end

dlen = length(ipt_seg);
% remove very close peaks
j = 1;
while j < length(potential_spikes)-1
    s = j+1;
    while s<=length(potential_spikes)
        if abs(potential_spikes(j)-potential_spikes(s))<10
            potential_spikes(s) = [];
        else
            s = s + 1;
        end
    end
    j = j + 1;
end

% keep candidates with ipt higher than 60% of the max possible ipt
potential_spikes_ipts = ipt_seg(potential_spikes);
potential_ipt_max = potential_spikes_ipts(1);
potential_spikes = potential_spikes(potential_spikes_ipts>=potential_ipt_max*0.5);
potential_spikes_ipts = potential_spikes_ipts(potential_spikes_ipts>=potential_ipt_max*0.5);

% sort the candidates using firing rate
fr_pre = 2048./(potential_spikes+5);
fr_post = 2048./(dlen-potential_spikes+5);
fr_mean = 2048*2/dlen;
if fr_mean > 3
    fr_error = abs(fr_pre-fr_post)./(2*fr_mean);
    fr_error(fr_error<0.25) = 0;
    [fr_error, s_order] = sort(fr_error);
    potential_spikes = potential_spikes(s_order);
else
    ipt_error = (potential_ipt_max - potential_spikes_ipts);
    [ipt_error, s_order] = sort(ipt_error);
    potential_spikes = potential_spikes(s_order);
end
% 
% ipt_error = (potential_ipt_max - potential_spikes_ipts)./potential_ipt_max;
% [fr_error, s_order] = sort(fr_error+ipt_error);
if length(potential_spikes) >= tops
    potential_spikes = potential_spikes(1:tops);
end

potential_spikes = potential_spikes + spike_cur + 9;
if isempty(potential_spikes)
    potential_spikes = spike_cur;
end
end
