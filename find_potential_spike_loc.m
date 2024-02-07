function potential_spikes = find_potential_spike_loc(ipt, spike_cur, spike_post, loc)
if nargin == 1
    spike_cur = 1;
    spike_post = length(ipt);
elseif nargin == 3
	loc = 0.5;
end

ipt_seg = ipt(spike_cur+10:spike_post-10);
% sort the ipt by magnitude
[v, potential_spikes] = sort(ipt_seg, 'descend');
if isempty(potential_spikes)
    disp('Bad!');
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
% sort the candidates using firing rate
spike_locations = potential_spikes/dlen;
location_errors = abs(spike_locations-loc);

[fr_error, s_order] = sort(location_errors);
potential_spikes = potential_spikes(s_order);
if ~isempty(potential_spikes)
    potential_spikes = potential_spikes(1);
    potential_spikes = potential_spikes + spike_cur + 9;
else
    potential_spikes = spike_cur;
end
end
