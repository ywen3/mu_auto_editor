function [added_spikes, removed_spikes, fine_tuned] = compare_mupulses(spikes_old, spikes_new)
% Compare two pulse trains and identify the added spikes and removed spikes
% input:
%   spikes_old: the old spikes or pulses
%   spikes_new: the new spikes or pulses
%
% output:
%   added_spikes: added spikes to the old spikes
%   removed_spikes: removed spikes from the old spikes
%   
% Author: Yue Wen at Shirley Ryan AbilityLab, 2022.
% 

% construct pulse train with spike index
if max(spikes_old)>1
    % input is spike indexes
    dlen = max([spikes_old(:); spikes_new(:)])+100;
    pulses_old = zeros(1, dlen);
    pulses_new = pulses_old;
    pulses_old(spikes_old) = 1;
    pulses_new(spikes_new) = 1;
else
    % input is pulse train
    pulses_old = spikes_old;
    pulses_new = spikes_new;
end

% extend pulse train for iterative comparison
added_spikes = [];
removed_spikes = [];
fine_tuned = [];
pulses_old = [zeros(1, 10), pulses_old, zeros(1, 10)];
pulses_new = [zeros(1, 10), pulses_new, zeros(1, 10)];
spikes_new_b = pulses_new;

tolerance = 2;
dlen = min([length(pulses_new)-tolerance, length(pulses_old)-tolerance]);
k = tolerance+1;
% check removed spikes and fine-tuned spikes
while k <= dlen
    if pulses_old(k) == 1 && pulses_new(k) == 0
        pulse_cnt = sum(pulses_new(k-tolerance:k+tolerance));
        if pulse_cnt == 1
            fine_tuned = [fine_tuned, k];
            pulses_new(k-tolerance:k+tolerance) = 0;
            pulses_new(k) = 1;
            k = k + tolerance;
        elseif pulse_cnt > 1
            fine_tuned = [fine_tuned, k];
            spike_indice = find(pulses_new(k-tolerance:k+tolerance)==1);
            spike_indice = spike_indice + k-tolerance - 1;
            [v, i]= min(abs(spike_indice - k));
            pulses_new(spike_indice(i)) = 0;
            pulses_new(k) = 1;
            % k = k + tolerance;
            % added_spikes = [added_spikes, ]
        elseif pulse_cnt == 0
            removed_spikes = [removed_spikes, k];
        end
    end
    k = k + 1;
end

% check added spikes
for k = tolerance+1:dlen
    if pulses_old(k) == 0 && pulses_new(k) == 1
        added_spikes = [added_spikes, k];
    end
end

% adjust for the initial 10 zeros
added_spikes = added_spikes - 10;
removed_spikes = removed_spikes - 10;

old_cnt = sum(pulses_old==1);
new_cnt = sum(spikes_new_b==1);
add_cnt = length(added_spikes);
removed_cnt = length(removed_spikes);
if old_cnt+add_cnt-removed_cnt ~= new_cnt
    disp('Compare_mupulses: the margin might be too large.');
end
end
