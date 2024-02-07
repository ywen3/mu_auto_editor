function [IPTs, MUPulses] = flip_ipts_spikes(IPTs, MUPulses)
if iscell(MUPulses)
    for mu_index = 1:length(MUPulses)
        ipt = IPTs(mu_index, :);
        spikes = MUPulses{mu_index};
        % flip signals
        ipt = flip(ipt);
        spikes = length(ipt) - spikes + 1;
        spikes = sort(spikes);
        % update IPTs and spikes
        IPTs(mu_index, :) = ipt;
        MUPulses{mu_index} = spikes;
    end
else
    % flip signals
    IPTs = flip(IPTs);
    MUPulses = length(IPTs) - MUPulses + 1;
    MUPulses = sort(MUPulses);
end
end