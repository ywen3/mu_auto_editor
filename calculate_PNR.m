function PNR = calculate_PNR(spikes, ipt, sampleRate)
% calculate PNR of the motor unit spike train
% Reference: "Accurate identification of motor unit discharge patterns from
% high-density surface EMG and validation with a novel signal-based performance metric"
%
% input: 
%       spikes: the index of spikes
%       ipt: innervated pulse train signal or the ICA source signal
% output:
%       PNR: pulse-to-noise-ratio
%
% Example:
%       spikes = MUPulses{1};
%       ipt = IPTs(1, :);
%       PNR = calculate_PNR(ipt, spikes)
%
% Author: Yue Wen at University of Central Florida, Nov. 20, 2023.

if nargin == 2
    sampleRate = 2048;
end

% the margin is to exclude the ipt around each spike
% the margin is 3 when sample rate is 2048 Hz
margin = round(sampleRate*3/2048);
marginIndex = -margin:margin;
spikeSurrounds = spikes + marginIndex';
signalIndex = zeros(size(ipt));
signalIndex(spikeSurrounds(:)) = 1;

nonSpikeIpt = ipt(signalIndex==0);
spikeIpt = ipt(spikes);
PNR = 10*log10(mean(spikeIpt.^2)/mean(nonSpikeIpt.^2));

% illustrate the spikes and non-spikes
% figure();
% plot(ipt);
% hold on;
% plot(spikes, ipt(spikes), 'ro');
% plot(find(signalIndex==0), ipt(signalIndex==0), 'bo');
end
