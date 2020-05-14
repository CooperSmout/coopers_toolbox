
% duration in seconds
% sampling rate in samples/second

function mat = puretone(Frequency,Duration,SamplingRate,varargin)
    if numel(varargin)
        volume = varargin{1};
    else
        volume=1;
    end
    mat = sin((1:Duration*SamplingRate)*2*pi*Frequency/SamplingRate)*volume;
return