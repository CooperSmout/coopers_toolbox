function sendTriggerPTB(ioObject, ioAddress, triggerNum, length)

if IsWindows
    io64(ioObject, ioAddress, triggerNum); % sends trigger number triggerNum
    WaitSecs('UntilTime', GetSecs+length/1000);
    % wait(length);
    io64(ioObject, ioAddress, 0); % eeg trigger reset to 0; so the trigger numbers are not added together
elseif isempty(ioObject)
    fprintf([num2str(triggerNum) ' '])
end
