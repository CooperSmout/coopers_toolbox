% Sends eeg trigger to Biopac - make sure Cogent is loaded or will
% experience delays on first call

function sendTrigger(ioObject, ioAddress, triggerNum, length)

io64(ioObject, ioAddress, triggerNum); % sends trigger number triggerNum
wait(length);
io64(ioObject, ioAddress, 0); % eeg trigger reset to 0; so the trigger numbers are not added together
