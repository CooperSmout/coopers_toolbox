function trig = readTrigger(ioObject, ioAddress, portNum)

trig = io64(ioObject, ioAddress) ~= portNum; % receives trigger and subtracts 
