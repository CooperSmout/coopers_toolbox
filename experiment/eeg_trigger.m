function eeg_trigger(io_obj, io_address, trigger, eeg_trigger_length)

io32(io_obj, io_address, trigger); % trigger to send = trigger (trigger number)
wait(eeg_trigger_length);
io32(io_obj, io_address, 0); % eeg trigger reset to 0; so the trigger numbers are not added together