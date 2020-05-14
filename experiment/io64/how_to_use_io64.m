
% VARIABLES YOU NEED TO DEFINE BEFORE WE BEGIN ----------------------------

trigger_value = 1;
% This is the value that will be sent to the parallel port. It can be any number up to 255.

pulse_duration = 2;
% This will be how long the trigger is sent for, in milliseconds.

port_address_1 = hex2dec('D050'); % this is the address of the parallel port on the computer in Psyc
%port_address_1 = hex2dec('378'); % this is the address of the parallel port on the computer in QBI
% This converts the hexadecimal string which identifies your port to a number, that io64 will use as the address of your parallel port.
% To find out the hexadecimal string of your port:
%   Go to your computer's "Control Panel".
%   Go to "Hardware" (sometimes called "Hardware and Sound" or something like that).
%   Go to "Device Manager".
%   Click the + sign next to "Ports (COM & LPT)".
%   Double-click on the port you want to find the address of.
%   Click on the "Resources" tab.
%   The hexadecimal string of that port will be the first bit (before the "-" symbol) in the "Setting" column next to "I/O Range".
% If you want to talk to more than one port, just define more than one port address.

send_triggers = 1; % if you set this variable to "0", the program will skip
                   % over the io lines of code when you are programming or
                   % piloting on a different computer. Otherwise, the
                   % program will encounter an error when it tries to talk
                   % to a non-existent parallell port.


% CODE TO USE IO64 --------------------------------------------------------

% Before you can send or read information from a parallel port, you have to put in these lines of code
if send_triggers
    ioObj = io64; % create an instance of the io64 object
    status = io64(ioObj); % initialize the inpout64 system driver
end


% At the start of your experiment, you should also set any ports you are using to have a value of zero
if send_triggers
    io64(ioObj,port_address_1,0); % set the port called port_address_1 to 0
end


% This is how you send information to a parallel port
if send_triggers
    io64(ioObj,port_address_1,trigger_value); % set the port called port_address_1 to your desired trigger value
    wait(pulse_duration); % wait a while
    io64(ioObj,port_address_1,0); % set the port called port_address_1 back to zero
end


% This is how you read information from a parallel port
if send_triggers
    port_value = io64(ioObj,port_address_1); % port_value will be the current value of the port called port_address_1
    % (in this case it will be zero because we just set it to zero)
end