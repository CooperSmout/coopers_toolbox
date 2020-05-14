
% VARIABLES YOU NEED TO DEFINE BEFORE WE BEGIN ----------------------------

trigger_value = 1;
% This is the value that will be sent to the parallel port. It can be any number up to 255.

pulse_duration = 2;
% This will be how long the trigger is sent for, in milliseconds.

port_address_1 = hex2dec('378');
% This converts the hexadecimal string which identifies your port (in this case '378') to a number, that io32 will use as the address of your parallel port.
% To find out the hexadecimal string of your port:
%   Go to your computer's "Control Panel".
%   Go to "Hardware" (sometimes called "Hardware and Sound" or something like that).
%   Go to "Device Manager".
%   Click the + sign next to "Ports (COM & LPT)".
%   Double-click on the port you want to find the address of.
%   Click on the "Resources" tab.
%   The hexadecimal string of that port will be the first bit (before the "-" symbol) in the "Setting" column next to "I/O Range".
%   The hexadecimal strings for the first three parallel ports will normally be: '378', 'CCC0', and 'CCD0'.
% If you want to talk to more than one port, just define more than one port address.


% CODE TO USE IO32 --------------------------------------------------------

% Before you can send or read information from a parallel port, you have to put in these 2 lines of code
ioObj = io32; % create an instance of the io32 object
status = io32(ioObj); % initialize the inpout32 system driver


% At the start of your experiment, you should also set any ports you are using to have a value of zero
io32(ioObj,port_address_1,0); % set the port called port_address_1 to 0


% This is how you send information to a parallel port
io32(ioObj,port_address_1,trigger_value); % set the port called port_address_1 to your desired trigger value
wait(pulse_duration); % wait a while
io32(ioObj,port_address_1,0); % set the port called port_address_1 back to zero


% This is how you read information from a parallel port
port_value = io32(ioObj,port_address_1); % port_value will be the current value of the port called port_address_1 
                                            % (in this case it will be zero because we just set it to zero)

