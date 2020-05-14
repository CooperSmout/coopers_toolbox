

config_display(1, 2, [0 0 0], [1 1 1], 'Helvetica', 25, 4 )
config_keyboard
start_cogent

t1 = time;

tic
for t = 1:200
    waituntil(t1+t*5);
    times(t)= time;
    times2(t) = toc;
end

stop_cogent



times'