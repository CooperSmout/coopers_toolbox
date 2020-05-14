

function plotFFT(y,Fs,lim)
        
    figure
    if nargin==3
        xlims = lim;
    else
        xlims = [];
    end    
    
    L = length(y);
    NFFT = 2^nextpow2(L);
    Y = fft(y,NFFT)/L;
    f = Fs/2*linspace(0,1,NFFT/2+1);
    plot(f,2*abs(Y(1:NFFT/2+1)))
    
    if ~isempty(xlims)
        xlim(xlims)
    end
    
    title('Single-Sided Amplitude Spectrum of y(t)')
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
    
    
end