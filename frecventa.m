% SCPI/SIeC Laboratory
%
% frecventa(SIGNAL, FS, TYPE, plotAxisT, plotAxisF)
% 
% This script plots the frequency representation of a sampled with Fs (Hz). 
%
% SIGNAL = is any vector containing a time domain signal
% FS = is the sampled frequency of the signal given (in Hz)
% [TYPE] = 's' standard view, 'p' for periodogram
% [plotAxisT]= [XMIN XMAX YMIN YMAX] axes limits for the time domain plot
% [plotAxisF]= [XMIN XMAX YMIN YMAX] axes limits for the Freq domain plot
% Legend [...] - optional argument
%
% Author: Valentin Tanasa
% Date: 8.10.2010
% Last Update: 8.03.2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=frecventa(signal,Fs,type, plotAxisT, plotAxisF)
%%%%%%  INPUT DATA  %%%%%%
number_of_samples = length(signal);
final_time=(number_of_samples-1)/Fs;
signal=signal(:);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 1. Plot continuous and sampled signal in time domain %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Define continuous time interval
    t = linspace(0, final_time, number_of_samples);
     %Take fourier transform - four 'continuous' "signal"
    fftSignal = fft(signal);
    
    %apply fftshift to put it in the form we are used to
    fftSignal = fftshift(fftSignal);
    
    %Next, calculate the frequency axis, which is defined by the sampling rate
    T = t(2)-t(1);
    fs = 1/T;
    f = fs/2*linspace(-1,1,number_of_samples);
if (nargin<=2)||((nargin>=3)&&(type~='p')) 
    figure;
    %Plot the "signal" - in time domain
    subplot(211),plot(t, signal);hold on;
    title(['Time-Domain signal, Fs=',num2str(Fs),' Hz']);
    if (nargin>=4)&&(length(plotAxisT)==4)
        axis(plotAxisT);
    end
    ylabel('Amplitude');
    xlabel('Time (s)');
       
      
    %Since the signal is complex, we need to plot the magnitude to get it to
    %look right, so we use abs (absolute value)
    subplot(212), plot(f, ((abs(fftSignal))));
    title('magnitude FFT of signal"');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude ');
    if (nargin==5)&&(length(plotAxisF)==4)
        axis(plotAxisF);
    end
    
    %legend('Continuous')
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 3. Plot sampled signal in frequency domain %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure;
    subplot(211),plot(t, signal);hold on;
    title(['Time-Domain signal, Fs=',num2str(Fs),' Hz']);
    ylabel('Amplitude');
    xlabel('Time (s)');
    if (nargin>=4)&&(length(plotAxisT)==4)
        axis(plotAxisT);
    end
    %subplot(312),plot(2*f/fs, abs(dfftSignal),'b');
    subplot(212), periodogram(signal,[],'onesided',[],Fs);
    if (nargin==5)&&(length(plotAxisF)==4)
         axis((plotAxisF));
    end
end
V=axis;
text(V(1),V(3)-((V(4)-V(3))/(4)),datestr(now))
set(gcf, 'PaperPosition', [0 0 8 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [8 5]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'FrecventaImg', 'pdf') %Save figure
%
% xc = xcorr(dsignal,'biased');
% xcft = fft(xc);
% ff = fs*linspace(0,1,2*dnumber_of_samples-1);
% subplot(313), plot(ff,abs(xcft));
% title('magnitude FFT of "dsignal"');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude ');
%
