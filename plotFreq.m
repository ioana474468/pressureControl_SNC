% SCPI/SIeC Laboratory
%
% frecventa(SIGNAL, FS, TYPE, plotAxisT, plotAxisF)
% 
% This script plots the frequency representation of a sampled with Fs (Hz). 
% The results are also saved into a PDF file: FrecventaImg.pdf !
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
% Major Rework: 10.03.2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotFreq(signal,Fs,type, plotAxisT, plotAxisF)
if nargin < 2
    error(message('MATLAB: plotFreq: At least first 2 inputs !'))
end
if nargin==2
    frecventa(signal,Fs);
end
if nargin==3
    frecventa(signal,Fs,type);
end
if nargin==4
    frecventa(signal,Fs,type,plotAxisT);
end
if nargin==5
    frecventa(signal,Fs,type, plotAxisT, plotAxisF);
end