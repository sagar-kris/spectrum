function [pks, locs] = tuner()

clf reset;
close all;
clear all;
%clearvars -EXCEPT recObj myRecording
clc;

% % % Record your voice or instrument.
recordingTime = 3; %(seconds)
sampleRate=20000; % samples per second (make even); must be at least twice the largest frequency of interest
recObj = audiorecorder(sampleRate,8,1); % 8 bits, 1 channel
disp('Play a tone.')
recordblocking(recObj, recordingTime);
disp('End of recording.');

% Play back the recording.
play(recObj);

% Store data in double-precision array.
myRecording = getaudiodata(recObj);

% ---------------------------------------------------
T=1/sampleRate;
N=round(recordingTime/T); if 2*round(N/2)~=round(N), N=N+1; end
t=[0:1:N-1]*T; recordingTime=max(t);
subplot(2,2,1); plot(t,myRecording); xlabel('time [sec]');
% ---------------------------------------------------
temp=fft(myRecording);
temp(1:(10*N/sampleRate))=0; % remove frequencies below 10Hz
F=fftshift(temp);
freq=[-N/2:1:N/2-1]*sampleRate/N;   %map to frequencies in the fft
% ---------------------------------------------------
maxAmp=max(abs(F)); % rescale the frequency amps for convenience
subplot(2,2,3); plot(freq,abs(F)/maxAmp); xlabel('frequency [Hz]'); ylabel('relative amplitude');
axis(subplot(2,2,3),[0 inf -inf inf])
% ---------------------------------------------------


% FIND PEAKS
A=abs(F)/maxAmp; % A is the scaled version of fft_shift (max amp == 1)
A(1:[length(A)/2])=0; % remove negative frequencies - A goes from 0-10k Hz?
%first find rough est peak
fEstMax = sampleRate/40;%temp val, increase later
fEstDelta = 10; % must be a multiple of sampleRate/N; affects speed, jmg?
fTemp = zeros(fEstMax/fEstDelta,1);
value = zeros(fEstMax/fEstDelta,1);
nPeaksToCheck = 10;
for fEst = 10:fEstDelta:fEstMax
    i = fEst/fEstDelta;
    fTemp(i) = fEst;
    [pksTemp,locsTemp] = findpeaks(A,'MinPeakDistance',round(fEst*.9*N/sampleRate)); % added round, jmg?
    %NOTE: code would be much more stable if MaxPeakDistance of 1.1*fEst
    %was also set as requirement, but findpeaks() does not have this functionality.
    if length(pksTemp)<nPeaksToCheck
        value(i) = sum(pksTemp);
    else
        value(i) = sum(pksTemp(1:nPeaksToCheck));
    end
end
% at this point, fTemp is 10:20:...
% value is an array of the running sum of pksTemp from
% 1:min(length(pksTemp),nPeaksToCheck)
[dummy, index] = max(value);
fBestEst = fTemp(index); %initial rough est of fundamental frequency
threshold=mean(A)+3*std(A); %test with different levels of noise
%threshold=.01;
k = find(A<threshold); %find indicies of values below threshold
tempA = A; tempA(k) = 0;
% jnd=10*N/sampleRate; %peaks separated by less than just-noticable-difference of 10Hz will sound as one
% [pks,locs] = findpeaks(tempA,'MinPeakDistance',jnd);
[pks,locs] = findpeaks(tempA,'MinPeakDistance',round(fBestEst*.9*N/sampleRate)); % added round, jmg?

% DETERMINE FUNDAMENTAL FREQUENCY FROM MAJOR PEAKS
%arbitrarily use any peaks greater than X% of max peak to make a weighted
%average estimate of the perceived fundamental frequency
pksUsed=find(pks>.5);
mult=round(freq(locs)/fBestEst);
fundFreq = freq(round(sum((pks(pksUsed) .* (locs(pksUsed)-(N/2+1)) ./ mult(pksUsed)'))/ sum(pks(pksUsed)))+(N/2+1));


% CALCULATE AND OUTPUT THE FINAL SOLUTIONS
%fundFreq = 441; %debug
[refFreq, notes] = xlsread('note frequencies.xlsx','B1:M10'); %imports note frequency chart
position = 0:8;
[dummy, iMins] = min(abs(refFreq-fundFreq));
[error, noteIndex] = min(min(abs(refFreq-fundFreq)));
positionIndex=iMins(noteIndex);
targetFreq = refFreq(positionIndex,noteIndex); %return closest note frequency
centsError = 1200/log(2)*log(fundFreq/targetFreq);
disp(['fundamental frequency = ' num2str(round(fundFreq)) ' Hz'])
if centsError>0
    flatOrSharp = 'sharp';
else
    flatOrSharp = 'flat';
end
    
disp([flatOrSharp,' of ',notes{noteIndex},num2str(position(positionIndex)),' by ',num2str(round(abs(centsError))),' cents'])

% ADD OVERTONE FREQUENCIES TO PLOT
subplot(2,2,3);
%plot(freq,A,freq(locs),pks,'or',freq,threshold*ones(length(freq),1))
plot(freq,A,freq(locs),pks,'or')
xlabel('frequency [Hz]'); ylabel('scaled');
%axis(subplot(2,1,2),[0 inf -inf inf])
axis(subplot(2,2,3),[0 (freq(locs(1))*10.5) -inf inf]) 
for i = 1:length(locs)
    [dummy, iMins] = min(abs(refFreq-freq(locs(i))));
    [error, noteIndex] = min(min(abs(refFreq-freq(locs(i)))));
    positionIndex=iMins(noteIndex);
    label(i) = strcat(notes(noteIndex),num2str(position(positionIndex)));
end

% print output of findInstrument.m
likelyInstr = findInstrument(pks);
disp(['The instrument is most likely a ', likelyInstr]);

text(freq(locs),pks,label,'HorizontalAlignment','center','VerticalAlignment','bottom');
% make all text in the figure to size 14 and bold
set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold')

end