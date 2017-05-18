sampleRate=1000;                 % samples per second (make even); must be at least twice the largest frequency of interest
timeMax=1;                       % seconds
stdNoise=1;                      % times Gaussian random
f1=200; f2=300;                  % Hz
amp1=1; amp2=2;                  % amplitude of the components
% ---------------------------------------------------
T=1/sampleRate;
N=round(timeMax/T); if 2*round(N/2)~=round(N), N=N+1; end
t=[0:1:N-1]*T; timeMax=max(t);
% ---------------------------------------------------
% ---------------------------------------------------
Nyquist=2*max([f1 f2]);
disp(['Sample Rate = ' num2str(sampleRate) ' Hz, Nyquist = ' num2str(Nyquist) ' Hz'])
if sampleRate<Nyquist, warning('need larger sampleRate'); disp('pause...'); pause; end
% ---------------------------------------------------
s=zeros(size(t));
noise=stdNoise*randn(size(t)); % Gaussian noise
phi1=rand(1,1)*2*pi; % random phase  
phi2=rand(1,1)*2*pi;
s=amp1*cos(2*pi*f1*t+phi1)+amp2*cos(2*pi*f2*t+phi2)+noise;
    
subplot(2,1,1); plot(t,s); xlabel('time [sec]'); 
 
% ---------------------------------------------------
F=fftshift(fft(s));
freq=[-N/2:1:N/2-1]*sampleRate/N;
% ---------------------------------------------------
maxAmp=max(abs(F)); % rescale the frequency amps for convenience
subplot(2,1,2); plot(freq,abs(F)/maxAmp); xlabel('frequency [Hz]'); ylabel('scaled');

