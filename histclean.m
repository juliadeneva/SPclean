%infile = 'best.t_sorted';
%infile = 'J1928+15.sift.dat.t_sorted';
%infile = 'B0301+19.sift.dat.t_sorted';
%infile = 'sift.dat.t_sorted0';
infile = 'p2030_53497_28594_0013_G35.78-01.42.N_0.wapp+.pulse_0.sift.dat';

delayfile = 'delays';

[idm2,zero,nsamp] = textread(delayfile,'%d%d%d');
[idm,zero,zero,it,snr,o1,o2,o3,o4,o5,o6] = textread(infile,'%d%d%d%d%f%f%f%f%d%f%f');
%[idm,zero,nsm,it,snr,o1,o2] = textread(infile,'%d%d%d%d%f%f%f');

n = length(idm);
ndm = 1272;
dt = 0.000064;
fradar = 1.0/12;
nit = 2;
nsigma = 7;
it_corrected = zeros(n,1);
ttol = 0.2;

for i=1:n
    it_corrected(i) = it(i) - round(nsamp(idm(i)+1)/2);
end

%Get indices for events with idm=0
dm0indices = find(idm<10);
%Bin events by 0.1s
[n,tbin] = hist(it_corrected(dm0indices)*dt,2000);

%Find sigma, removing peaks over several iterations
sigma = std(n)
for i=1:nit    
    nlow = n(find(n/sigma < nsigma));
    %length(nlow)
    sigma = std(nlow)
end

%Find the true peaks
ipeaks = find(n/sigma > nsigma);
tpeaks = tbin(ipeaks);
npeaks = n(ipeaks);

%FFT them
df = 0.001;
maxf = 0.1;
f = 0.01:df:maxf;
dft = [];

for i=1:length(f)
    dftterms = exp(-2*pi*sqrt(-1)*tpeaks*f(i));
    dft = [dft sum(dftterms)];
end

dft2 = dft.*conj(dft);
dft2 = dft2/dft2(1);
%dft2 = dft2(2:length(dft2));
%dft = dft(2:length(dft));
%f = f(2:length(f));
%plot(f,dft2),hold on,plot(f,dft2,'bo');

%See if there's a peak
rms = std(dft2);
[dft2peak,ipeak] = max(dft2/rms);
f(ipeak)

ampall = sqrt(dft.*conj(dft));
phaseall = log(dft./ampall)/(-sqrt(-1));

%Get the phase from the radar fundamental in the DFT
%dftpeak = dft(ipeak-2:ipeak+2);
dftpeak = dft(ipeak);
amp = sqrt(dftpeak.*conj(dftpeak));
weights = amp/sum(amp);
phase = sum(weights.*log(dftpeak./amp)/(-sqrt(-1)))
%fradarw = sum(weights.*f(ipeak-2:ipeak+2))
fradarw = sum(weights.*f(ipeak))
tradar = phase/(fradarw*2*pi)

%Now see where the zapping windows are
figure, subplot(2,1,1);
plot(it_corrected*dt,idm,'bo');
hold on;
i = 0;
tzap = tradar;
while(tzap < max(it_corrected*dt))
    plot([tzap tzap],[0 1272],'r');
    ii = find(abs(it_corrected*dt-tzap)<ttol);
    %plot(it_corrected(ii)*dt,idm(ii),'ro');
    idm(ii) = -1;
    i = i+1;
    tzap = tradar + i/fradarw;    
end
ylabel('iDM');
title('DFT cleaning');

subplot(2,1,2);
igood = find(idm > -1);
plot(it_corrected(igood)*dt,idm(igood),'bo');
xlabel('t(s)');
ylabel('iDM');

figure;
subplot(2,1,1);
plot([fradarw fradarw],[min(ampall) max(ampall)],'r'),hold on;
plot(f,ampall,'bo');
xlabel('f(Hz)');
ylabel('Amp');
subplot(2,1,2);
plot([fradarw fradarw],[min(real(phaseall)) max(real(phaseall))],'r'),hold on;
plot(f,phaseall,'bo');
xlabel('f(Hz)');
ylabel('Phase');
