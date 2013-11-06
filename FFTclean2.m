%infile = 'best.t_sorted';
infile = 'sift.dat.t_sorted';
delayfile = 'delays';

[idm2,zero,nsamp] = textread(delayfile,'%d%d%d');
[idm,zero,zero,it,snr,o1,o2,o3,o4,o5,o6] = textread(infile,'%d%d%d%d%f%f%f%f%d%f%f');
%[idm,zero,nsm,it,snr,o1,o2] = textread(infile,'%d%d%d%d%f%f%f');

n = length(idm);
ndm = 1272;
dt = 0.000064;
fradar = 1.0/12;

evcounts = zeros(ndm,1);
it_corrected = zeros(n,1);
ittol = 1000;

for i=1:n
    it_corrected(i) = it(i) - round(nsamp(idm(i)+1)/2);
end

%Get indices for events with idm=0
dm0indices = find(idm<2);
%DFT the supposed radar events
df = 0.01;
maxf = 1;
f = 0:df:maxf;
dft = [];
for i=1:length(f)
    iti = it_corrected(dm0indices);
    amp = snr(dm0indices);
    dftterms = exp(-2*pi*sqrt(-1)*iti*f(i));
    dft = [dft sum(dftterms)];
end

dft2 = dft.*conj(dft);
dft2 = dft2/dft2(1);
dft2 = dft2(2:length(dft2)/2-1);
dft = dft(2:length(dft)/2-1);
f = maxf/length(dft2):maxf/length(dft2):maxf;
plot(f,dft2),hold on,plot(f,dft2,'bo');

%See if there's a peak
rms = std(dft2);
[dft2peak,ipeak] = max(dft2/rms);
f(ipeak)
%dft2peak

%Is the found f(peak) a radar (sub)harmonic?
if(f(ipeak) > fradar)
    if(mod(f(ipeak),fradar) < 0.001)
        factor = round(f(ipeak)/fradar);
        fradar = f(ipeak)/factor;
        ipeak1 = round(ipeak/factor);
    else
        ipeak1 = ipeak;
    end
else
    if(mod(fradar,f(ipeak)) < 0.001)
        factor = round(fradar/f(ipeak));
        fradar = fradar/factor;
        ipeak1 = ipeak/factor;
    else
        ipeak1 = ipeak;
    end
end
f(ipeak1)

%Get the phase from the radar fundamental in the DFT
dftpeak = dft(ipeak1);
amp = sqrt(dftpeak*conj(dftpeak));
phase = log(dftpeak/amp)/(-sqrt(-1))
tradar = phase/(f(ipeak1)*2*pi)

%Now see where the zapping windows are
figure;
plot(it_corrected*dt,idm,'ro');
hold on;
i = 0;
tzap = tradar;
while(tzap < max(it_corrected*dt))
    tzap = tradar + i/f(ipeak1);
    plot([tzap tzap],[0 1272],'b');
    i = i+1;
end
