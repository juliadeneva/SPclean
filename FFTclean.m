%infile = 'best.t_sorted';
infile = 'sift.dat.t_sorted';
delayfile = 'delays';

[idm2,zero,nsamp] = textread(delayfile,'%d%d%d');
[idm,zero,zero,it,snr,o1,o2,o3,o4,o5,o6] = textread(infile,'%d%d%d%d%f%f%f%f%d%f%f');
%[idm,zero,nsm,it,snr,o1,o2] = textread(infile,'%d%d%d%d%f%f%f');

n = length(idm);
ndm = 1272;
dt = 0.000064;
tbin = 0.1; %for binning events
tbinsamp = round(tbin/dt);
fradar = 1.0/12;

evcounts = zeros(ndm,1);
it_corrected = zeros(n,1);
ittol = 1000;

%Remove DM-dependent slanting
it_corrected = it - round(nsamp(idm+1)/2);
%Bin events into 0.1s-bins
it_corrected2 = round(it_corrected/tbinsamp)*tbinsamp;

%Get indices for events with idm=0
dm0indices = find(idm<1);
%Fake time series
t = 0:dt:max(it_corrected)*dt;
ts = zeros(max(it_corrected)+1,1);
for i=1:length(dm0indices)
    ii = dm0indices(i);
    %iti = it_corrected(ii);
    iti = it_corrected2(ii);
    ts(iti) = snr(ii);
end

%Get the spectrum bin corresponding to radar fundamental
df = 1/max(t);
f = 0:df:0.5/dt-df;
spec = fft(ts);
spec = spec(1:length(f));
spec2 = spec.*conj(spec);
plot(spec2);

binradar = round(fradar/df)+1;
%sumbins = sum(spec2(binradar-1:binradar+1));
%fradar = spec2(binradar-1)/sumbins*f(binradar-1) + ...
%    spec2(binradar)/sumbins*f(binradar) + ...
%    spec2(binradar+1)/sumbins*f(binradar+1)

dft = spec(binradar);
%dft = mean(spec(binradar-1:binradar+1));
amp = sqrt(dft*conj(dft));
phase = log(dft/amp)/sqrt(-1);
%tradar = abs(phase/(2*pi*fradar))
tradar = abs(phase/(2*pi*f(binradar)))

%Now see where the zapping windows are
figure;
plot(it_corrected2*dt,idm,'ro');
hold on;
i = 0;
tzap = tradar; 
while(tzap < max(it_corrected*dt))
    tzap = tradar + i/fradar;
    plot([tzap tzap],[0 1272],'b');
    i = i+1;
end
