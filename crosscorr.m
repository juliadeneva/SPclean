delayfile = 'delays';
infile0 = 'sift.dat.t_sorted0';
infile1 = 'sift.dat.t_sorted6';

[idmdelay,zero,nsamp] = textread(delayfile,'%d%d%d');
[idm0,zero,zero,it0,snr0,o10,o20,o30,o40,o50,o60] = textread(infile0,'%d%d%d%d%f%f%f%f%d%f%f');
[idm1,zero,zero,it1,snr1,o11,o21,o31,o41,o51,o61] = textread(infile1,'%d%d%d%d%f%f%f%f%d%f%f');

n0 = length(idm0);
ndm = 1272;
dmcutoff = 100;
dt = 0.000064;
tbinw = 0.01;
tobs = 135;
tbins = tbinw/2:tbinw:round(tobs/tbinw)*tbinw;

it_corrected0 = it0 - round(nsamp(idm0+1)/2);
it_corrected1 = it1 - round(nsamp(idm1+1)/2);


subplot(2,1,1);
plot(it_corrected0*dt,idm0,'ro');
subplot(2,1,2);
plot(it_corrected1*dt,idm1,'ro');

%Get indices for events with idm~0 and bin events
dm0indices0 = find(idm0<dmcutoff);
[counts0] = hist(it_corrected0(dm0indices0)*dt,tbins);
dm0indices1 = find(idm1<dmcutoff);
[counts1] = hist(it_corrected1(dm0indices1)*dt,tbins);

%Add SNRs of events in the same bin
snrsum0 = zeros(length(tbins),1);
snrsum1 = zeros(length(tbins),1);
for i=1:length(tbins)
    if(counts0(i) ~= 0)
        ii = find(idm0 < dmcutoff & it_corrected0*dt > tbins(i)-tbinw/2 & ...
            it_corrected0*dt < tbins(i)+tbinw/2);
        snrsum0(i) = sum(snr0(ii));
    end
    
    if(counts1(i) ~= 0)
        ii = find(idm1 < dmcutoff & it_corrected1*dt > tbins(i)-tbinw/2 & ...
            it_corrected1*dt < tbins(i)+tbinw/2);
        snrsum1(i) = sum(snr1(ii));
    end
end

figure;
subplot(3,1,1),plot(tbins,snrsum0,'b');
subplot(3,1,2),plot(tbins,snrsum1,'r');
subplot(3,1,3),plot(tbins, (snrsum0 & snrsum1),'g')

% [c,lags] = xcorr(snrsum0,snrsum1);
% %[c,lags] = xcorr(counts0,counts1);
% figure;
% plot(lags,c);