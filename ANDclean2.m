ndm = 1272;
dmcutoff = 50;
dt = 0.000064;
tbinw = 0.1;
tobs = 135;
tbins = tbinw/2:tbinw:round(tobs/tbinw)*tbinw;
ttol = tbinw/2;
maxpenalty = 2;
nfiles = 7;

penalties = zeros(length(tbins),7);
zapranges = zeros(length(tbins),7);
%lastbeam = zeros(length(tbins),7)-1;

delayfile = 'delays';
[idmdelay,zero,nsamp] = textread(delayfile,'%d%d%d');

%Construct zap ranges for all beams
for ninfile = 0:nfiles-1
    infile0 = sprintf('%s%d','sift.dat.t_sorted',ninfile)
    [idm0,zero,zero,it0,snr0,o10,o20,o30,o40,o50,o60] = textread(infile0,'%d%d%d%d%f%f%f%f%d%f%f');
    %infile0 = sprintf('%s%d','best.match.t_sorted',ninfile)
    %[idm0,zero,zero,it0,snr0,o10,o20] = textread(infile0,'%d%d%d%d%f%f%f');

    it_corrected0 = it0 - round(nsamp(idm0+1)/2);

    %Add SNRs of events in the same bin
    for i=1:length(tbins)
        ii = find(idm0 < dmcutoff & it_corrected0*dt > tbins(i)-tbinw/2 & ...
            it_corrected0*dt < tbins(i)+tbinw/2);
        zapranges(i,ninfile+1) = sum(snr0(ii));
    end

end


%Assign penalties to time bins depending on how many beams an events was
%detected in
for abeam = 0:nfiles-1
    for bbeam=abeam+1:nfiles-1
        andzap = zapranges(:,abeam+1) & zapranges(:,bbeam+1);
        if abeam == 0
            penalty = 1;
        elseif bbeam-abeam > 3
            penalty = 6-(bbeam-abeam);
        else
            penalty = bbeam-abeam;
        end
        
        ibad = find(andzap == 1);
        penalties(ibad,abeam+1) = penalties(ibad,abeam+1) + penalty;
        penalties(ibad,bbeam+1) = penalties(ibad,bbeam+1) + penalty;
    end
end

%surf(penalties),shading flat;

%Read pulse files again, zap events within zap ranges & with high penalties
for ninfile = 0:nfiles-1
    infile0 = sprintf('%s%d','sift.dat.t_sorted',ninfile)
    [idm0,zero,zero,it0,snr0,o10,o20,o30,o40,o50,o60] = textread(infile0,'%d%d%d%d%f%f%f%f%d%f%f');
    nevents = length(idm0);

    it_corrected0 = it0 - round(nsamp(idm0+1)/2);
    subplot(nfiles,1,ninfile+1), plot(it_corrected0*dt,idm0,'bo'),hold on;
    
    ibad = find(penalties(:,ninfile+1) > maxpenalty);
    for i=1:length(ibad)
        tzap = tbins(ibad(i));
        %subplot(nfiles,1,ninfile+1),plot([tzap-ttol tzap-ttol],[0 1272],'r');
        %subplot(nfiles,1,ninfile+1),plot([tzap+ttol tzap+ttol],[0 1272],'r');
        
        ii = find(it_corrected0*dt > tzap-ttol & it_corrected0*dt < tzap+ttol);
        subplot(nfiles,1,ninfile+1),plot(it_corrected0(ii)*dt,idm0(ii),'ro');
        idm0(ii) = -1;
        axis([0 tobs 0 ndm]);
    end
    
    ylabel('iDM');
    title(sprintf('%s%d %s%d %s%d','(BEFORE)Beam:',ninfile,'Tot.events:',length(idm0),...
        'Good events:',length(idm0)-length(find(idm0 == -1))));
end


% figure;
% for ninfile = 0:nfiles-1
%     subplot(nfiles,1,ninfile+1),plot(tbins,penalties(:,ninfile+1));
% end
