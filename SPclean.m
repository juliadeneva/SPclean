%infile = 'J1928+15.best.t_sorted';
infile = 'J1928+15.sift.dat.t_sorted';
%infile = 'B0301+19.sift.dat.t_sorted';
%infile = 'sift.dat.t_sorted';
outfile = 'pulses.good';
delayfile = 'delays';

[idm2,zero,nsamp] = textread(delayfile,'%d%d%d');
[idm,zero,zero,it,snr,o1,o2,o3,o4,o5,o6] = textread(infile,'%d%d%d%d%f%f%f%f%d%f%f');
%[idm,zero,nsm,it,snr,o1,o2] = textread(infile,'%d%d%d%d%f%f%f');

n = length(idm);
ndm = 1272;
dt = 0.000064;

evcounts = zeros(ndm,1);
it_corrected = zeros(n,1);
ttol = 0.1;
fractol = 0.3;

%plot(it,idm,'bo'),hold on;
for i=1:n
    it_corrected(i) = it(i) - nsamp(idm(i)+1)/2;
end
t_corrected = it_corrected*dt;
%t_corrected = it*dt;

subplot(2,1,1);
plot(t_corrected,idm,'ro'), hold on;
ylabel('iDM');
title('Single-blast cleaning');

bookmark = 0;
for i=1:n
   if(i<n & (t_corrected(i+1)-t_corrected(i) < ttol))
       evcounts(idm(i)+1) = evcounts(idm(i)+1)+1;
       %evcounts(idm(i)+1) = 1;
       if(bookmark == 0)
           bookmark = i; 
       end
   else
       if(sum(evcounts)/ndm > fractol)
            for j=bookmark:i
               idm(j) = -1;
            end
            plot([t_corrected(bookmark) t_corrected(bookmark)],[0 1272],'b');
            plot([t_corrected(i) t_corrected(i)],[0 1272],'g');
       end
       bookmark = 0;       
       evcounts = zeros(ndm,1);
   end
end

igoods = find(idm >= 0);
tgoods = t_corrected(igoods);
idmgoods = idm(igoods);
subplot(2,1,2);
plot(tgoods,idmgoods,'ro');
xlabel('t(s)');
ylabel('iDM');
