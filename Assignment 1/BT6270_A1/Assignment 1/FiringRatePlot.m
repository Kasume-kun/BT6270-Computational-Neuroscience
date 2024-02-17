close all;
stepSize=0.01; % Kept low for faster computation
freq=zeros(1,60);
ct=1;

for ImpCur=0:stepSize:1
    gkmax=.36;
    vk=-77;
    gnamax=1.20;
    vna=50;
    gl=0.003;
    vl=-54.387;
    cm=.01;
    
    dt=0.01;
    niter=50000;
    t=(1:niter)*dt;
    iapp=ImpCur*ones(1,niter);
    v=-64.9964;
    m=0.0530;
    h=0.5960;
    n=0.3177;
    
    gnahist=zeros(1,niter);
    gkhist=zeros(1,niter);
    vhist=zeros(1,niter);
    mhist=zeros(1,niter);
    hhist=zeros(1,niter);
    nhist=zeros(1,niter);
    
    
    for iter=1:niter
      gna=gnamax*m^3*h;
      gk=gkmax*n^4;
      gtot=gna+gk+gl;
      vinf = ((gna*vna+gk*vk+gl*vl)+ iapp(iter))/gtot;
      tauv = cm/gtot;
      v=vinf+(v-vinf)*exp(-dt/tauv);
      alpham = 0.1*(v+40)/(1-exp(-(v+40)/10));
      betam = 4*exp(-0.0556*(v+65));
      alphan = 0.01*(v+55)/(1-exp(-(v+55)/10));
      betan = 0.125*exp(-(v+65)/80);
      alphah = 0.07*exp(-0.05*(v+65));
      betah = 1/(1+exp(-0.1*(v+35)));
      taum = 1/(alpham+betam);
      tauh = 1/(alphah+betah);
      taun = 1/(alphan+betan);
      minf = alpham*taum;
      hinf = alphah*tauh;
      ninf = alphan*taun;
      m=minf+(m-minf)*exp(-dt/taum);
      h=hinf+(h-hinf)*exp(-dt/tauh);
      n=ninf+(n-ninf)*exp(-dt/taun);
      vhist(iter)=v; mhist(iter)=m; hhist(iter)=h; nhist(iter)=n;
    end
    j=1;
    peaks=zeros;
    [mag,loc]=findpeaks(vhist);
    for temp=1:length(mag)
        if mag(temp) >=10 % Min value for a spike to be considered a AP
            peaks(j)=mag(temp);
            j=j+1;
        end
    end
    if peaks ~= 0
        no_of_peaks(ct)=length(peaks);
    else
        no_of_peaks(ct)=0;
    end
    ct=ct+1;
end
I1=0.0224;   % Calculated at a stepsize of 0.0001
I2=0.0625;   
I3=0.4578;   
freq=no_of_peaks*1000/(niter/100);
plot(0:stepSize:1,freq,'Color','#21DECC','LineWidth',1.2);
xlim([0 0.7]);
xlabel('I_{Ext}');
ylabel('No of Spikes per second');
hold on;
I1=0*freq+I1;
plot(I1,freq,'Color','#DE2133','LineWidth',1.2);
text(I1(1),-3,'I1','Color','#DE2133');
I2=0*freq+I2;
plot(I2,freq,'Color','#DE2133','LineWidth',1.2);
text(I2(1),-3,'I2','Color','#DE2133');
I3=0*freq+I3;
plot(I3,freq,'Color','#DE2133','LineWidth',1.2);
text(I3(1),-3,'I3','Color','#DE2133');