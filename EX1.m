clear; close all; clc;

load ('S1');
load ('S2');
Si = S1;
%%
fs = 10000; % sampling rate
dt = 1/fs; % time step
N = length (Si); % number of components in Si
t = dt*(1:N); % time vector
segment_length=0.3;
current_length=0.2;
dc_number=(N/fs)/segment_length; %number of segments
%%

TH=-30; %Threshold 
SaTH=zeros(1,N); 
for i=1:+1:N %finding when the Voltage is above the threshold
    if Si(i)>TH
        SaTH(i)=1;
    else
        SaTH(i)=0;
    end
end

SaTH_Changes=diff(SaTH);
L2H=find(SaTH_Changes==1); %crossing the threshold, low to high
H2L=find(SaTH_Changes==-1); %crossing the threshold, high to low

for i=1:+1:length(H2L) %loop for every spike
    [LocalMax,LocalMaxIndex]=max(Si((L2H(i):H2L(i))));
     LM(i)=LocalMax; %saves the value of each spike
     MaxIndices(i)=L2H(i)+LocalMaxIndex-1; %save the index of each spike
     SpikeTime(i)=t(MaxIndices(i)); %saves the time of each spike
end

%%
for i=1:+1:dc_number %loop for each segment
dc_time=i*segment_length; %the end time of each segment
SC(i)=length(find((SpikeTime <= dc_time) & (SpikeTime > (dc_time-segment_length)))); %counts the number of spikes in the segment
R(i)=SC(i)/current_length; %calculates the firing rate
end
%% graph
plot(t,Si)
hold on
title('intra-cellular recordings in a current-clamped cell')
grid on
set(gca,'FontSize',12);
xlabel('t (sec)')
ylabel('V (mV)')
ylim([-70,max(Si)+4]);
xlim([0,max(t)]);

plot(t(L2H),Si(L2H),'go')
plot(t(H2L),Si(H2L),'ro')
plot(SpikeTime,LM,'ko')
for i=1:+1:dc_number
text (i*segment_length-0.2, max(Si)+1.5, "R= " + num2str(R(i)) + " Hz ",'FontSize',8)
end
legend({'SI','L2H','H2L ','LM'},'Location','bestoutside')

