%% Process ultrasound A-mode signal for accuracy/repeatibility tests
% input: A-mode ultrasound signal

clear all; clc;
clear figure;

filename='signalch2.xlsx';
A=xlsread(filename);
A1=A(:,2);
t=A(:,1);
N=size(A,1);

%% %% define signal speed (ss) in medium 
% need to calculate this first before performing TGC - to determine ss 
% distance = (signal propagation time * signal speed in medium)/2

% define mode
% mode = input('1 - P/E, 2 - T. Define mode =  ');
mode = 1

% define medium
% medium = input('1 - water, 2 - others. Define medium =  ');
medium = 1

if medium==1;
    ss=1484; % signal speed in 20deg water = 1484m/s
else
% ave signal speed ssa = s_fat*(% of fat) + s_muscle*(% of muscle)
ssm=1630;    % signal speed in muscle = 1630m/s
ssf=1458;    % signal speed in fat = 1458m/s

sf = input('% of fat = ');
sm = input('% of muscle = ');

if sf+sm~=1;
    warning('percentage not equals 1');    
end;

ss = ssf*sf + ssm*sm;
end;

%%  Signal smoothing - low pass filter butterworth? 
% result B

freq = 2000000   ;
[b,a]=butter(2,50000/freq,'low');  
%[b,a]=butter(2,100000/freq,'high'); 
%[b,a]=butter(2,[500000 1500000]/freq); 
B=filtfilt(b,a,A1);

figure(1);
plot(A1,'y'); hold on; plot(B); hold off; ylim([-1,1]);
title('Low pass filter');

%% %% create time gain compensation (TGC) function based on attenuation value
% result A_f

% transmit frequency, and round trip distance
tgc_alpha = 0.10;		% [dB/(MHz cm)]
freq = 2000000   ;          % [2 MHz]
% r = 0.15        ;               % [m] wrt to time

% distance = (time*speed)/2
% ss=1484;  % speed of sound in water at 20deg
i=1;
while i<=size(A,1)
r(i,1)=(A(i,1)*ss)/2;
i=i+1;
end;

% tgc = exp(2*tgc_alpha*freq*(r));
tgc = exp(tgc_alpha*freq/1e6*r*100);

j=1;
while j<=size(B,1)
A_f(j,1)=B(j)*tgc(j,1);  %A_f - TGCed signal to be used
j=j+1;
end;

figure(2);
subplot(3,1,1); plot(A(:,1),A(:,2)); title('Raw Signal'); ylim([-2.0 2.0]);
subplot(3,1,2); plot(tgc); ylabel('TGC');
subplot(3,1,3); plot(A(:,1),A_f(:,1),'k'); title('Resulting Signal'); ylim([-2.0 2.0]); hold on;


%% Determine start point - zero crossing
 
% find_zero = diff(sign(y')) %pointless hilbert transform does not cross
% zero

find_zero = diff(sign(A_f));
leng = length(find_zero);
b=1;

for b=1:1:leng;
    if find_zero(b) == -2;
        dist = (t(b)*1484)/2;          % distance from ultrasound signal
    end;
end;

%% pick the right range
er_zone=0.05;  % in m, 2 cm
er_time=er_zone*2/1484;
er_time = str2num(sprintf('%.8f',er_time)); %round to 8 decimal places
loc = find(A(:,1) == er_time);

find_zero1 = diff(sign(A_f(loc:end)));

for b=loc:1:leng;
    if find_zero1(b) == 2;
        dist1 = (t(b)*1484)/2+er_zone;
        disp(dist1);% *first* distance from ultrasound signal
        break;
    end;
end;

figure(3);
subplot(3,1,1); plot(t(loc:end,1),A(loc:end,2)); title('The Originals'); ylim([-0.2 0.2]);
subplot(3,1,2); plot(t(loc:end,1),A_f(loc:end,1),'r'); title('Filtered Signal'); ylim([-0.2 0.2]);
subplot(3,1,3); plot(t((loc:leng),1),find_zero(loc:end,1),'k'); title('Zeros-Pos-Neg'); ylim([-2.0 2.0]); 
title('Zero crossing');

