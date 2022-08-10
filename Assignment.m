clc
close all
clear

% Load the data sample
load('Sample_1.mat');

% Divide each sample data by 250 to make the data smaller
ECGsignal = (Orig_Sig ./ 250);

% Set the time to 10 seconds with the length of the data sample
x = linspace(1, 10, length(ECGsignal));

% Calculate the root mean square of the data sample
T = length(x);
f = ECGsignal .^ 2;
RMS = rms(f, T);

% Filter the output sample using a low and high filter with chosen Hz
[b, a] = butter(1, 0.022, 'Low');
y1 = filtfilt(b, a, RMS);
[b, a] = butter(1, 0.01, 'High');
y = filtfilt(b, a, y1);

% Plot the graph
grid on
grid minor
hold on
plot(x,y)

% Set the variable to the length of x
Sample_Index_Number = 1 : length(x);

% Finds the mean of the peaks in data sample to find the actual R peak
MEANPEAK = mean(findpeaks(y, Sample_Index_Number));
AvgPeakHeight = MEANPEAK + 0.5;

% Find the R peaks
[PKS,LOCS] = findpeaks(y, Sample_Index_Number);

LOCS_PKS_Index = find(PKS > AvgPeakHeight);

% Find R Wave
Rate = zeros(1, length(LOCS_PKS_Index));
for L = 1 : length(LOCS_PKS_Index)
    Rate(L) = x(LOCS(LOCS_PKS_Index(L)));
    plot(x(LOCS(LOCS_PKS_Index(L))), PKS(LOCS_PKS_Index(L)), 'rv', 'MarkerFaceColor', 'r');
end

% Find S Wave
TempIndexS = zeros(length(LOCS_PKS_Index), 1);
RecordS = zeros(length(LOCS_PKS_Index), 1);
for L = 1 : length(LOCS_PKS_Index)
    TempIndexS(L) = LOCS((LOCS_PKS_Index(L)));
    
    FR = y(TempIndexS(L) + 1);
    BR = y(TempIndexS(L));
    
    while 1
        if (TempIndexS(L) + 1) < length(y)
            if (BR > FR)
                TempIndexS(L) = TempIndexS(L) + 1;
                BR = FR;
                FR = y(TempIndexS(L) + 1);
            else
                break;
            end
        else
            break;
        end
    end
    
    plot(x(TempIndexS(L)), BR, 'k*', 'MarkerFaceColor', 'k')
    
    RecordS(L) = x(TempIndexS(L));
end

% Find Q Wave
TempIndexQ = zeros(length(LOCS_PKS_Index), 1);
RecordQ = zeros(length(LOCS_PKS_Index), 1);
for L = 1 : length(LOCS_PKS_Index)
    TempIndexQ(L) = LOCS((LOCS_PKS_Index(L)));
    
    FR = y(TempIndexQ(L) - 1);
    BR = y(TempIndexQ(L));
    
    while 1
        if (TempIndexQ(L) - 1) > 0
            if (BR > FR)
                TempIndexQ(L) = TempIndexQ(L) - 1;
                BR = FR;
                FR = y(TempIndexQ(L) - 1 );
            else
                break;
            end
        else
            break;
        end
    end
    
    plot(x(TempIndexQ(L)), BR, 'og', 'MarkerFaceColor', 'g')
    
    RecordQ(L) = x(TempIndexQ(L));
end

% Find peak after P Wave
TempIndexPf = zeros(length(TempIndexQ), 1);
for L = 1 : length(TempIndexQ)
    TempIndexPf(L) = TempIndexQ(L);
    
    FR = y(TempIndexPf(L) - 1);
    BR = y(TempIndexPf(L));
    
    while 1
        if (TempIndexPf(L) - 1) > 0
            if (BR < FR)
                TempIndexPf(L) = TempIndexPf(L) - 1;
                BR = FR;
                FR = y(TempIndexPf(L) - 1);
            else
                break;
            end
        else
            break;
        end
    end
end

% Find P Wave
RecordP = zeros(length(TempIndexPf), 1);
for L = 1 : length(TempIndexPf)
    TempIndexP = TempIndexPf(L);
    
    FR = y(TempIndexP - 1);
    BR = y(TempIndexP);
    
    while 1
        if (TempIndexP - 1) > 0
            if (BR > FR)
                TempIndexP = TempIndexP - 1;
                BR = FR;
                FR = y(TempIndexP - 1);
            else
                break;
            end
        else
            break;
        end
    end
    
    plot(x(TempIndexP), BR, 'hb', 'MarkerFaceColor', 'b')
    
    RecordP(L) = x(TempIndexP);
end

% Find peak before T Wave
TempIndexTp = zeros(length(TempIndexS), 1);
for L = 1 : length(TempIndexS)
    TempIndexTp(L) = TempIndexS(L);
    
    FR = y(TempIndexTp(L) + 1);
    BR = y(TempIndexTp(L));
    
    while 1
        if (TempIndexTp(L) + 1) < length(y)
            if (BR < FR)
                TempIndexTp(L) = TempIndexTp(L) + 1;
                BR = FR;
                FR = y(TempIndexTp(L) + 1);
            else
                break;
            end
        else
            break;
        end
    end
end

% Find T Wave
RecordT = zeros(length(TempIndexTp), 1);
for L = 1 : length(TempIndexTp)
    TempIndexT = TempIndexTp(L);
    
    FR = y(TempIndexT + 1);
    BR = y(TempIndexT);
    
    while 1
        if (TempIndexT + 1) < length(y)
            if (BR > FR)
                TempIndexT = TempIndexT + 1;
                BR = FR;
                FR = y(TempIndexT + 1);
            else
                break;
            end
        else
            break;
        end
    end
    
    plot(x(TempIndexT), BR, 'py', 'MarkerFaceColor', 'y')
    
    RecordT(L) = x(TempIndexT);
end

% BPM Calculation
BPM = mean(diff(Rate)) * 6 * length(LOCS_PKS_Index);

% Time Interval Calculation
PR_difference = RecordQ - RecordP;
PR_Interval = mean(PR_difference);

QT_difference = RecordT - RecordQ;
QT_Interval = mean(QT_difference);

QS_difference = RecordS - RecordQ;
QRS_duration = mean(QS_difference);

% Disolay BPM in command window
fprintf("BPM: %f beats per minute!\n", BPM);

% Check is PR Interval is normal
if (PR_Interval >= 0.12 && PR_Interval <= 0.2)
    fprintf("PR Interval: %f - Normal\n", PR_Interval);
else
    fprintf("PR Interval: %f - Abnormal\n", PR_Interval);
end

% Check is QT Interval is normal
if (QT_Interval <= 0.38)
    fprintf("QT Interval: %f - Normal\n", QT_Interval);
else
    fprintf("QT Interval: %f - Abnormal\n", QT_Interval);
end

% Check is QRS Duration is normal
if (QRS_duration <= 0.1)
    fprintf("QRS duration: %f - Normal\n", QRS_duration);
else
    fprintf("QRS duration: %f - Abnormal\n", QRS_duration);
end