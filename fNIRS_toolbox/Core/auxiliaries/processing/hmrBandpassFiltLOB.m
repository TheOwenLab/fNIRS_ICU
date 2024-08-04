% y2 = hmrBandpassFilt( y, fs, hpf, lpf )
%
% UI NAME:
% Bandpass_Filter
%
% y2 = hmrBandpassFilt( y, fs, hpf, lpf )
% Perform a bandpass filter
%
% INPUT:
% y - data to filter #time points x #channels of data
% fs - sample frequency (Hz). If length(fs)>1 then this is assumed to be a time
%      vector from which fs is estimated
% hpf - high pass filter frequency (Hz)
%       Typical value is 0 to 0.02.
% lpf - low pass filter frequency (Hz)
%       Typical value is 0.5 to 3.
%
% OUTPUT:
% y2 - filtered data

function [y2,ylpf] = hmrBandpassFilt( y, fs, hpf, lpf )



% convert t to fs
% assume fs is a time vector if length>1
if length(fs)>1
    fs = 1/(fs(2)-fs(1));
end



% low pass filter
FilterType = 1;
FilterOrder = 3;
%[fa,fb]=butter(FilterOrder,lpf*2/fs);
if FilterType==1 | FilterType==5
    [fb,fa] = MakeFilter(FilterType,FilterOrder,fs,lpf,'low');
elseif FilterType==4
%    [fb,fa] = MakeFilter(FilterType,FilterOrder,fs,lpf,'low',Filter_Rp,Filter_Rs);
else
%    [fb,fa] = MakeFilter(FilterType,FilterOrder,fs,lpf,'low',Filter_Rp);
end
% Reorganize data with NaN so that filter can work 
[rk,ma] = fillmissing(y,'linear');
ylpf=filtfilt(fb,fa,rk);
for hb = 1:size(y,3)
    for chn = 1:size(y,2)
        clear lst
        lst = find( ma(:,chn,hb) == 1 );
        ylpf(lst,chn,hb) = NaN;
    end
end



% high pass filter
FilterType = 1;
FilterOrder = 5;
if FilterType==1 | FilterType==5
    [fb,fa] = MakeFilter(FilterType,FilterOrder,fs,hpf,'high');
elseif FilterType==4
%    [fb,fa] = MakeFilter(FilterType,FilterOrder,fs,hpf,'high',Filter_Rp,Filter_Rs);
else
%    [fb,fa] = MakeFilter(FilterType,FilterOrder,fs,hpf,'high',Filter_Rp);
end

if FilterType~=5
    % Reorganize data with NaN so that filter can work
    clear rk ma
    [rk,ma] = fillmissing(ylpf,'linear');
    y2=filtfilt(fb,fa,rk);
    for hb = 1:size(y,3)
        for chn = 1:size(y,2)
            clear lst
            lst = find( ma(:,chn,hb) == 1 );
            y2(lst,chn,hb) = NaN;
        end
    end
else
    y2 = ylpf;
end
