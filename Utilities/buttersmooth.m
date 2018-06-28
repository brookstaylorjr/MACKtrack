function traj_out = buttersmooth(traj_in,cutoff, butterOrder)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% traj_out  = buttersmooth(traj_in,cutoff, butterOrder)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% BUTTERSMOOTH performs Butterworth filtering on an input (e.g. expression) trajectory.
%
% traj_in     input trajectory (1-D vector - starting and ending NaNs will be ignored
% cutoff      low pass (if numel==1) or bandpass (if numel==2) frequency limits of Butterworth filter. Default = 0.1
% butterOrder order of Butterworth filter. default = 1.
%
% Code somewhat stolen from Michael Zhao. MICHAEL I TOOK YOUR COOOOOOOOODE
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin<3
    butterOrder = 2;
    if nargin<2
        cutoff = 0.05;
    end
end

% Get boundary NaN indices
firstpoint=find(~isnan(traj_in),1,'first');
lastpoint=find(~isnan(traj_in),1,'last');
temp=traj_in(firstpoint:lastpoint);


% Interpolate over any remaining NaNs
if sum(isnan(temp))>0
    t_tmp = find(~isnan(temp));
    temp = interp1(t_tmp,temp(~isnan(temp)),1:length(temp));
end


% Smooth with butterworth filter with specified order and cutoff
[B,A] = butter(butterOrder,cutoff,'low');
smoothtemp = filtfilt(B,A,temp);


% Return smoothed data (replace non-NaNs in original)
traj_out = traj_in;
traj_out(firstpoint:lastpoint)=smoothtemp;