function [OUT]=Get_power_freq_PSD(csc_name,cleanLFP,dwnsmpl)
%saves artifact threshold- you could probably make a globals file to save
%this in
%Goal is to create a starter function that can be used to itterate and
%generate an average hippocampal PSD for saline and MK801 animals
if ischar(csc_name)
    region=csc_name;
    reg_name=sprintf('*%s*.ncs',region);
    file=find_files(reg_name); %Idea is to identify hippocampus ones or if not then refer to animal if they are standardized
    
else
    file={csc_name}; %CSC 4 thalamus, CSC 12 MSN, CSC 8 HIPP, 19 MPFC in acute dataset
    
    
end
%will want to downsample and detrend
%then remove artifacts
%then plot the PSD in an OUT file along with rat number, drug conditon
% clear all
PLOT_IT=1; %change value to 1 if you want to see the plots


%%%%%%%% Get Power for cannoical bands as well as peak frequency %%%%
%%% includes artifact rejection and visual inspection
if nargin<3
    dwnsmpl=0;
end
if nargin<2
    cleanLFP=0;
end
OUT.aborted=false; %This will become true if something is off in the data

% These are based on Nancy's labeling but we can adjust depending on the project or dataset


    
    %A channel translation table would useful here
    parts=strsplit(pwd,'\');
    OUT.animal=parts{9}(1:4); %This stuff will change depending on your directory, right now it is set up for the ACUTE MK801
    OUT.drug=parts{7};
    OUT.day=parts{8};
    
    %THIS is how it worked for Nancy's stuff
    file=find_files('*HIPP*.ncs'); %Idea is to identify hippocampus ones or if not then refer to animal if they are standardized
    if isempty(file) && str2double(OUT.animal)>=1042
        file=find_files('CSC8*.ncs'); %nancy says for all animals 1042 and on that csc 8 is hippocampus, 16 is mpfc
    elseif isempty(file) && str2double(OUT.animal)<=1042
        OUT.aborted=true;
        disp('no regional file found')
        return
    end



if dwnsmpl==1
    downsample_fq=1000; %I think the files are often recorded at 32k samples/sec, we can lower that
    
    [LFP,sFreq]=convert_dwnspl_detrend(file{1},downsample_fq); %second input is desired downsampled freq, if you want to convert to seconds make last arg 1
    
    endtime=LFP(end,1)-(2*60); %2 min before end
    starttime=LFP(end,1)-(7*60);  %make 5 min range
    ix=LFP(:,1)>starttime & LFP(:,1)<endtime;
    % LFP_i=LFP(ix,:);
    % figure
    % plot(LFP(ix,1),LFP(ix,2))
    % [~,y]=ginput(1);
    % OUT.art_thresh=y;
else
    [LFP,sFreq]=convert_dwnspl_detrend(file{1});
end
if cleanLFP==1
    [BIX,artifact_times_usec] = LD_Clean_LFP(LFP,[],6e4,sFreq); %Third argument is the threshold for artifact rejection
    perc_bad = sum(BIX)/length(BIX);
    fprintf('BAD percent: %2.2f\n',perc_bad*100)
    if perc_bad > .3 %if its grabbing 10 minutes this means that ill get 6
        disp('Too much bad data')
        OUT.aborted=true;
        return
    end
    
    newLFP=LFP(~BIX,:);
    clear LFP
    LFP=newLFP;
end
% THIS REMOVE 60hz- note this leaves a bit dip in the PSD at 60, thats fine
% but in case it looks weird that is why
d=Notch_filter_cowen(sFreq);
notched_LFP=filtfilt(d,LFP(:,2));
LFP(:,2)=notched_LFP;
clear notched_LFP

disp('A notch filter was applied to remove 60Hz')

if PLOT_IT==1
    figure
    plot(LFP(:,1),LFP(:,2))
    xlabel('Time')
    ylabel('microvolts')
    if isfield(OUT,'animal')
        title(sprintf('Animal %s %s %s',OUT.animal, OUT.drug, OUT.day))
    else
        title('Raw data')
    end
end

Theta_filt = designfilt('bandpassiir','FilterOrder',12, ...  %YOU CAN ADJUST THE BANDS ON HERE- IE THETA IS from 5 to 12 (second and third white entry) but can be changed
    'HalfPowerFrequency1',5, 'HalfPowerFrequency2',12, ...
    'SampleRate',sFreq,'designmethod', 'butter');
HighGamma = designfilt('bandpassiir','FilterOrder',12, ...
    'HalfPowerFrequency1',130, 'HalfPowerFrequency2',180, ... %numbers are based on acute mK801
    'SampleRate',sFreq,'designmethod', 'butter');
LowGamma = designfilt('bandpassiir','FilterOrder',12, ...
    'HalfPowerFrequency1',30, 'HalfPowerFrequency2',50, ...
    'SampleRate',sFreq,'designmethod', 'butter');
Delta=designfilt('bandpassiir','FilterOrder',12, ...
    'HalfPowerFrequency1',0.5, 'HalfPowerFrequency2',3, ...
    'SampleRate',sFreq,'designmethod', 'butter');
Broad_gamma=designfilt('bandpassiir','FilterOrder',12, ...
    'HalfPowerFrequency1',30, 'HalfPowerFrequency2',120, ...
    'SampleRate',sFreq,'designmethod', 'butter');


% se=pentropy(newLFP(:,2),downsample_fq,'FrequencyLimits',[30 150]);


hgpow=abs(hilbert(filtfilt(HighGamma,LFP(:,2))));
lgpow=abs(hilbert(filtfilt(LowGamma,LFP(:,2))));
thetapow=abs(hilbert(filtfilt(Theta_filt,LFP(:,2))));
deltapow=abs(hilbert(filtfilt(Delta,LFP(:,2))));
broadpow=abs(hilbert(filtfilt(Broad_gamma,LFP(:,2))));


OUT.broad_pow=nanmean(broadpow);  %anytime you want to save anything you just add it to the OUT structure
OUT.hg_delt=nanmean(hgpow./deltapow);
OUT.lg_delt=nanmean(lgpow./deltapow);
OUT.theta_delt=nanmean(thetapow./deltapow);
OUT.broad_delt=nanmean(broadpow./deltapow);
OUT.raw_theta=nanmean(thetapow);
OUT.raw_hg=nanmean(hgpow);
OUT.raw_lg=nanmean(lgpow);



% Here is a way to have a look at the PSD- this is worth investigating to
% spot any irregularities
[pxx,f] =pmtm(LFP(:,2),5,[1:0.5:120],sFreq); %The middle argument is how you chose the frequency range you want - format: you want it starting at X in increments of X to X
dbpxx=10*log10(pxx); %turns to decibels, reduces variablity but still you really want this relative to some baseline

OUT.psd=pxx;
OUT.dbpsd=dbpxx;
OUT.freqs_forPSD=f;

if PLOT_IT==1  %note if you don't downsample and have a long recording this will take a long time to plot
    figure
    plot(f,dbpxx)
    ylabel('Decibels')
    xlabel('Frequency')
    ylim([10 80])%% this set the scale. Advantage is all the figures you plot will have same scale- disadvantage is if some are out of the range, you can change it if so
    if isfield(OUT,'animal')
        title(sprintf('Animal %s %s %s',OUT.animal, OUT.drug, OUT.day))
    else
        title('Power Spectral Density (PSD)')
    end
    
end

% FINDING PEAK FREQUENCIES
[thetaFreq,t]=instfreq(LFP(:,2),sFreq,'FrequencyLimits',[5 10]); %instafreq seems to want to work in seconds
thetafreq=mean(thetaFreq);
OUT.thetafreq=thetafreq;

[lowgammafreq,t]=instfreq(LFP(:,2),sFreq,'FrequencyLimits',[30 50]); %instafreq seems to want to work in seconds
lowgammaf=mean(lowgammafreq);
OUT.lowgammafrex=lowgammaf;
%low gamma= 30-50
%high gamma= 130-180
[highgammafreq,t]=instfreq(LFP(:,2),sFreq,'FrequencyLimits',[130 180]); %instafreq seems to want to work in seconds
highgammaf=mean(highgammafreq);
OUT.highgammafrex=highgammaf;



