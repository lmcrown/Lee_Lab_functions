%% MK801_Question_Itterator_nonfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Analysis_function=@Get_power_freq_PSD;

% this will change based on what functions are available and how
% userfriendly it is and we want it to be
%
fun1=Analysis_function;
fun2=@Get_power_freq_PSD;
if isequal(fun1,fun2)
    r=inputdlg('enter target region (must match a string in the file- like HIPP or THAL');
    region=r{1};
    c=inputdlg('do you want to get rid of artifacts? Enter 1 for yes, 0 for no');
    cleanLFP=c{1};
    d=inputdlg('do you want to speed up processing by reducing sample rate?');
    dwnsmpl=d{1};
    input1=region;
    input2=logical(cleanLFP);
    input3=logical(dwnsmpl);
end
% you can also set these permanent by just putting them in here like
% "inputs=[___,___,___) and commneting out the dialog box

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

root_results_directory = 'F:\Keck Medicine of USC\MK801_Results\Acute';
specific_results_directory = fullfile(root_results_directory, func2str(Analysis_function));
if ~exist(specific_results_directory,'dir')
    mkdir(specific_results_directory)
end
root_data_directory = 'F:\Keck Medicine of USC\Zepeda, Nancy - 20946 ALL EEG FILES AGAIN\20946-Acute\EEG\';
cd(root_data_directory)

time=dir(root_data_directory);
good_ix = [];
for ii = 1:length(time)
    if time(ii).isdir && ~any(strfind(time(ii).name,'.') )
        good_ix = [good_ix, ii];
    end
end
time2=time(good_ix);
for itime=1:length(time2)
    tdir=fullfile(root_data_directory,time2(itime).name);
    s=strsplit(time2(itime).name,'-');
    tpoint=s{3};
    cd(tdir)
    drugs=dir(tdir);
    good_ix = [];
    for ii = 1:length(drugs)
        if drugs(ii).isdir && ~any(strfind(drugs(ii).name,'.') )
            good_ix = [good_ix, ii];
        end
    end
    drug2 = drugs(good_ix); %only rat files that actually have crap in them
    % now you just have the 2 drug folders
    for idir=1:length(drug2)
        drugdir=fullfile(tdir,drug2(idir).name);
        cd(drugdir)
        day=dir(fullfile(tdir,drug2(idir).name));
        good_ix = [];
        for ii = 1:length(day)
            if day(ii).isdir && ~any(strfind(day(ii).name,'.') )
                good_ix = [good_ix, ii];
            end
        end
        day2 = day(good_ix);
        for iday=1:length(day2)
            daydir=fullfile(drugdir,day2(iday).name);
            da=day2(iday).name;
            cd(daydir)
            
            good_ix2 = [];
            rats=dir(daydir);
            for ii = 1:length(rats)
                if rats(ii).isdir && ~any(strfind(rats(ii).name,'.') )
                    good_ix2 = [good_ix2, ii];
                end
            end
            rats2=rats(good_ix2);
            
            for irat=1:length(rats2)
                cd(fullfile(daydir,rats2(irat).name))
                [Dset]=Analysis_function(input1,input2,input3);
                if Dset.aborted==false
                    dayresults=fullfile(root_results_directory, func2str(Analysis_function),da);
                    if ~exist(dayresults,'dir')
                        mkdir(dayresults)
                    end
                    save(fullfile(dayresults,[Dset.animal '_' tpoint]) ,'Dset') % i want to make a folder for the day
                end
                cd ..
            end
            cd ..
            cd ..
        end
    end
end
msgbox('All done')
