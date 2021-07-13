%% LFP processing instructions

%Possible Functions:
    %Get_power_freq_PSD
        %Possible arguments in:(csc_name,cleanLFP,dwnsmpl) %first is a CSC
        %name if you have just one or something, otherwise if itterating
        %through a large dataset you can set it to automatically find a
        %region of interest--- this i need to work on a little still
        %second is if you want to do artifact rejection and exclude any
        %datasets with too much artifact
        %third is if you want to reduce the sample rate, this will reduce
        %the file size and will speed up things if you want to process a
        %long bit of data
        
        %Support functions (you  must have these to run it with all its
        %functions)
            %convert_dwnsmpl_detrend - this will do the downsampling, it
            %will also change the .ncs file into a .mat file
            %by default it detrends which just means if the signal drifted
            %at all it will keep it at a mean of 0- this is useful when
            %using threshold and is pretty standard
            
            %LD_Clean_LFP- this performs the artifact rejection and will
            %reject an entire dataset if it contains too much artifact
            
            %Notch_filter_cowen- my old PI made this and its a little bit
            %better than the regular notch_filer- but to be fair all you
            %have to do to notch something is design a filter from 59.5 or
            %so to 60.5
            
            %nlx2matCSC_Matrix(EEG_filename,original_interval_usec)- this
            %is andaptation of a neuralynx function, you can use whatever
            %but they have a bunch of functions to interface with matlab
          