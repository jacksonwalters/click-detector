clear all; clc;

Program_folder=pwd;

% twowhales
% example1
% SW1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dialog boxes:

answer_Default_flag = questdlg('Open dialog box?:','','Yes','No','Other');
% Handle response
switch answer_Default_flag
    case 'Yes'
        Default_flag = 1;
    case 'No'
        Default_flag = 0;
end

if Default_flag
    
    prompt = {'Recording name:','Recording format:'};
    dlgtitle = 'Input'; 
    dims = [1 35];
    definput = {'channelABCD_2023-02-25_10-27-38','flac'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    files=dir([cell2mat(answer(1)) '.' cell2mat(answer(2))]);
    
    
    answer_Recording_Type = questdlg('Choose reciever type', ...
        'Recording type', ...
        'Bouy','Tag','Other');
    % Handle response
    switch answer_Recording_Type
        case 'Bouy'
            Tag_flag = 0;
        case 'Tag'
            Tag_flag = 1;
    end

    answer_Detector_Type = questdlg('Pick detector:','Detector type','Echolocation','Coda','Other');
    % Handle response
    switch answer_Detector_Type
        case 'Echolocation'
            Detector_flag = 0;
        case 'Coda'
            Detector_flag = 1;
    end

    answer_plot = questdlg('Visualize results?:','Plot activation','Yes','No','Other');
    % Handle response
    switch answer_plot
        case 'Yes'
            Plot_flag = 1;
        case 'No'
            Plot_flag = 0;
    end
else
    %% Pick Detector
    Detector_flag=0; % customized for buoy recievers: 1- apply coda detector | 0- apply echolocation clicks detector
    Plot_flag=1;     % 1- show detection figures | 0- dont show
    Tag_flag=0;    % 1- customized for Dtags: set to 1 for recordings from tags.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Check filename:
    
    Error_flag=length(files);
    if Error_flag==0
        error('File not found. Please check the file name (line 4)');        
    end
        
    %% Determine operational parameters (Global parameters)

    F_low_coda = 3e3;
    F_high_coda = 7e3;
    F_low_echo = 2e3;
    F_high_echo = 22e3;
    FsAnalyze = 48e3;
    T_sec=10; %[sec]                         % Define duration of window for analysis in seconds
    %%

    if Detector_flag
        F_low=F_low_coda;
        F_high=F_high_coda;
    else
        F_low=F_low_echo;
        F_high=F_high_echo;
    end
    
        
    for file_ind=1:length(files)
        
        filename=files(file_ind).name;
        FLAC=strfind(filename, '.flac');
        if ~isempty(FLAC)
            Name_save=filename(1:FLAC-1);
        else
            WAV=strfind(filename, '.wav');
            if ~isempty(FLAC)
                Name_save=filename(1:WAV-1);
            end
        end
        
        filename=files(file_ind).name;
        if ~Detector_flag
            Det_insert={'Buffer index','ToA[sec]','Whale index within buffer'};
            writecell(Det_insert,[Name_save '-ecolocation.csv'],'WriteMode','append');
            % filesave=[files(file_ind).folder filesep filename '-ecolocation.csv'];
        end
        file_to_read = [files(file_ind).folder,filesep,filename];
%         filesave=[Program_folder '\' Name '.xls'];
        Audio_name=[strfind(filename, 'flac') strfind(filename, 'wav')];
        if ~isempty(Audio_name)
            Timer=0;
            [y,Fs] = audioread(filename);                 % load recordings
            Y=y(:,1);                                     % Choose chanel one from the WRU
            if Fs < FsAnalyze
                S_factor = 1;
            elseif Detector_flag
                S_factor=floor(Fs/FsAnalyze);             % Define factor for resampling to 48khz
            else
                S_factor=floor(Fs/FsAnalyze);
            end        
            File_duration=(1/Fs)*(length(Y)-1);           % Calculate duration of the loaded recording
            Y_decimated = decimate(Y,S_factor);           % Resample recording to 48khz
            F_ds=Fs/S_factor;                             % Sample frequency of the decimated recording (48khz by default)       
            T=F_ds*T_sec;                                 % Define duration of window for analysis in samples
            T_raw=Fs*T_sec;
            NOI=floor(File_duration/T_sec);               % Calculate the number of windows in the current recording
            buffer_index=0; 
            Gather_TOA=[];
            Detected_SW=[];
            TOA_other_whale=[];
            TOA_other=[];
            tic
            for Buffer_ind=1:NOI
                Y_filtered=bandpass(Y_decimated(int32((Buffer_ind-1)*T+1):int32((Buffer_ind-1)*T+T)),[F_low, F_high],F_ds);     % Aply band pass filter and extract buffer                             
                if Tag_flag
                    [~,~,TOA,TOA_other]=Run_Detector_server(Y_filtered,Fs,F_ds,Detector_flag,Plot_flag,Tag_flag);  % Run detection
                else
                    [ToA_separated,TOA,~,~]=Run_Detector_server(Y_filtered,Fs,F_ds,Detector_flag,Plot_flag,Tag_flag);  % Run detection   
                    Gather(Buffer_ind).ToA_separated=ToA_separated;
                     NOW=length(Gather(Buffer_ind).ToA_separated);
                     for now=1:NOW
                        ToA_ind=Gather(Buffer_ind).ToA_separated{now} + 10*(Buffer_ind-1);
                        Lw=length(ToA_ind);
                        for w_ind=1:Lw
                           Det_insert={Buffer_ind,ToA_ind(w_ind),['SW' num2str(now)]};
                           if Detector_flag == 1
                               writecell(Det_insert,[Name_save '-codas.csv'],'WriteMode','append');
                           else
                               writecell(Det_insert,[Name_save '-ecolocation.csv'],'WriteMode','append');
                           end
                        end
                     end
                end

                Gather_TOA=[Gather_TOA (Buffer_ind-1)*T_sec+sort(TOA)];                
                TOA_other_whale=[TOA_other_whale (Buffer_ind-1)*T_sec+sort(TOA_other)];
                
                % DISPLAY THE NUMBER OF DETECTED SW

%                 Detected_SW = [Detected_SW length(Gather(Buffer_ind).ToA_separated)];
%           
%                 Gather_TOA=[Gather_TOA (Buffer_ind-1)*T_sec+sort(TOA)]; 
%                 Detected_SW = [Detected_SW nan(1,length(Gather_TOA)-length(Detected_SW))];
            end
            toc
          if ~isempty(Gather_TOA)
            if Tag_flag
                  Name_save=[Program_folder '/' 'Tag_detected_' filename '.xls'];
                  writematrix(Gather_TOA',Name_save);
                  filesave=[Program_folder '/' 'Other_whale_detected_' filename '.xls'];
                  writematrix(TOA_other_whale',Name_save);                                                    
            else                
                writematrix(Gather_TOA',Name_save);  % Save time of arrivals of the detected clicks
%                   result=table(Gather_TOA', Detected_SW','VariableNames',{'Detected Clicks','Detected SW'});
%                   writetable(result,filesave, 'WriteVariableNames', true);  % Save time of arrivals of the detected clicks                                                            
            end
          end
        end
        
    end
    




    




