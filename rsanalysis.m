%% File List
ft_defaults;
datadir = 'C:\Users\k1513504\Documents\resting-state-data\';
files = dir([datadir, '\*RS.bdf']);
n = numel(files);

%% Analysis Variables
trialdur = 60;
discard.start = 1;
discard.end = 1;

nepoch = trialdur*4/2;

for subject = 1:1
    file = fullfile(datadir, files(subject).name);
    
    %% Trial Definition
    cfgtrldef = [];
    cfgtrldef.trialdef.eventtype = 'STATUS';
    cfgtrldef.trialfun = 'ft_trialfun_general';
    cfgtrldef.dataset = file;
    cfgtrldef.trialdef.prestim = -discard.start;
    cfgtrldef.trialdef.poststim = trialdur -discard.end;
    cfgtrldef.trialdef.eventvalue = 100:256;
    try
        cfgtrldef = ft_definetrial(cfgtrldef);
    catch define_trial_error
        disp(cfgtrldef.trialdef.eventvalue);
        disp(cfgtrldef.dataset);
        rethrow(define_trial_error);
    end
    cfgtrldef.trl = cfgtrldef.trl(1, :);
    % assign trials based on trigger
    for iTrial = 1:3
        cfgtrldef.trl(iTrial+1, :) = cfgtrldef.trl(1, :) +...
                                    [(iTrial)*60*cfgtrldef.trl(1, 3),...
                                    (iTrial)*60*cfgtrldef.trl(1, 3),...
                                    0, mod(iTrial, 2)];
    end
    % only take the second and fourth trial (eyes closed)
    cfgtrldef.trl = cfgtrldef.trl([2, 4], :);
    
    %% Preprocessing
    % Preprocess data
    cfgpreproc = cfgtrldef;
    cfgpreproc.channel = 1:72;
    cfgpreproc.continuous = 'yes';
    cfgpreproc.detrend = 'yes';
    cfgpreproc.demean = 'yes';
    cfgpreproc.reref = 'no';
    preprocdata = ft_preprocessing(cfgpreproc);
    
    %% Create Bipolar Channels
    preprocdata.label{73} = 'EOG';
    preprocdata.label{74} = 'JAW';
    for trial = 1:numel(preprocdata.trial)
        preprocdata.trial{trial}(73, :) = ...
            preprocdata.trial{trial}(ismember(preprocdata.label, 'EXG4'), :) - ...
            preprocdata.trial{trial}(ismember(preprocdata.label, 'EXG3'), :);
        preprocdata.trial{trial}(74, :) = ...
            preprocdata.trial{trial}(ismember(preprocdata.label, 'EXG5'), :) - ...
            preprocdata.trial{trial}(ismember(preprocdata.label, 'EXG6'), :);
    end
    
    %% Re-sample data for handling
    cfgresample = [];
    cfgresample.resamplefs = 512;
    cfgresample.outputfile = fullfile(datadir, 'preproc-resample', strrep(files(subject).name, 'bdf', 'mat'));
    resampledata = ft_resampledata(cfgresample, preprocdata);
    
    %% Artifact detection
    cfgartifact = [];
    cfgartifact.trl = cfgpreproc.trl;
    cfgartifact.continuous = 'yes';
    % EOG specific
    cfgartifact.artfctdef.eog.channel = 'EOG';
    cfgartifact.artfctdef.eog.trlpadding = 0;
    cfgartifact.artfctdef.eog.interactive = 'no';
    [cfgartifact, artifact.eog] = ft_artifact_eog(cfgartifact, resampledata);
    % Muscle specific
    cfgartifact.artfctdef.muscle.channel = 'JAW';
    cfgartifact.artfctdef.muscle.trlpadding = 0;
    cfgartifact.artfctdef.muscle.interactive = 'yes';
    [cfgartifact, artifact.muscle] = ft_artifact_muscle(cfgartifact, resampledata);
    
    % Reject Artifacts
    
    
    %% Electrode rejection
%     cfgvisual = [];
%     cfgvisual.inputfile = fullfile(datadir, 'preproc-resample', strrep(files(subject).name, 'bdf', 'mat'));
%     cfgvisual.method = 'trial';
%     cfgvisual.keepchannel = 'no';
%     reject_data = ft_rejectvisual(cfgvisual);
    
    
    
    
    %% Chop everything into 2s trials, reject blinks
    
    
    
end