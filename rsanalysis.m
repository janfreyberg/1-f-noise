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

for subject = 1:n
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
    cfgpreproc.outputfile = fullfile(datadir, 'preproc-resample', strrep(files(subject).name, 'bdf', 'mat'));
    if ~exist(cfgpreproc.outputfile, 'file')
        preprocdata = ft_preprocessing(cfgpreproc);
    else
        preprocdata = getfield(load(cfgpreproc.outputfile), 'data');
    end
    
    
    %% Create Bipolar Channels
    for trial = 1:numel(preprocdata.trial)
        preprocdata.trial{trial}(73, :) = ...
            preprocdata.trial{trial}(68, :) - ... %68 is ext. 4
            preprocdata.trial{trial}(67, :); %67 is ext. 3
        preprocdata.trial{trial}(74, :) = ...
            preprocdata.trial{trial}(69, :) - ... %69 is ext. 5
            preprocdata.trial{trial}(70, :); %70 is ext. 6
    end
    preprocdata.label{73} = 'EOG'; % Above and below the eye
    preprocdata.label{74} = 'JAW'; % Left and Right Mastoid
    
    
    %% Re-sample data for handling
    cfgresample = [];
    cfgresample.resamplefs = 512;
    cfgresample.outputfile = fullfile(datadir, 'preproc-resample', strrep(files(subject).name, 'bdf', 'mat'));
    if ~exist(cfgresample.outputfile, 'file')
        resampledata = ft_resampledata(cfgresample, preprocdata);
    end
    
    
    %% Artifact detection
    cfgartifact = [];
    cfgartifact.trl = cfgpreproc.trl;
    cfgartifact.continuous = 'yes';
    % EOG specific
    % This shouldn't reject many artifacts since we were using the
    % eyes-closed trials here, but just to be on the safe side.
    cfgartifact.artfctdef.eog.channel = 'EOG';
    cfgartifact.artfctdef.eog.trlpadding = 0;
    cfgartifact.artfctdef.eog.interactive = 'no';
    [cfgartifact, artifact.eog] = ft_artifact_eog(cfgartifact, resampledata);
    % Muscle specific
    % This should reject any artifacts that stem from head movements and
    % jaw clenching!
    cfgartifact.artfctdef.muscle.channel = 'JAW';
    cfgartifact.artfctdef.muscle.trlpadding = 0;
    cfgartifact.artfctdef.muscle.interactive = 'no';
    [cfgartifact, artifact.muscle] = ft_artifact_muscle(cfgartifact, resampledata);
    
    % Reject Artifacts
    cfgartifact.artfctdef.reject = 'partial';
    cfgartifact.outputfile = fullfile(datadir, 'autoreject-artifacts', strrep(files(subject).name, 'bdf', 'mat'));
    if ~exist(cfgartifact.outputfile, 'file')
        rejectartifactdata = ft_rejectartifact(cfgartifact, resampledata);
    end
    
    % Electrode rejection
%     cfgvisual = [];
%     cfgvisual.inputfile = fullfile(datadir, 'autoreject-artifacts', strrep(files(subject).name, 'bdf', 'mat'));
%     cfgvisual.outputfile = fullfile(datadir, 'visualreject-electrodes', strrep(files(subject).name, 'bdf', 'mat'));
%     cfgvisual.method = 'trial';
%     cfgvisual.metric = 'var';
%     cfgvisual.channel = 1:64;
%     cfgvisual.keepchannel = 'repair';
%     cfgvisual.keepchannel = 'no';
%     reject_data = ft_rejectvisual(cfgvisual);
    
    
    %% Frequency Analysis
    cfgfreq = [];
    cfgfreq.method = 'mtmfft';
    cfgfreq.output = 'pow';
    cfgfreq.channel = 1:64;
    cfgfreq.foilim = [2, 24];
    cfgfreq.tapsmofrq = 0.2;
    cfgfreq.inputfile = cfgartifact.outputfile;
    cfgfreq.outputfile = fullfile(datadir, 'fft', strrep(files(subject).name, 'bdf', 'mat'));
    
    if ~exist(cfgfreq.outputfile, 'file')
        freqdata = ft_freqanalysis(cfgfreq);
    end
    
    
    %% Linear Regression
end