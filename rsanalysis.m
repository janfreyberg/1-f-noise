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
    
    % Fix labels: if labels are A1-B32, replace with 1020
    preprocdata = utils.fixlabels(preprocdata);
    
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
    resampledata = ft_resampledata(cfgresample, preprocdata);
    
    
    
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
    else
        rejectartifactdata = getfield(load(cfgartifact.outputfile), 'data');
    end
    
    % Remove Short Trials (below 2 seconds)
    trials = 1:numel(rejectartifactdata.trial);
    keep = true(size(trials));
    for trial = trials
        keep(trial) = peak2peak(rejectartifactdata.time{trial}) > 4;
    end
    rejectartifactdata = ft_selectdata(rejectartifactdata, 'rpt', trials(keep));
    
    
    
    %% Frequency Analysis
    cfgfreq = [];
    cfgfreq.method = 'mtmfft';
    cfgfreq.output = 'pow';
    cfgfreq.channel = 1:64;
    cfgfreq.foi = 2:0.25:48;
    cfgfreq.tapsmofrq = 0.5;
%     cfgfreq.inputfile = cfgartifact.outputfile;
    cfgfreq.outputfile = fullfile(datadir, 'fft', strrep(files(subject).name, 'bdf', 'mat'));
    
    if ~exist(cfgfreq.outputfile, 'file')
        freqdata = ft_freqanalysis(cfgfreq, rejectartifactdata);
    else
        freqdata = getfield(load(cfgfreq.outputfile), 'freq');
    end
    
    
    %% Linear Regression
    % First, fit a line for each electrode
    frequencyband = freqdata.freq < 7 | freqdata.freq > 14;
    x = freqdata.freq(frequencyband)';
    X = [ones(size(x)), x];
    for electrode = 1:64
        y = log10(freqdata.powspctrm(electrode, frequencyband))';
        regression(subject).byelectrode(electrode, :) = X\y;
    end
    
    % Identify outliers based on intercept and +ve slope:
    outliers = regression(subject).byelectrode(:, 2) > mean(regression(subject).byelectrode(:, 2)) + 2*std(regression(subject).byelectrode(:, 2)) | ...
                regression(subject).byelectrode(:, 2) < mean(regression(subject).byelectrode(:, 2)) - 2*std(regression(subject).byelectrode(:, 2)) | ...
                regression(subject).byelectrode(:, 1) >= 0;
    
    % Fit line to meaned data:
    y = log10(mean(freqdata.powspctrm(~outliers, frequencyband)))';
    regression(subject).average = X\y;
end


% prepare summary data for spreadsheet
for subject = 1:n
    diagnosis{subject, 1} = files(subject).name(4:6);
    id{subject, 1} = files(subject).name(4:11);
    intercept(subject, 1) = regression(subject).average(1);
    slope(subject, 1) = regression(subject).average(2);
end

writetable(table(id, diagnosis, intercept, slope), fullfile(pwd, 'results.csv'))