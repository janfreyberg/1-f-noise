%% File List
ft_defaults;
datadir = 'C:\Users\k1513504\Documents\resting-state-data\';
files = dir([datadir, '\*RS.bdf']);
n = numel(files);

% Group properties
groups(1).label = 'CTR';
groups(2).label = 'ASC';
for group = 1:numel(groups)
    groups(group).index = not(cellfun('isempty', strfind({files(:).name}, groups(group).label)));
    groups(group).n = sum(groups(group).index);
end

%% Analysis Variables
trialdur = 60;
discard.start = 1;
discard.end = 1;

nepoch = trialdur*4/2;


for subject = 1:n
    file = fullfile(datadir, files(subject).name);
    if groups(1).index(subject)
        group = 1;
    elseif groups(2).index(subject)
        group = 2;
    end
    %% Trial Definition
    cfgtrldef = [];
    cfgtrldef.trialdef.eventtype = 'STATUS';
    cfgtrldef.trialfun = 'ft_trialfun_general';
    cfgtrldef.dataset = file;
    cfgtrldef.trialdef.prestim = -discard.start;
    cfgtrldef.trialdef.poststim = trialdur -discard.end;
    cfgtrldef.trialdef.eventvalue = 100:256;
    
    if ~exist(fullfile(datadir, 'preproc-resample', strrep(files(subject).name, 'bdf', 'mat')), 'file')
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
    end
    
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
        preprocdata = utils.fixlabels(preprocdata);
        % Resample here
        cfgresample = [];
        cfgresample.resamplefs = 512;
        cfgresample.outputfile = cfgpreproc.outputfile;
        resampledata = ft_resampledata(cfgresample, preprocdata);
    else
        preprocdata = getfield(load(cfgpreproc.outputfile), 'data');
        resampledata = preprocdata;
    end
    clear preprocdata
    
    %% Create Bipolar Channels
    for trial = 1:numel(resampledata.trial)
        resampledata.trial{trial}(73, :) = ...
            resampledata.trial{trial}(68, :) - ... %68 is ext. 4
            resampledata.trial{trial}(67, :); %67 is ext. 3
        resampledata.trial{trial}(74, :) = ...
            resampledata.trial{trial}(69, :) - ... %69 is ext. 5
            resampledata.trial{trial}(70, :); %70 is ext. 6
    end
    resampledata.label{73} = 'EOG'; % Above and below the eye
    resampledata.label{74} = 'JAW'; % Left and Right Mastoid
    
    
    %% Artifact detection
    cfgartifact = [];
    cfgartifact.trl = resampledata.cfg.previous.trl;
        cfgartifact.trl(:, 1) = ceil(cfgartifact.trl(:, 1)/2);
        cfgartifact.trl(:, 3) = discard.start * cfgresample.resamplefs;
        cfgartifact.trl(:, 2) = cfgartifact.trl(:, 1) + cellfun(@numel, resampledata.time)' - 1;
    resampledata.cfg.previous = [];
    cfgartifact.continuous = 'yes';
    % EOG specific
    % This shouldn't reject many artifacts since we were using the
    % eyes-closed trials here, but just to be on the safe side.
    cfgartifact.artfctdef.eog.channel = 'EOG';
    cfgartifact.artfctdef.eog.trlpadding = 0;
    cfgartifact.artfctdef.eog.interactive = 'no';
    
    % Muscle specific
    % This should reject any artifacts that stem from head movements and
    % jaw clenching!
    cfgartifact.artfctdef.muscle.channel = 'JAW';
    cfgartifact.artfctdef.muscle.trlpadding = 0;
    cfgartifact.artfctdef.muscle.interactive = 'no';
    
    
    % Reject Artifacts
    cfgartifact.artfctdef.reject = 'partial';
    cfgartifact.outputfile = fullfile(datadir, 'autoreject-artifacts', strrep(files(subject).name, 'bdf', 'mat'));
    if ~exist(cfgartifact.outputfile, 'file')
        [cfgartifact, artifact.eog] = ft_artifact_eog(cfgartifact, resampledata);
        [cfgartifact, artifact.muscle] = ft_artifact_muscle(cfgartifact, resampledata);
        rejectartifactdata = ft_rejectartifact(cfgartifact, resampledata);
    else
        rejectartifactdata = getfield(load(cfgartifact.outputfile), 'data');
    end
    
    % Remove Short Trials (below 2 seconds)
    if ~exist(fullfile(datadir, 'fft', strrep(files(subject).name, 'bdf', 'mat')), 'file');
    trials = 1:numel(rejectartifactdata.trial);
    keep = true(size(trials));
    for trial = trials
        keep(trial) = peak2peak(rejectartifactdata.time{trial}) > 10;
    end
    rejectartifactdata = ft_selectdata(rejectartifactdata, 'rpt', trials(keep));
    end
    
    
    %% Frequency Analysis
    cfgfreq = [];
    cfgfreq.method = 'mtmfft';
    cfgfreq.output = 'pow';
    cfgfreq.channel = 1:64;
    cfgfreq.foi = 2:0.25:148;
    cfgfreq.foi((cfgfreq.foi>46&cfgfreq.foi<54) | (cfgfreq.foi>96&cfgfreq.foi<104)) = [];
    cfgfreq.tapsmofrq = 0.5;
    cfgfreq.outputfile = fullfile(datadir, 'fft', strrep(files(subject).name, 'bdf', 'mat'));
    
    if ~exist(cfgfreq.outputfile, 'file')
        freqdata = ft_freqanalysis(cfgfreq, rejectartifactdata);
    else
        freqdata = getfield(load(cfgfreq.outputfile), 'freq');
    end
    
    
    %% Linear Regression
    % First, fit a line for each electrode
    frequencyband = (freqdata.freq < 6 | freqdata.freq > 14) & freqdata.freq < 24;
%     frequencyband = (freqdata.freq > 14) & freqdata.freq < 24;
    hifrequencyband = freqdata.freq > 40;
    x = freqdata.freq(frequencyband)';
    X = [ones(size(x)), x];
    xhi = freqdata.freq(hifrequencyband)';
    Xhi = [ones(size(xhi)), xhi];
    for electrode = 1:64
        y = log10(freqdata.powspctrm(electrode, frequencyband))';
        regression(subject).byelectrode(electrode, :) = X\y;
    end
    
    % Identify outliers based on intercept and +ve slope:
    outliers(subject).electrodes =...
                regression(subject).byelectrode(:, 2) > mean(regression(subject).byelectrode(:, 2)) + 2*std(regression(subject).byelectrode(:, 2)) | ...
                regression(subject).byelectrode(:, 2) < mean(regression(subject).byelectrode(:, 2)) - 2*std(regression(subject).byelectrode(:, 2)) | ...
                regression(subject).byelectrode(:, 1) >= 0;
    
    % Fit line to meaned data:
    y = log10(mean(freqdata.powspctrm(~outliers(subject).electrodes, frequencyband)))';
    yhi = log10(mean(freqdata.powspctrm(~outliers(subject).electrodes, hifrequencyband)))';
    regression(subject).average = X\y;
    regression(subject).averagehi = Xhi\yhi;
    
    
    %% Data for plotting
    if subject == 1
        for i = 1:2
        groups(i).average.powspctrm = zeros(size(freqdata.freq));
        groups(i).average.freq = freqdata.freq;
        end
    end
    groups(group).average.powspctrm = groups(group).average.powspctrm + mean(freqdata.powspctrm(~outliers(subject).electrodes, :), 1)./groups(group).n;
end


%% Prepare summary data for spreadsheet
for subject = 1:n
    diagnosis{subject, 1} = files(subject).name(4:6);
    id{subject, 1} = files(subject).name(4:11);
    intercept(subject, 1) = regression(subject).average(1);
    slope(subject, 1) = regression(subject).average(2);
    intercept_high(subject, 1) = regression(subject).averagehi(1);
    slope_high(subject, 1) = regression(subject).averagehi(2);
end

writetable(table(id, diagnosis, intercept, slope, intercept_high, slope_high), fullfile(pwd, 'results.csv'));


%% Prepare a plot of the average spectrum
figure; hold on;
colors = {[83 148 255]/255, [255 117 117]/255};
for group = 1:2
    plot(groups(group).average.freq, log10(groups(group).average.powspctrm), 'Color', colors{group});
    ShadePlotForEmpahsis([min(freqdata.freq(frequencyband)), max(freqdata.freq(frequencyband))], 'k', 0.2);
    ShadePlotForEmpahsis([min(freqdata.freq(hifrequencyband)), max(freqdata.freq(hifrequencyband))], 'k', 0.2);
end

