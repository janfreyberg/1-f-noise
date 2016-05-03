%% File List
ft_defaults;
datadir = 'C:\Users\k1513504\Documents\resting-state-data\';
files = dir([datadir, '\*RS.bdf']);
n = numel(files);

%% Analysis Variables
trialdur = 60;
discard.start = 1;
discard.end = 1;


for subject = 1:n
    file = fullfile(datadir, files(subject).name);
    
    %% Trial Definition
    cfg_trldef = [];
    cfg_trldef.trialdef.eventtype = 'STATUS';
    cfg_trldef.trialfun = 'ft_trialfun_general';
    cfg_trldef.dataset = file;
    cfg_trldef.trialdef.prestim = -discard.start;
    cfg_trldef.trialdef.poststim = trialdur -discard.end;
    cfg_trldef.trialdef.eventvalue = 100:256;
    try
        cfg_trldef = ft_definetrial(cfg_trldef);
    catch define_trial_error
        disp(cfg_trldef.trialdef.eventvalue);
        disp(cfg_trldef.dataset);
        rethrow(define_trial_error);
    end
    cfg_trldef.trl = cfg_trldef.trl(1, :);
    for iTrial = 1:3
        cfg_trldef.trl(iTrial+1, :) = cfg_trldef.trl(1, :) +...
                                    [(iTrial)*60*cfg_trldef.trl(1, 3),...
                                    (iTrial)*60*cfg_trldef.trl(1, 3),...
                                    0, mod(iTrial, 2)];
    end
    
    %% Preprocessing
    % Preprocess data
    cfg_preproc = cfg_trldef;
    cfg_preproc.channel = 1:72;
    cfg_preproc.continuous = 'yes';
    cfg_preproc.detrend = 'yes';
    cfg_preproc.demean = 'yes';
    cfg_preproc.reref = 'no';
    
    preprocdata = ft_preprocessing(cfg_preproc);
    
    % Re-sample data for handling
    cfg_resample = [];
    cfg_resample.resamplefs = 512;
    cfg_resample.outputfile = fullfile(datadir, 'preproc-resample', strrep(files(subject).name, 'bdf', 'mat'));
    
    ft_resampledata(cfg_resample, preprocdata);
    
    
    %% Electrode rejection
    cfg_visual = [];
    cfg_visual.inputfile = fullfile(datadir, 'preproc-resample', strrep(files(subject).name, 'bdf', 'mat'));
    cfg_visual.method = 'trial';
    cfg_visual.keepchannel = 'no';
    
    reject_data = ft_rejectvisual(cfg_visual);
end