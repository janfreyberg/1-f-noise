function data = fixlabels(data)
% Fixes labels (1020 layout)
if any(ismember(data.label, 'B1'))
    data.label(1:72) = getfield(load(fullfile(pwd, '72labels.mat')), 'labels');
end
end
