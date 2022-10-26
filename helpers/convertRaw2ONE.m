dataPath = 'D:\ephys_results\processedData\audioVis';
rawFiles = dir(fullfile(dataPath,'*raw.mat'));
for d = 1:numel(rawFiles)
    dat = load(fullfile(rawFiles(d).folder,rawFiles(d).name));
    expRef = rawFiles(d).name(1:end-8);
    save2ONE(dataPath,expRef,dat)
end