function mkdir_CB(dirName)
    
    if ~exist(dirName,'dir')
        mkdir(dirName)
    end