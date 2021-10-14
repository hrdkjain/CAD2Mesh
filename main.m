%% Title: CAD2Mesh
% Author: Hardik Jain
% Date: 2019
% Link: https://github.com/hrdkjain/CAD2Mesh

%% Main source to convert CAD model to manifold mesh

if exist('../Matlab-Functions', 'dir') ~= 7
    fprintf("Error: Matlab-Functions not found, make sure it is placed in same directory as CAD2Mesh");
    return
end

addpath('../Matlab-Functions', '../Matlab-Functions/MeshUtils', '../Matlab-Functions/ProgressBar', 'Skeleton3D', 'Voxelization')
warning('off')
include  % header file include.m

% list files
inputFileList = {subdir(srcDir, srcExt).name};
inputFileListLength = length(inputFileList);
inputFileList = reshape(inputFileList, [inputFileListLength, 1]);
outputFileList = cell(inputFileListLength, 1);
outputFldList = {};
for i=1:inputFileListLength
    inputFile = inputFileList{i};
    outputFileList{i} = strrep(inputFile, srcDir, dstDir);
    outputFileDir = fileparts(outputFileList{i});
    if ~any(strcmp(outputFldList,outputFileDir))
        outputFldList{end+1,1} = outputFileDir;
    end
end

% create output directories
for i=1:length(outputFldList)
    [status, msg, msgID] = mkdir(outputFldList{i});
    if ~status
        print(msg);
    end
end
    
% Log file to hold the processed file names and time taken
fid = parallel.pool.Constant(@() fopen(fullfile(dstDir,'Report_CAD2Mesh.txt'),'a'),@fclose);
spmd(1)
    fprintf(fid.Value, "Started: %s\n", datestr(now));
end
hbar = ParforProgressbar(inputFileListLength, 'title', srcDir, 'parpool', {'local', numcores});
parfor i=1:inputFileListLength
    hbar.increment();
    inputFile = inputFileList{i};
    outputFile = outputFileList{i};
    if exist(outputFile,'file')==2
        continue;
    end
    
    fprintf('%s: ', inputFile);
    try
        [V,F] = read_mesh(inputFile);
    catch err
        printError(err);
        continue;
    end
    
    % voxelize and make genus zero
    [facen,vertn] = voxelelize_genus(V,F,sizen,1,1);
    vertn = perform_mesh_smoothing(facen,vertn,vertn);
    [vertn, facen] = withTS(vertn, facen);
    log = save_manifold_mesh(meshlabserverpath, vertn, facen, outputFile);
    if ~isempty(log)
        log = [ sprintf('%s:',outputFile), log];
        fprintf(fid.Value, '%s\n', log);
    end
    fprintf('\n');
end

spmd(1)
    fprintf(fid.Value, "Finished: %s\n", datestr(now));
end

close(hbar);
delete(gcp);
