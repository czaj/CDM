function repo = cdm_test_paths
% cdm_test_paths adds the CDM package (and Tools, if found alongside) to the path.

repo = fileparts(fileparts(mfilename('fullpath')));
addpath(repo);

toolsDir = fullfile(fileparts(repo),'Tools');
if exist(toolsDir,'dir') == 7
    addpath(genpath(toolsDir));
end
end
