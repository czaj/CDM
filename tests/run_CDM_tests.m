function nFail = run_CDM_tests
% run_CDM_tests runs the whole CDM test suite.
%
%   run_CDM_tests            prints a report and errors if anything failed
%   nFail = run_CDM_tests    returns the number of failures instead
%
% Checks tagged KNOWN refer to open issues from the code audit and are reported
% without failing the suite. When a KNOWN check starts passing, remove its entry.

cdm_test_paths;

tests = {@test_CDM_gradients,@test_CDM_options,@test_CDM_recovery};
nFail = 0;

for i = 1:numel(tests)
    name = func2str(tests{i});
    fprintf('\n=== %s ===\n',name);
    try
        nFail = nFail + tests{i}(1);
    catch ME
        nFail = nFail + 1;
        fprintf('  FAIL %s errored: %s\n',name,ME.message);
    end
end

fprintf('\n==================================================\n');
if nFail == 0
    fprintf('CDM test suite: all checks passed\n');
else
    fprintf('CDM test suite: %d failure(s)\n',nFail);
end

if nargout == 0 && nFail > 0
    error('run_CDM_tests:Failed','%d check(s) failed.',nFail);
end
end
