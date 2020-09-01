close all
clear
suite = matlab.unittest.TestSuite.fromFolder('./test');
fprintf('Number of tests found: %i\n', length(suite))
w = warning();
warning('off', 'all');
result = run(suite);
warning(w)