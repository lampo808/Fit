suite = matlab.unittest.TestSuite.fromFolder('./test');
suite = suite(1:end-1);
w = warning();
warning('off', 'all');
result = run(suite);
warning(w)