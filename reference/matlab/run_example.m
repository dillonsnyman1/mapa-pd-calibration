% RUN_EXAMPLE  Run the MAPA pipeline on the bundled example dataset.
%
%   cd reference/matlab
%   octave --no-gui run_example.m

if exist('OCTAVE_VERSION', 'builtin')
    try, pkg load datatypes; catch, end
end

this_dir     = fileparts(mfilename('fullpath'));
fixtures_dir = fullfile(this_dir, '..', 'fixtures');

obs = read_csv_table(fullfile(fixtures_dir, 'raw_observations.csv'));
obs = table(double(obs.score), double(obs.bad), 'VariableNames', {'score', 'bad'});

result = run_pipeline(obs, 10, 50, 10);

fprintf('\nCalibrated bands:\n');
fprintf('%-12s %-12s %-8s %-8s %-10s\n', 'score_min', 'score_max', 'n_obs', 'n_bads', 'pd');
fprintf('%s\n', repmat('-', 1, 54));
for i = 1:height(result.bands)
    b = result.bands(i, :);
    fprintf('%-12g %-12g %-8g %-8g %-10.4f\n', ...
        b.score_min, b.score_max, b.n_obs, b.n_bads, b.pd);
end

fprintf('\nPD at selected scores:\n');
test_scores = [400, 500, 600, 650, 700, 750, 800];
for i = 1:numel(test_scores)
    s = test_scores(i);
    fprintf('  score %d -> PD = %.4f\n', s, result.pd_for_score(s));
end
