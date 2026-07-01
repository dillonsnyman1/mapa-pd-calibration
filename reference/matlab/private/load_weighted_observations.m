function obs = load_weighted_observations(fixtures_dir)
    t = read_csv_table(fullfile(fixtures_dir, 'raw_observations_weighted.csv'));
    obs = table(double(t.score), double(t.bad), double(t.weight), ...
                'VariableNames', {'score', 'bad', 'weight'});
end
