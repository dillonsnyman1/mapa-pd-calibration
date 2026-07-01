function obs = load_observations(fixtures_dir)
    t = read_csv_table(fullfile(fixtures_dir, 'raw_observations.csv'));
    obs = table(double(t.score), double(t.bad), 'VariableNames', {'score', 'bad'});
end
