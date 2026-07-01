function t = read_csv_table(filepath)
% Read a numeric CSV file into a table.
% Uses readtable when available (MATLAB), falls back to dlmread (Octave).
    if exist('readtable') > 0
        t = readtable(filepath);
    else
        fid = fopen(filepath, 'r');
        header_line = fgetl(fid);
        fclose(fid);
        col_names = strtrim(strsplit(strtrim(header_line), ','));
        data = dlmread(filepath, ',', 1, 0);
        col_cell = num2cell(data, 1);
        t = table(col_cell{:}, 'VariableNames', col_names);
    end
end
