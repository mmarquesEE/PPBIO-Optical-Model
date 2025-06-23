function data_dict = load_experimental_data(dir_path)
    % Reads all .xlsx files from a directory into a struct.
    data_dict = struct();
    % Get the list of files as a struct array
    file_list_struct = dir(fullfile(dir_path, '*.xlsx')); 
    
    % Get just the names into a cell array and sort them
    file_list_names = sort({file_list_struct.name});
    
    for i = 1:length(file_list_names)
        file_name = file_list_names{i}; % Get the filename string
        file_path = fullfile(dir_path, file_name);
        try
            data = readtable(file_path);
            if ~all(ismember({'Time', 'RefractiveIndex'}, data.Properties.VariableNames))
                warning('File %s is missing "Time" or "RefractiveIndex" and will be skipped.', file_name);
                continue;
            end
            field_name = matlab.lang.makeValidName(file_name);
            data_dict.(field_name) = data;
        catch ME
            warning('Could not read file %s. Error: %s', file_name, ME.message);
        end
    end
end