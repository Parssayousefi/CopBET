function tbl = CreateTBL(csvFolderPath)
    files = dir(fullfile(csvFolderPath, '*.csv'));
    numFiles = length(files);
    if numFiles == 0
        error('No CSV files found.');
    end
    
    % Prepare the table structure as expected by CopBET
    dataCells = cell(numFiles, 1);
    subjects = cell(numFiles, 1);
    tasks = cell(numFiles, 1);
    sessions = cell(numFiles, 1);
    
    for i = 1:numFiles
        filename = fullfile(csvFolderPath, files(i).name);
        csvData = readtable(filename);
        
        % Assuming the first few columns are time series data
        % and the last four columns are Subject, Session, Task, Dataset
        dataCells{i} = table2array(csvData(:, 1:end-4));
        subjects{i} = csvData.Subject(1);
        tasks{i} = csvData.Task(1);
        sessions{i} = csvData.Session(1);
    end
    
    % Create the final table
    tbl = table(dataCells, subjects, tasks, sessions, 'VariableNames', {'TimeSeries', 'Subject', 'Task', 'Session'});
end

