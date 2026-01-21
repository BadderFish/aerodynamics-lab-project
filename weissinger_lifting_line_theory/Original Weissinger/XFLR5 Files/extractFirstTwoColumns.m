function extractFirstTwoColumns(inputFile, outputFile)
% EXTRACTFIRSTTWOCOLUMNS reads a text file with 4 columns and
% writes only the first 2 columns to a new text file.
%
% Usage:
%   extractFirstTwoColumns('input.txt', 'output.txt')

    % Read the full data from the input file
    data = readmatrix(inputFile);

    % Check if it has at least 2 columns
    if size(data, 2) < 2
        error('The input file must have at least 2 columns.');
    end

    % Extract the first two columns
    data2col = data(:, 1:2);

    % Open the output file for writing
    fid = fopen(outputFile, 'w');
    if fid == -1
        error('Could not open output file for writing.');
    end

    % Write header (optional, edit as needed)
    fprintf(fid, 'CD       T2-VLM1\n');

    % Write the two-column data
    for i = 1:size(data2col, 1)
        fprintf(fid, '%0.8f        %0.8f\n', data2col(i,1), data2col(i,2));
    end

    fclose(fid);
    fprintf('File "%s" created with first 2 columns from "%s".\n', outputFile, inputFile);
end