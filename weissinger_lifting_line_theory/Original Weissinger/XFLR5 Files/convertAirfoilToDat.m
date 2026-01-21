function convertAirfoilToDat(inputFile, outputFile)
    % Function to read airfoil coordinate data from a text file
    % and write it to a .dat file in XFLR5-compatible format.

    % Open the input file
    fid_in = fopen(inputFile, 'r');
    if fid_in == -1
        error('Could not open input file: %s', inputFile);
    end

    % Read the first line (airfoil name)
    header = fgetl(fid_in);

    % Read remaining lines (coordinate data)
    coords = fscanf(fid_in, '%f %f', [2 Inf])';
    fclose(fid_in);

    % Open the output file
    fid_out = fopen(outputFile, 'w');
    if fid_out == -1
        error('Could not open output file: %s', outputFile);
    end

    % Write header
    % fprintf(fid_out, '%s\n', strtrim(header));

    % Write coordinate data
    for i = 1:size(coords, 1)
        fprintf(fid_out, '% .8f % .8f\n', coords(i,1), coords(i,2));
    end

    fclose(fid_out);
    fprintf('Successfully wrote to %s\n', outputFile);
end