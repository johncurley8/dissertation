% Define the number of firms
N = 250;

% Define file name
mod_filename = 'het_TFP.mod';
%mod_filename = './OccBin Solver/Jermann_Quadrini_2012_RBC_Occ.mod';

% Read the entire contents of the .mod file
fid = fopen(mod_filename, 'r');
mod_content = fread(fid, '*char')';
fclose(fid);

% Define the marker where the insertion should occur
marker = '%%%% MARKET CLEARING CONDITIONS (APPENDED BY mkt_clearing.m) %%%%';

% Find the position of the marker in the file
marker_position = strfind(mod_content, marker);

if isempty(marker_position)
    error('Marker for market clearing conditions not found in the .mod file.');
end

% Start with a properly formatted string
market_conditions = sprintf('\n[name=''Capital market clearing'']\n');
market_conditions = [market_conditions, 'k = '];

for i = 1:N
    if i < N
        market_conditions = [market_conditions, sprintf('k%d + ', i)];
    else
        market_conditions = [market_conditions, sprintf('k%d;\n\n', i)];
    end
end

% Labor Market Clearing
market_conditions = [market_conditions, sprintf('[name=''Labor market clearing'']\n')];
market_conditions = [market_conditions, 'n = '];

for i = 1:N
    if i < N
        market_conditions = [market_conditions, sprintf('n%d + ', i)];
    else
        market_conditions = [market_conditions, sprintf('n%d;\n\n', i)];
    end
end

% Bond Market Clearing
market_conditions = [market_conditions, sprintf('[name=''Bond market clearing'']\n')];
market_conditions = [market_conditions, 'b = '];

for i = 1:N
    if i < N
        market_conditions = [market_conditions, sprintf('b%d + ', i)];
    else
        market_conditions = [market_conditions, sprintf('b%d;\n\n', i)];
    end
end

% Dividend Market Clearing
market_conditions = [market_conditions, sprintf('[name=''Dividend market clearing'']\n')];
market_conditions = [market_conditions, 'd = '];

for i = 1:N
    if i < N
        market_conditions = [market_conditions, sprintf('d%d + ', i)];
    else
        market_conditions = [market_conditions, sprintf('d%d;\n\n', i)];
    end
end

% Output Market Clearing
market_conditions = [market_conditions, sprintf('[name=''Output market clearing'']\n')];
market_conditions = [market_conditions, 'y = '];

for i = 1:N
    if i < N
        market_conditions = [market_conditions, sprintf('y%d + ', i)];
    else
        market_conditions = [market_conditions, sprintf('y%d;\n\n', i)];
    end
end

% Ensure correct newline formatting
market_conditions = strrep(market_conditions, '\n', newline);

% Replace the marker with the updated market conditions
updated_mod_content = strrep(mod_content, marker, [marker newline market_conditions]);

% Write the updated content back to the .mod file
fid = fopen(mod_filename, 'w');
fwrite(fid, updated_mod_content);
fclose(fid);

fprintf('Market clearing conditions successfully formatted and inserted into %s for N = %d firms.\n', mod_filename, N);
