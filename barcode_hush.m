clear all
clc
close all

fileID = fopen('Pls3.hush9.fa','r');
assert(fileID>0);

tline = fgetl(fileID);
while ischar(tline)
    n = str2double(tline(6:end)); %new oligomer they are not ordered
    tline = fgetl(fileID); % first value is coming from the same chromosome
    tline = split(tline, ', ');
    sequence(n, :) = tline{1};
    hits(n) = str2num(tline{2});
    
    tline = fgets(fileID); % Read line by line
end
fclose(fileID);


% Off = sum(hits > 20)
% position = find(hits > 20)
% NumberOligos = length(hits)-Off

%% Open the mouse barcodes file

barcodes = fileread('../mouse barcodes.txt');
barcodes = char(strsplit(barcodes,'\n'));

%% Insert barcodes and print the final oligos

RWdonor = 1; 
RWacceptor = 2; 
FW = 11; % 8 Ogt  9 Magix   10 Kdm5c   11 Pls3
i = 1; % to not count the removed oligos for high mismatches in order to define the dye
fileID = fopen('all_final.txt', 'a');
for each = 1:n
    if hits(each) < 20
        if rem(i, 2) == 0
            fprintf(fileID, [lower(barcodes(FW, 1:20)) lower(sequence(each, :)) ...
                lower(barcodes(RWdonor, 1:20)) '\n']);
            i = i + 1;
        else
            fprintf(fileID, [lower(barcodes(FW, 1:20)) lower(sequence(each, :)) ...
                lower(barcodes(RWacceptor, 1:20)) '\n']);
            i = i + 1;
        end
    end
end
fclose(fileID);