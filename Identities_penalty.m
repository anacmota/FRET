clear all
clc

%% Open the files for homologies filtered
cd('Pls3')
fileID = fopen('homologiesFiltAll.txt','r');
assert(fileID>0);

%% Sum scores in the file

thresholdHomologies = 48; % 80% of the length of the oligomer: 48 if it's a 60-mer

load('Identities.mat');
tline = fgetl(fileID);
while ischar(tline)
    
    if strncmp(tline, 'mer_', 4) == 1 %when mer_ is encountered means that the score sum of a oligomer is over
        n = str2double(tline(5:end)); %new oligomer they are not ordered
        tline = fgetl(fileID); % first value is coming from the same chromosome
        if sum(tline(1:5) == '60/60') ~= 5
            Identities(n) = 1000; 
        end
    else
        if str2double(tline(1:2)) > thresholdHomologies
            Identities(n) = Identities(n) + 1;   
        end
    end
    
    tline = fgets(fileID); % Read line by line
end
fclose(fileID);

save('IdentitiesH', 'Identities')