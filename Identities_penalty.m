clear all
clc

%% Open the files for homologies filtered
cd('Xist_transcript')
fileID = fopen('homologiesFilt.txt','r');
assert(fileID>0);

%% Sum scores in the file

load('Identities.mat');
tline = fgetl(fileID);
while ischar(tline)
    
    if strncmp(tline, 'mer_', 4) == 1 %when mer_ is encountered means that the score sum of a oligomer is over
        n = str2double(tline(5:end)); %new oligomer they are not ordered
        tline = fgetl(fileID); % first value is coming from the same chromosome
        if sum(tline(1:5) == '30/30') ~= 5
            Identities(n) = 1000; 
        end
    else
        if str2double(tline(1:2)) > 19
            Identities(n) = Identities(n) + 1;   
        end
    end
    
    tline = fgets(fileID); % Read line by line
end
fclose(fileID);

save('IdentitiesH', 'Identities')