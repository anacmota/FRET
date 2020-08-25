%% Correction values read from files
clear all
close all
clc

Correction = input('Correct images? (y/n)\n', 's');

filenum = '436'; % USER INPUT!!! 387 436
ID = readtable(['iAM', filenum,'.csv']);

% 432 is Magix
% 436 is Pls3
% 443 is Ogt
% 444 is Kdm5c

color = input('Which dyes? (C for Cy3-Cy5, A for a488-a594 or D for deconvolution of a488-a594)\n', 's');
if strcmpi(color, 'C')
    color={'Cy5'; 'tmr'; 'tmrfret'};
elseif strcmpi(color, 'A')
    color={'a594'; 'a488'; 'a594fret'};
elseif strcmpi(color, 'D')
    color={'dw_a594'; 'dw_a488'; 'dw_a594fret'};
end

InactInt = [];
ActInt = [];
SplitInactFRET = [];
SplitActFRET = [];
SplitInactDonor = [];
SplitActDonor = [];
I = [];
A = [];
NSNR_F_I = [];
NSNR_F_A = [];

% the label 0 are indistinguishable clusters
Position0 = find(ID.Label == 0);
ID(Position0, :) = [];

% Reading the number of the last field from the long extension
LastIndex = strsplit(ID.File{length(ID.File)},'/');
Extension = strjoin(LastIndex(1:length(LastIndex)-1), '/');
LastIndex = strsplit(LastIndex{length(LastIndex)},'.');

for Fields = 1:str2double(LastIndex{1})
    
    if Fields < 10
        extension = [Extension, '/00', num2str(Fields), '.NM'];
    elseif Fields > 9 && Fields < 100
        extension = [Extension, '/0', num2str(Fields), '.NM'];
    else
        extension = [Extension, num2str(Fields), '.NM'];
    end
    
    PositionFile =  find(strcmp(extension, ID.File));
    
    %% Identify the Nuclei to be analysed
    nuclei = unique(ID.Nuclei(PositionFile));
    
    for n = 1:length(nuclei)
        PositionNuclei = PositionFile(eq(nuclei(n), ID.Nuclei(PositionFile)));

        if (sum(ID.Label(PositionNuclei(strcmp(ID.Channel(PositionNuclei), color(1))))==1) == sum(ID.Label(PositionNuclei(strcmp(ID.Channel(PositionNuclei), color(2))))==1)) ...
                && (sum(ID.Label(PositionNuclei(strcmp(ID.Channel(PositionNuclei), color(1))))==2) == sum(ID.Label(PositionNuclei(strcmp(ID.Channel(PositionNuclei), color(2))))==2)) ...
                && any(strcmp(ID.Channel(PositionNuclei), color(2))) && all(ID.Value(PositionNuclei) < 10^4)
            
            Inact = PositionNuclei(ID.Label(PositionNuclei) == 2);
            Act = PositionNuclei(ID.Label(PositionNuclei) == 1);

            if length(Inact) == 2 && length(Act) == 3
                InactInt = [InactInt; [ID.Value(Inact)' 0]];
                ActInt = [ActInt; ID.Value(Act)'];
                
            elseif length(Act) == 2 && length(Inact) == 3
                InactInt = [InactInt; ID.Value(Inact)'];
                ActInt = [ActInt; [ID.Value(Act)' 0]];
                
            elseif length(Inact) == 3 && length(Act) == 3
                InactInt = [InactInt; ID.Value(Inact)'];
                NSNR_F_I = [NSNR_F_I; ID.NSNR(Inact(strcmp(ID.Channel(Inact), color(3))))'];
                ActInt = [ActInt; ID.Value(Act)'];
                NSNR_F_A = [NSNR_F_A; ID.NSNR(Act(strcmp(ID.Channel(Inact), color(3))))'];
                
            elseif length(Inact) >= 3 && length(Act) >= 3
        
                [Var1, Var2, Var3, Var4, Var5] = splitDotsSort(ID, Inact, color);
                NSNR_F_I = [NSNR_F_I; Var5];
                InactInt = [InactInt; Var1]; % [A, D, F] Dot chosen based on Max FRET intensity from split dots
                SplitInactFRET = [SplitInactFRET; Var2]; % FRET intensity from all split dots
                SplitInactDonor = [SplitInactDonor; Var3]; % Donor intensity from all split dots
                I = [I; Var4]; % If there is one dot (not splitted), Var4 is 0. If there is more than 1 dot, Var4 is 1.
                
                [Var1, Var2, Var3, Var4, Var5] = splitDotsSort(ID, Act, color);
                NSNR_F_A = [NSNR_F_A; Var5];
                ActInt = [ActInt; Var1];
                SplitActFRET = [SplitActFRET; Var2];
                SplitActDonor = [SplitActDonor; Var3];
                A = [A; Var4];
            end
            
        elseif any(strcmp(ID.Channel(PositionNuclei), color(2))) && all(ID.Value(PositionNuclei) < 10^4)
            disp(['missing field ' num2str(Fields) ' nuclei ' num2str(nuclei(n))])
        end
        
    end
end

% Correction of channels

if strcmpi(Correction, 'y')
    binrng = -0.5:0.1:1;             % Create Bin Ranges
    if isempty(setdiff(color, {'a594'; 'a488'; 'a594fret'}))
        ActIntCorr = ActInt(:, 3) - ActInt(:, 2)*0.085 - ActInt(:, 1)*0.056;
        InactIntCorr = InactInt(:, 3) - InactInt(:, 2)*0.085 - InactInt(:, 1)*0.056;
    elseif isempty(setdiff(color, {'dw_a594'; 'dw_a488'; 'dw_a594fret'}))
        ActIntCorr = ActInt(:, 3) - ActInt(:, 2)*0.107 - ActInt(:, 1)*0.10713;
        InactIntCorr = InactInt(:, 3) - InactInt(:, 2)*0.107 - InactInt(:, 1)*0.10713;
    elseif isempty(setdiff(color, {'Cy5'; 'tmr'; 'tmrfret'}))
        ActIntCorr = ActInt(:, 3) - ActInt(:, 2)*0.278 - ActInt(:, 1)*0.09;
        InactIntCorr = InactInt(:, 3) - InactInt(:, 2)*0.278 - InactInt(:, 1)*0.09;
    end
else
    ActIntCorr = ActInt(:, 3);
    InactIntCorr = InactInt(:, 3);
    binrng = 0:0.05:1;
end

% Calculate FRET ratio

InactEfficiencyFRET = InactIntCorr./(InactIntCorr+InactInt(:, 2));
ActEfficiencyFRET = ActIntCorr./(ActIntCorr+ActInt(:, 2));

length(InactInt(:, 3))
length(ActInt(:, 3))





%% Split dots measurements

if ~isempty(SplitInactFRET)
    InactEffFRET = SplitInactFRET./(SplitInactFRET+SplitInactDonor);
    ActEffFRET = SplitActFRET./(SplitActFRET+SplitActDonor);
    
    figure('Position', [0 0 1200 600])
    sgtitle('Dot splitting occurence and distribution')
    
    subplot(2, 4, 1)
    I = logical(I);
    barplot(InactEffFRET(I), 'b', binrng)
    xlabel('Inactive allele FRET efficiency (split dot)')
    
    subplot(2, 4, 2)
    barplot(InactEffFRET(~I), 'b', binrng)
    xlabel('Inactive allele FRET efficiency (not split dot)')
    
    subplot(2, 4, 3)
    A = logical(A);
    barplot(ActEffFRET(~A), 'y', binrng)
    xlabel('Active allele FRET efficiency (not split dot)')
    %[h,p] = ttest2(InactEffFRET(I),ActEffFRET(~A));
    
    subplot(2, 4, 4)
    barplot(ActEffFRET(A), 'y', binrng)
    xlabel('Active allele FRET efficiency (split dot)')
    
    subplot(2, 4, 5)
    binrngD = 0:2:20;
    barplot(SplitInactDonor(I), 'b', binrngD)
    xlabel('Inactive allele donor intensity (split dot)')
    
    subplot(2, 4, 6)
    barplot(SplitInactDonor(~I), 'b', binrngD)
    xlabel('Inactive allele donor intensity (not split dot)')
    
    subplot(2, 4, 7)
    barplot(SplitActDonor(~A), 'y', binrngD)
    xlabel('Active allele donor intensity (not split dot)')
    
    subplot(2, 4, 8)
    barplot(SplitActDonor(A), 'y', binrngD)
    xlabel('Active allele donor intensity (split dot)')
end






function [Var1, Var2, Var3, Var4, Var5] = splitDotsSort(ID, Position, color)
% InactInt = [InactInt; Var1];               % [A, D, F] Dot chosen based on Max FRET intensity from split dots
% SplitInactFRET = [SplitInactFRET; Var2];   % FRET intensity from all split dots
% SplitInactDonor = [SplitInactDonor; Var3]; % Donor intensity from all split dots
% I = [I; Var4];                             % If there is one dot (not splitted), Var4 is 0. If there is more than 1 dot, Var4 is 1.

positionXYZ = [ID.x*0.13, ID.y*0.13, ID.z*0.3];

[~, i] = sort(sum(positionXYZ(Position, :), 2));
Position = Position(i);

Acceptor = ID.Value(Position(strcmp(ID.Channel(Position), color(1))));
Donor = ID.Value(Position(strcmp(ID.Channel(Position), color(2))));
FRET_Position = Position(strcmp(ID.Channel(Position), color(3)));
FRET = ID.Value(FRET_Position);

if length(Position) == 3
    Var1 = [Acceptor, Donor, FRET];
    Var4 = false;
    Var5 = ID.NSNR(FRET_Position);
elseif length(Position) > 3
    [~,i] = max(FRET);
    Var1 = [Acceptor(i), Donor(i), FRET(i)];
    Var4 = true(length(FRET), 1);
    Var5 = ID.NSNR(FRET_Position(i));
end

FRET = FRET - Donor*0.107 - Acceptor*0.10713; % This step might be automatised with the table
Var2 = FRET;
Var3 = Donor;
end

function barplot(Var, color, binrng)
counts = histc(Var, binrng);
bar(binrng, counts, color)
ylabel('# nuclei')
title(['median = ' num2str(median(Var))])
end