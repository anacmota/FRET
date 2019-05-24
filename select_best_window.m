clear all
close all
clc

%% Open the mouse barcodes file

barcodes = fileread('mouse barcodes.txt');
barcodes = char(strsplit(barcodes,'\n'));

%% Variables

spaceO = input('Space in between the oligomers same color: (in case of alternating: 0)\n');
spaceC = input('Space in between the oligomers different colors:\n');
numberColor = input('Number of oligomers alternating between different colors:\n');
numberOligos = input('Number of oligomers as output:\n');

namePrint = ['Ssmall' num2str(spaceO) '_Sbig' num2str(spaceC) '_OinI' num2str(numberColor) '_Ototal' num2str(numberOligos)];

cd('Xist_transcript')
oligoSize = 30;
TOLERANCE = 0.3;
OligosIslands = ceil(numberColor*0.7); % 4 oligos minimum for 5 and 8 for 10
homLimit = 5;
load('IdentitiesH.mat')

%% Open the files for oligo list

OligoListID = fileread('Oligos_List.fa');
OligoListID = strsplit(OligoListID,' ');
mer = char(OligoListID(2:4:end));
OligoList(str2double(mer(:, 5:end))) = OligoListID(4:4:end);

%% Create all possible windows

mkdir(namePrint)
cd(namePrint)
stepBig = spaceC + oligoSize;
VarBig = floor(spaceC*TOLERANCE);
stepSmall = spaceO + oligoSize;
VarSmall = ceil(spaceO*TOLERANCE);


position = ones(3000, numberOligos)*99999;
homology = ones(3000, numberOligos)*99999;
a = 1;
current = a*VarBig;

while (current + (stepBig + (stepSmall+VarSmall)*(numberColor-1))*numberOligos/numberColor) + VarBig*7 < length(Identities) % removed the +VarBig
    i = 1;
    while i <= numberOligos
        
        b = 1;
        for currentTemp = current-VarBig:current+VarBig % Variability from the large steps
            
            for c = 1:numberColor % The small steps are according with nr oligos inside the islands
                [ScoreTemp, IndexTemp] = min(Identities(currentTemp - VarSmall + stepSmall : currentTemp + VarSmall + stepSmall));
                currentTemp = currentTemp-VarSmall+stepSmall+IndexTemp-1;
                positionTemp(b, c) = currentTemp;
                homologyTemp(b, c) = ScoreTemp;
            end
            b = b + 1;
        end
         
        I = [1:floor(length(homologyTemp)/2), ceil(length(homologyTemp)/2):-1:1]/length(homologyTemp);
        [~, currentTemp] = max((sum(homologyTemp, 2)==min(sum(homologyTemp, 2)))+I');
        
        homology(a, i:i+numberColor-1) = homologyTemp(currentTemp, :);
        position(a, i:i+numberColor-1) = positionTemp(currentTemp, :);
        current = position(a, i+numberColor-1) + stepBig - stepSmall; % the small step is included even in the first oligo which should not have
        i = i + numberColor;
    end
    a = a + 1;
    current = a*2*VarBig;
end

%% Select the best window for oligomers by score, perfect oligos or oligos matches below a threshold

[Score, Index] = min( sum( homology, 2));

[BadNaN, IndexNaN] = min( sum( homology~=0, 2));

[BadN, IndexN] = min( sum( homology > homLimit, 2));

%% Plot the good and bad oligomers

figure(1)
sgtitle('Oligos below 80% homology and with 35%-80% GC content')
subplot(3, 1, 1)
Bad = sum(homology(Index, :)~=0);
area(homology(Index, :)==0, 'LineStyle', 'none')
hold on
area(homology(Index, :)~=0, 'LineStyle', 'none')
axis([1 numberOligos 0.49 0.51])
set(gca,'YTick', [])
xlabel('Oligomer number')
legend(['good  n: ' num2str(numberOligos-Bad)], ['bad    n: ' num2str(Bad)])
title(['minimum score (homologies and GC penalties)    window: ' num2str(Index)])

subplot(3, 1, 2)
area(homology(IndexNaN, :)==0, 'LineStyle', 'none')
hold on
area(homology(IndexNaN, :)~=0, 'LineStyle', 'none')
axis([1 numberOligos 0.49 0.51])
set(gca,'YTick', [])
xlabel('Oligomer number')
legend(['good  n: ' num2str(numberOligos-BadNaN)], ['bad    n: ' num2str(BadNaN)])
title(['min bad oligos with no homologies    window: ' num2str(IndexNaN)])

subplot(3, 1, 3)
area(homology(IndexN, :)<=homLimit, 'LineStyle', 'none')
hold on
area(homology(IndexN, :)>homLimit, 'LineStyle', 'none')
axis([1 numberOligos 0.49 0.51])
set(gca,'YTick', [])
xlabel('Oligomer number')
legend(['good  n: ' num2str(numberOligos-BadN)], ['bad    n: ' num2str(BadN)])
title(['min bad oligos with max ' num2str(homLimit) ' homologies    window: ' num2str(IndexN)])

saveas(gcf, 'Oligos fit distribution.png')

%% Select the best window for islands

Islands = reshape( homology', numberColor, [])';
Islands = reshape( sum( Islands <= homLimit, 2), numberOligos/numberColor, [])'; % 1 is good and 0 is bad

Islands = Islands';
BadIslandsPosition = find(Islands < OligosIslands);
InBetweenBad = find(BadIslandsPosition( 2:length(BadIslandsPosition)) - BadIslandsPosition( 1:(length(BadIslandsPosition)-1))==2);
Islands(BadIslandsPosition(InBetweenBad)+1) = 0;
Islands = Islands';

[BadIslands, IndexIslands] = min( sum( Islands < OligosIslands, 2)); % 5 is all good and 0 is all bad

GoodIslandsPos = repmat(Islands(IndexIslands, :)>=OligosIslands, numberColor, 1);
FinalOligos = position(IndexIslands, GoodIslandsPos(:));
FinalOligos(homology(IndexIslands, GoodIslandsPos(:))>homLimit) = NaN;

%% Plot the good and bad islands

figure(2)
area(Islands(IndexIslands, :)>=OligosIslands, 'LineStyle', 'none')
hold on
area(Islands(IndexIslands, :)<OligosIslands, 'LineStyle', 'none')
axis([1 numberOligos/numberColor 0.49 0.51])
set(gca,'YTick', [])
xlabel('Island number')
legend(['good  i: ' num2str(numberOligos/numberColor-BadIslands) '  o: ' ...
    num2str(sum(homology(IndexIslands, GoodIslandsPos(:))<=homLimit))], ['bad    i: ' num2str(sum(BadIslands)) ...
    '  o: ' num2str(numberOligos-sum(homology(IndexIslands, GoodIslandsPos(:))<=homLimit))])
title(['min bad oligos with max ' num2str(homLimit) ' homologies    window: ' num2str(IndexIslands)])
saveas(gcf, 'Islands fit distribution.png')

%% Plot the homologies and identities from the final oligos

figure(3)
hist(position(IndexIslands, 2:numberOligos)-position(IndexIslands, 1:(numberOligos-1))-oligoSize, 30)
title(['total length ' num2str(position(IndexIslands, numberOligos)-position(IndexIslands, 1)) ' bp'])
xlabel('Distance between consecutive oligos (bp)')
ylabel('# Oligos')
saveas(gcf, 'Oligos distance.png')

figure(4)
hist(Identities(FinalOligos(~isnan(FinalOligos))), 30)
title('Homologies distribution in the selected probe')
xlabel('Homologies per oligo')
ylabel('# Oligos')
saveas(gcf, 'Homologies per oligo.png')

%% Insert barcodes and print the final oligos

RWdonor = 1; 
RWacceptor = 2; 
FW = 20; % From 8 to 26 since from 3 to 7 was used for the spotting probes as FW 8-12 13-17 18-20
i = 1;
fileID = fopen('result.txt', 'w');
for elec = 1:length(FinalOligos)/numberColor
        
    if rem( elec, 2) == 0 
        for n = 1:numberColor
            
            if isnan(FinalOligos(i))==0
                fprintf(fileID, ['> mer_' num2str(i) '\n' lower(barcodes(FWdonor, 1:20)) OligoList{FinalOligos(i)} ...
                lower(barcodes(RW, 1:20)) '\n']);
            end
            i = i + 1;
        end
    
    elseif rem( elec, 2) == 1
        for n = 1:numberColor
            if isnan(FinalOligos(i))==0
                fprintf(fileID, ['> mer_' num2str(i) '\n' lower(barcodes(FWacceptor, 1:20)) OligoList{FinalOligos(i)} ...
                lower(barcodes(RW, 1:20)) '\n']);
            end
            i = i + 1;
        end
    end
    
end
fclose(fileID);


fileID = fopen('Oligos_List_information.txt', 'w');
% just integrate in the list the oligomers and removal of spaces
fprintf(fileID, [' Oligo Size: %d \n Space in between the oligomers same color: %d \n Space in between the' ...
' oligomers different colors: %d \n Number of oligomers alternating between different colors: %d \n ' ...
'Number of oligomers in a probe: %d \n Homology limit accepted: %d \n Ratio of tolerance given to the distances: %d \n' ...
' Minimum number of oligos per island: %d \n'], oligoSize, spaceO, spaceC, numberColor, numberOligos, homLimit, TOLERANCE, OligosIslands);
fclose(fileID);

% melt_duplex /Users/anamota/Desktop/FRET_probe_design/result.txt -C -f 25 -n 0.3
% tm = readtable('tm.txt');
% Tm = tm{:,{'Tm'}};
% hist(Tm, 30)
% xlabel('Melting temperature (°C)')
% ylabel('Counts #')