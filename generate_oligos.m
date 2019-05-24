clear all
clc
% To download sequence of interest go to http://genome.ucsc.edu/cgi-bin/das/mm10/dna?segment=chrX:101640630,101742858)
% http://genome.ucsc.edu/cgi-bin/das/mm10/dna?segment=chrX:101742858,101760000
% Xist transcript
% https://www.ncbi.nlm.nih.gov/nuccore/NR_001463.3?report=fasta
% All transcripts
% https://www.gencodegenes.org/mouse/release_M10.html

geneSeq = fileread('Xist_transcript_short.txt'); % reference sequence from mouse
geneSeq = regexprep(geneSeq,'[^\w'']',''); % remove white spaces

oligoSize = 30;
Oligo = [];
mer = 1;
namePrint = 'Xist_transcript_short'; % name for the folder

mkdir(namePrint)
cd(namePrint)
fileID = fopen('Oligos_List.fa', 'w'); % name for the oligos list


while length(geneSeq) > mer + oligoSize
    % Calculate the GC content
    Oligo = geneSeq((1:oligoSize)+mer);
    num_G = nnz(Oligo == 'G');
    num_C = nnz(Oligo == 'C');
    gc_percentage = round( (num_G + num_C) / oligoSize * 100 );
    
    if gc_percentage < 80 && gc_percentage > 35
        % Remove the oligos with 7 consecutive nucleic acid that are the same
        A = 'AAAAAAA';
        G = 'GGGGGGG';
        C = 'CCCCCCC';
        T = 'TTTTTTT';
        RemoveCons = contains(string(Oligo), A) + contains(string(Oligo), G) + contains(string(Oligo), C) + contains(string(Oligo), T);
        
        if RemoveCons == 0
            fprintf(fileID, '> mer_%d \n %s \n', mer, Oligo); % save the filtered oligomers
        else
            Identities(mer) = 1000; % defined penalty for 7 consecutive nucleotides
        end
        
    else
        Identities(mer) = 1000; % defined penalty for too high or too low GC content
    end
    mer=1+mer;
end

fclose(fileID);
Identities(length(Identities)+1:mer-1) = 0;
save('Identities', 'Identities')

% Run on terminal for mouse
% makeblastdb -in RefGenome/GCF_000001635.26_GRCm38.p6_genomic.fna -out RefGenome/whole_mouse_genome.fa -dbtype nucl
% blastn -query /Users/anamota/Desktop/FRET_probe_design/Oligos1/Oligos_List.txt -db RefGenome/whole_mouse_genome.fa -out /Users/anamota/Desktop/FRET_probe_design/Oligos1/homologies.txt -evalue 10 -word_size 11 -ungapped -perc_identity 80 -qcov_hsp_perc 80 -num_threads 4
% awk '$1 ~ /^Query=/ {print $2} /^ Identities/ {print $3}' < /Users/anamota/Desktop/FRET_probe_design/test/homologies.txt > /Users/anamota/Desktop/FRET_probe_design/test/homologiesFilt.txt

% sed -i 's#60/60^$##g' homologiesFilt.txt
% cat *.txt > homologiesFiltAll.txt