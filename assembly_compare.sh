#!/bin/bash

##Ok - to compare assemblies/align with nucmer and make dotplot?

for i in 1 2 9 12
do
    qseq=canu/${i}/${i}.nanopore.contigs.fasta
    
    ~/Code/mummer-4.0.0beta/nucmer -p wga/${i}.canu NC_016845.fasta ${qseq}
    ~/Code/mummer-4.0.0beta/mummerplot --png -p wga/${i}.layout.canu wga/${i}.canu.delta -l -R NC_016845.fasta -Q ${qseq}
    ~/Code/mummer-4.0.0beta/mummerplot --png -p wga/${i}.canu wga/${i}.canu.delta 

    ~/Code/mummer-4.0.0beta/mummerplot --postscript -p wga/${i}.layout.canu wga/${i}.canu.delta -l -R NC_016845.fasta -Q ${qseq}
    ~/Code/mummer-4.0.0beta/mummerplot --postscript -p wga/${i}.canu wga/${i}.canu.delta 

    lastal -P 8 /mithril/Data/NGS/Reference/kpneumo/kpneumo ${qseq} >wga_last/${i}.canu.maf
    last-dotplot wga_last/${i}.canu.maf wga_last/${i}.canu.pdf
    last-dotplot wga_last/${i}.canu.maf wga_last/${i}.canu.png

    
done

##spades

for i in 1 2 9 12
do
    qseq=spades/${i}/contigs.fasta
    
    ~/Code/mummer-4.0.0beta/nucmer -p wga/${i}.spades NC_016845.fasta ${qseq}
    ~/Code/mummer-4.0.0beta/mummerplot --png -p wga/${i}.layout.spades wga/${i}.spades.delta -l -R NC_016845.fasta -Q ${qseq}
    ~/Code/mummer-4.0.0beta/mummerplot --png -p wga/${i}.spades wga/${i}.spades.delta 

    ~/Code/mummer-4.0.0beta/mummerplot --postscript -p wga/${i}.layout.spades wga/${i}.spades.delta -l -R NC_016845.fasta -Q ${qseq}
    ~/Code/mummer-4.0.0beta/mummerplot --postscript -p wga/${i}.spades wga/${i}.spades.delta 

    lastal -P 8 /mithril/Data/NGS/Reference/kpneumo/kpneumo ${qseq} >wga_last/${i}.spades.maf
    last-dotplot wga_last/${i}.spades.maf wga_last/${i}.spades.pdf
    last-dotplot wga_last/${i}.spades.maf wga_last/${i}.spades.png


    
done
