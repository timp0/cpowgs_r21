#!/bin/bash

assfile=$1

genome=$2

case $genome in
    kpn )
	echo "kpn"
	reffasta=/mithril/Data/NGS/Reference/kpneumo/NC_016845.fasta
	;;
esac

outdir=$3/wga

prefix=$4

mkdir -p ${outdir}


##Ok - to compare assemblies/align with nucmer and make dotplot?

~/Code/mummer-4.0.0beta/nucmer -p ${outdir}/${prefix} ${reffasta} ${assfile}
~/Code/mummer-4.0.0beta/mummerplot --png -p ${outdir}/${prefix}.layout.canu ${outdir}/${prefix}.delta \
				   -l -R ${reffasta} -Q ${assfile}
~/Code/mummer-4.0.0beta/mummerplot --png -p ${outdir}/${prefix}.canu ${outdir}/${prefix}.delta 

~/Code/mummer-4.0.0beta/mummerplot --postscript -p ${outdir}/${prefix}.layout.canu ${outdir}/${prefix}.delta \
				   -l -R ${reffasta} -Q ${assfile}
~/Code/mummer-4.0.0beta/mummerplot --postscript -p ${outdir}/${prefix}.canu ${outdir}/${prefix}.delta 

lastal -P 8 /mithril/Data/NGS/Reference/kpneumo/kpneumo ${assfile} >${outdir}/${prefix}.canu.maf
last-dotplot ${outdir}/${prefix}.canu.maf ${outdir}/${prefix}.canu.pdf
last-dotplot ${outdir}/${prefix}.canu.maf ${outdir}/${prefix}.canu.png

    
