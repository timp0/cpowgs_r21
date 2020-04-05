##Ok - this code is to do a quick annotation of k-pneumo like assemblies output from canu
##Example use
##Rscript ~/Code/timp_nanopore/oxford/fq_check.R $PWD/170519_cdna_cel_smu1_111.csv.gz ~/Dropbox/Data/Nanopore/170528_worm/ 170519_cdna

library(tidyverse)

codedir="~/Code/timp_nanopore/"


##Default input variable
workdir="/atium/Data/Nanopore/Analysis/170512_2897_KLPN/"
assfile=file.path(workdir, "canu_assembly", "170512_2897_KLPN.contigs.fasta")
prefix="170512_2897_KLPN"
outdir="~/Dropbox/Data/CPOWGS/170512_2897_KLPN"


##This pulls in the command line stuff
args=(commandArgs(TRUE))

##Check if there are any variables, and if not, keep default
if(length(args)==0) {
    print("No arguments.")  
}else{
    csvfile=args[1]
    plotdir=args[2]
    modif=args[3]
}

##print(args)

##make annotation dir
if (!(dir.exists(file.path(outdir, "annot")))) {
    dir.create(file.path(outdir, "annot"))
}


##Setup blast/card info
blast.cols=c("query.id", "subject.id", "per.identity", "align.len", "mismatch", "gap.opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit.score")
card.dir="/atium/Data/Reference/CARD/"
card.tab=read_csv(file.path(card.dir, "aro.csv")) %>%
    mutate(Description=gsub('\n', '', Description))


setwd(workdir)

##Run Prokka
system(paste0("~/Code/prokka/bin/prokka --outdir ", workdir, "/", prefix, "_genus --genus Klebsiella --usegenus --prefix ",
              prefix, "_genus ", assfile))

file.copy(file.path(workdir, paste0(prefix, "_genus"), paste0(prefix, "_genus.gff")),
          file.path(outdir, "annot", paste0(prefix, "_genus.gff")))

##Run RGI (CARD)
system(paste0("python ~/Code/rgi_card/rgi.py -t contig -i ", assfile, " -n 8 -o ", prefix, ".rgicard"))
system(paste0("python ~/Code/rgi_card/convertJsonToTSV.py -i ", prefix, ".rgicard.json -o ", prefix, ".rgicard"))

rgi=read_tsv(file.path(workdir, paste0(prefix, ".rgicard.txt"))) %>%
    mutate(clabel=gsub('_[0-9]*$', '', ORF_ID)) %>%
    mutate(attrib=paste0("Name=",ARO_name, ";gene=", ARO_name, ";ARO_category=", AR0_category, ";cutoff=", CUT_OFF))
    
rgi.gff=rgi %>%
    mutate(source="rgicard", feature="CDS", frame="0") %>%
    select(seqname=clabel, source=source, feature=feature, start=START,
           end=STOP, score=Best_Identities, strand=ORIENTATION,
           frame=frame, attribute=attrib) 
    
write_tsv(rgi.gff, file.path(outdir, "annot", paste0(prefix, ".rgicard.gff")), col_names=F)
    
    
##Run CARD Blast
system(paste0("blastn -db /atium/Data/Reference/CARD/CARDblast -query ", assfile,
              " -evalue 1e-50 -outfmt 6 -best_hit_overhang .1 -best_hit_score_edge 1e-5 >", workdir, "/", prefix, ".cardblast.tsv"))
    
card.res=read_tsv(file.path(workdir, paste0(prefix, ".cardblast.tsv")), col_names=blast.cols) %>%
    mutate(aro=sapply(strsplit(subject.id, split="\\|"), function(x) {x[5]})) %>%
    mutate(aro.idx=pmatch(aro, card.tab$Accession, duplicates.ok=TRUE)) %>%
    mutate(aro.name=card.tab$Name[aro.idx]) %>%
    mutate(aro.desc=card.tab$Description[aro.idx])   
    
card.gff=card.res %>%
    mutate(source="cardblast", feature="CDS", frame="0") %>%
    mutate(strand=ifelse(s.end>s.start, yes="+", no="-")) %>%
    mutate(attrib=paste0("Name=", aro.name, ";gene=", aro.name, ";per.identity=", per.identity, ";bit.score=", bit.score,";evalue=",evalue, ";ARO=",aro, ";Description=",aro.desc)) %>%
    select(seqname=query.id, source=source, feature=feature, start=q.start, end=q.end, score=bit.score, strand=strand, frame=frame, attribute=attrib)

write_tsv(card.gff, file.path(outdir, "annot", paste0(prefix, ".cardblast.gff")), col_names=F)
    
    
    
##Look for specific gene mutations
mutgenedb="/atium/Data/Nanopore/Analysis/161223_kpneumo/genedb/mutdbase"
system(paste0("blastn -db ", mutgenedb, " -query ", assfile, " -outfmt 6 >", workdir, prefix, ".mutblast.tsv"))

mut.res=read_tsv(file.path(workdir, paste0(prefix, ".mutblast.tsv")), col_names=blast.cols)

mut.gff=mut.res %>%
    mutate(source="mutblast", feature="CDS", frame="0") %>%
    mutate(strand=ifelse(s.end>s.start, yes="+", no="-")) %>%
    mutate(attrib=paste0("Name=", subject.id, ";gene=", subject.id, ";per.identity=", per.identity, ";bit.score=", bit.score,";evalue=",evalue)) %>%
    select(seqname=query.id, source=source, feature=feature, start=q.start, end=q.end, score=bit.score, strand=strand, frame=frame, attribute=attrib)

write_tsv(mut.gff, file.path(outdir, "annot", paste0(prefix, ".mutblast.gff")), col_names=F)


##look for specific is recog sequences
isedb="/atium/Data/Nanopore/Analysis/161223_kpneumo/genedb/isecp1blast"
system(paste0("blastn -task blastn-short -db ", isedb, " -query ", assfile, " -outfmt 6 >", workdir, prefix, ".iseblast.tsv"))

ise.res=read_tsv(file.path(workdir, paste0(prefix, ".iseblast.tsv")), col_names=blast.cols)

ise.gff=ise.res %>%
    mutate(source="iseblast", feature="Mobile_element", frame="0") %>%
    mutate(strand=ifelse(s.end>s.start, yes="+", no="-")) %>%
    mutate(attrib=paste0("Name=", subject.id, ";gene=", subject.id, ";per.identity=", per.identity, ";bit.score=", bit.score,";evalue=",evalue)) %>%
    select(seqname=query.id, source=source, feature=feature, start=q.start, end=q.end, score=bit.score, strand=strand, frame=frame, attribute=attrib)

write_tsv(ise.gff, file.path(outdir, "annot", paste0(prefix, ".iseblast.gff")), col_names=F)


##MLST typing
system(paste0("~/Code/mlst/bin/mlst ", assfile, " >", workdir, prefix, ".mlst.tsv"))
    




