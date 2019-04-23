
##########################
message("#################################################################################################################")
message("#StrainR  v0.2")
message("#J Bisanz 18 Apr 2019")
message("#filters reads, then maps to preprocessed index and normalizes")
message("#################################################################################################################")
message(" ")

##########################
#Get arguments
suppressMessages(library(optparse))
option_list = list(
  make_option(c("-1", "--forward"), type="character", help="forward read in fastq.gz format", metavar="character"),
  make_option(c("-2", "--reverse"), type="character", help="reverse read in fastq.gz format", metavar="character"),
  make_option(c("-r", "--reference"), type="character", default="DB_StrainR", help="reference directory from PreProcessR.R", metavar="character"),
  make_option(c("-t", "--threads"), type="numeric", default="4", help="number of cpus to use /threads to launch [default= %default]", metavar="numeric"),
  make_option(c("-q", "--truncQ"), type="numeric", default="2", help="truncate reads after a quality score of this or lower [default= %default]", metavar="numeric"),
  make_option(c("-e", "--maxEE"), type="numeric", default="2", help="discard reads with more than this many expected errors [default= %default]", metavar="numeric"),
  make_option(c("-b", "--bbmap"), type="character", default="bbmap.sh", help="location of bbmap.sh if not in path [default= %default]", metavar="character"),
  make_option(c("-c", "--samtools"), type="character", default="samtools", help="location of samtools if not in path [default= %default]", metavar="character"),
  make_option(c("-m", "--mem"), type="numeric", default="12", help="amount of RAM to use for mapping [default= %default]", metavar="numeric"),
  make_option(c("-o", "--outdir"), type="character", default="Normalized_StrainR", help="output directory for normalized data", metavar="character"),
  make_option(c("-p", "--outprefix"), type="character", default="sample", help="prefix for output files/ a sample identifier", metavar="character"),
  make_option(c("-s", "--scratch"), type="character", default="tmp",help="directory for temporary files", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$forward)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

message(date(), "---> Loading libraries")
suppressMessages(library(dada2))
suppressMessages(library(Biostrings))
suppressMessages(library(tidyverse))
suppressMessages(library(ShortRead))

##########################
dir.create(opt$scratch, showWarnings = FALSE)
dir.create(opt$outdir, showWarnings = FALSE)
dir.create(paste0(opt$outdir,"/figures"), showWarnings = FALSE)
opt$filtforward<-gsub("..+\\/", paste0(opt$scratch,"/"), opt$forward)
opt$filtreverse<-gsub("..+\\/", paste0(opt$scratch,"/"), opt$reverse)
##########################


if(file.exists(paste0(opt$outdir, "/", opt$outprefix,".bam"))){ 
  message(date(), "---> Skipping read mapping for ", opt$outprefix," as bam file already found")
}else{

  message(date(), "---> Trimming reads for ", opt$outprefix)
  
  #trim reads
  filterAndTrim(
    fwd=opt$forward, filt=opt$filtforward, 
    rev=opt$reverse, filt.rev=opt$filtreverse,
    compress = TRUE,
    truncQ = opt$truncQ,
    truncLen = 0,
    trimLeft = 10,
    maxLen = Inf,
    minLen = 50,
    maxN = 0,
    minQ = 0,
    maxEE = opt$maxEE,
    rm.phix = TRUE,
    multithread = opt$threads,
    n=1e6,
    verbose = TRUE)
  
  message(date(), "--->Plotting quality profile for ", opt$outprefix)
  
  pdf(paste0(opt$outdir,"/figures/",opt$outprefix,"_readquality.pdf"), height=6, width=8, useDingbats=F)
    plotQualityProfile(c(opt$forward, opt$reverse))
  dev.off()
  
  ##########################
  #map reads
  message(date(), "--->Mapping reads for ", opt$outprefix)
  system(paste0(
    opt$bbmap,
    " in=", opt$filtforward,
    " in2=", opt$filtreverse,
    " ref=", opt$reference, "/BBindex/BBIndex.fasta",
    " out=", opt$outdir, "/", opt$outprefix,".sam",
    " rpkm=", opt$outdir, "/",opt$outprefix,".rpkm",
    " perfectmode=t ",
    " threads=",opt$threads,
    " local=f",
    " ambiguous=toss",
    " pairedonly=t",
    " -Xmx",opt$mem,"g",
    " nodisk=t"
  ))
  
  system(paste0(
    opt$samtools,
    " view",
    " --threads ", opt$threads,
    " -bS ", opt$outdir, "/", opt$outprefix,".sam",
    " | ",
    opt$samtools,
    " sort",
    " --threads ", opt$threads,
    " -o ", opt$outdir, "/", opt$outprefix,".bam"
  ))
  
  unlink(paste0(opt$outdir, "/", opt$outprefix,".sam"))
  unlink(opt$filtforward)
  unlink(opt$filtreverse)
}
##########################
#Get mappings and merge with fragment counts
message(date(), "--->Normalizing Data ", opt$outprefix)

mapped<-system(paste0("grep '#Mapped' ",opt$outdir, "/",opt$outprefix,".rpkm"), intern = TRUE) %>% gsub("#Mapped\t","", .) %>% as.numeric()

Norm<-
  read_tsv(paste0(opt$outdir, "/",opt$outprefix,".rpkm"), skip=4) %>% 
  dplyr::rename(FragmentID=`#Name`) %>%
  left_join(
    read_tsv(paste0(opt$reference,"/KmerContent.report"), col_types="ccccid")
  ) %>%
  mutate(Total_Mapped_Reads_In_Sample=mapped) %>%
  select(StrainID, ContigID, Start_Stop, Unique_Kmers=Nunique, Length_Contig=Length, Bases, Coverage, Mapped_Reads=Reads, Mapped_Frags=Frags, Total_Mapped_Reads_In_Sample) %>%
  mutate(FUKM=Mapped_Frags/(Unique_Kmers/1e3)/(Total_Mapped_Reads_In_Sample/1e6))

write_tsv(Norm, paste0(opt$outdir, "/", opt$outprefix,".abundances"))

message(date(), "---> Plotting normalized data for ", opt$outprefix)
pplot<-
  Norm %>%
  filter(Unique_Kmers>10) %>%
  ggplot(aes(x=StrainID, y=log10(FUKM), fill=StrainID)) +
  geom_boxplot(outlier.alpha=0) +
  geom_jitter(width=0.2, shape=21) +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ggtitle(opt$outprefix)
ggsave(paste0(opt$outdir, "/figures/", opt$outprefix,".pdf"),pplot, device="pdf", height=8, width=11, useDingbats=F )

unlink(opt$scratch, recursive = TRUE)
message(date(), "--->--->---> StrainR complete for ", opt$outprefix)
