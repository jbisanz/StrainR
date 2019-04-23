
##########################
message("#################################################################################################################")
message("#PreProcessR v0.2")
message("#J Bisanz 28 March 2019")
message("#fragments reference genomes and calculates unique kmers on a per fragment basis, then assembles")
message("#################################################################################################################")
message(" ")

##########################
#Get arguments
suppressMessages(library(optparse))
option_list = list(
  make_option(c("-i", "--indir"), type="character", help="folder containing reference genomes in fasta format. Names must not contain any other field separator than _ and only fasta files in the directory", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="StrainRDB", help="folder to build output database in [default= %default]", metavar="character"),
  make_option(c("-t", "--threads"), type="numeric", default="2", help="number of cpus to use /threads to launch [default= %default]", metavar="numeric"),
  make_option(c("-n", "--jobs"), type="numeric", default="4", help="number of jobs that can run in parallel, total processor usage is jobs x threads [default= %default]", metavar="numeric"),
  make_option(c("-k", "--kmersize"), type="numeric", default="81", help="kmer size", metavar="numeric"),
  make_option(c("-f", "--fragmentsize"), type="numeric", default="50000", help="size to fragment reference genomes to [default= %default]", metavar="numeric"),
  make_option(c("-e", "--excludesize"), type="numeric", default="10000", help="exclude fragments below this size [default= %default]", metavar="numeric"),
  make_option(c("-b", "--bbmap"), type="character", default="bbmap.sh", help="location of bbmap.sh if not in path [default= %default]", metavar="character"),
  make_option(c("-j", "--jellyfish"), type="character", default="jellyfish", help="location of jellyfish if not in path [default= %default]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$indir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

##########################

dir.create(opt$outdir)
#######################################################
# Load required libraries
message(date(), "---> Loading libraries")
suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))
suppressMessages(library(doParallel))
suppressMessages(library(foreach))
suppressMessages(library(data.table))
suppressMessages(library(Matrix))
suppressMessages(library(Matrix.utils))
suppressMessages(library(vegan))
suppressMessages(library(openssl))



# Load in all reference contigs, fragment and rename
# Assuming that the strain name does not have periods
# All fragments are labelled as Strain;ContigID;Start_Stop;Length
if(file.exists(paste0(opt$outdir,"/Fragments/"))){
  message(date(), "---> References already fragmented ...skipping")
} else {
  message(date(), "---> References being fragmented to ", opt$fragmentsize, "bp fragments.")
  dir.create(paste0(opt$outdir, "/Fragments"))
  
  Refs<-list.files(opt$indir)
  #FragmentLookup<-tibble()
  pclust<-makeForkCluster(nnodes = opt$jobs*opt$threads)
  registerDoParallel(pclust)
  lookup<-foreach(Ref=Refs) %dopar% {
    cltr<-tibble()
    #message("Fragmenting:", Ref)
    Strain<-gsub("\\...+","", Ref)
    Fragments<-readDNAStringSet(paste0(opt$indir,"/", Ref))
    Fragments<-lapply(Fragments, function(x) Views(x, successiveIRanges(rep(opt$fragmentsize, ceiling(length(x)/opt$fragmentsize)))))
    
    Fragments<-lapply(names(Fragments), function(x){
      suppressWarnings(y<-DNAStringSet(Fragments[[x]]))
      names(y)<-paste0(Strain,";",x,";",y@ranges@start,"_", y@ranges@start + y@ranges@width, ";",y@ranges@width)
      return(y)
    })
    
    Fragments<-do.call(c, Fragments)
    #message("--->Fragmented into :", length(Fragments), " fragments")
    Fragments<-Fragments[sapply(Fragments, length)>opt$excludesize] #remove fragments smaller than excludesize
    
    #message("--->After size filtering there are :", length(Fragments), " fragments")
    
    #writeXStringSet(Fragments[[Strain]], paste0(opt$outdir,"/Fragments/",Strain,".frags"))
    for(frag in names(Fragments)){
      tmp<-DNAStringSet(Fragments[[frag]])
      tmp<-c(tmp, reverseComplement(tmp))
      names(tmp)<-c(frag, paste0(frag," RC"))
      fname<-md5(frag)
      suppressWarnings(cltr<-bind_rows(cltr, tibble(Fragment=frag, md5=fname)))
      writeXStringSet(tmp, paste0(opt$outdir,"/Fragments/",fname,".frags"))
      rm(tmp)
    }
    rm(Fragments)
    return(cltr)
  }
  #rm(Ref, Strain)
  do.call(bind_rows, lookup) %>% write_tsv(paste0(opt$outdir,"/Frags.md5"))
  stopCluster(pclust)
  #write_tsv(FragmentLookup, paste0(opt$outdir, "/Fragments/FragmentLookup.txt"))
}

gc(verbose=TRUE)
message(date(), " Fragmentation complete")

#######################################################
# Now generate kmer profiles using Jellyfish 2

if(file.exists(paste0(opt$outdir,"/Kmers"))){
  message(date(), "---> Kmers already generated ...skipping")
} else {
  message(date(), "---> Canonical ", opt$kmersize, "mers being generated")
  dir.create(paste0(opt$outdir, "/Kmers"))
  
  Refs<-list.files(paste0(opt$outdir, "/Fragments")) %>% grep("\\.frags$", .,  value=T)
  Refs<-split(Refs, sample(1:opt$jobs, length(Refs), replace=T))
  pclust<-makeForkCluster(nnodes = opt$jobs)
  registerDoParallel(pclust)
  foreach(chunk=Refs) %dopar% {
    for(Ref in chunk){
      Strain<-gsub("\\.frags","", Ref)
      system(paste0(
        "jellyfish count -C",
        " --size=20M",
        " --mer-len=", opt$kmersize,
        " --threads=", opt$threads,
        " --output ",  opt$outdir,"/Kmers/",Strain, ".jf",
        " ",opt$outdir,"/Fragments/",Ref
      ))
      
      system(paste0(
        opt$jellyfish,
        " dump",
        " --output ",  opt$outdir,"/Kmers/",Strain, ".kmers",
        " ",opt$outdir,"/Kmers/",Strain, ".jf"
      ))
      
      system(paste0("rm ", opt$outdir,"/Kmers/",Strain, ".jf"))
    }
  }  
  stopCluster(pclust)
}

gc(verbose=TRUE)
message(date(), " Kmer generation complete")

#######################################################
# Now build a list of all observed kmers

if(file.exists(paste0(opt$outdir,"/Kmers.txt"))){
  message(date(), "---> Kmers list already built... skipping")
} else {
  message(date(), " Building Kmer list")
  system(paste0("sed -n '1~2!p' ", opt$outdir,"/Kmers/*kmers > ",opt$outdir,"/Kmers2rep.txt")) # discard fasta seqs
  system(paste0("sort --parallel=",opt$jobs*opt$threads," ",opt$outdir,"/Kmers2rep.txt > ",opt$outdir,"/Kmers2un.txt"))# Sort
  system(paste0("uniq ",opt$outdir,"/Kmers2un.txt > ",opt$outdir,"/Kmers.txt"))# Keep only unique Kmers
  system(paste0("rm ", opt$outdir,"/Kmers2rep.txt"))
  system(paste0("rm ", opt$outdir,"/Kmers2un.txt"))
}


#######################################################
# Now Generate Kmer Indicies

if(file.exists(paste0(opt$outdir,"/Kindex"))){
  message(date(), "---> Kmers indices already generated... skipping")
} else {
  dir.create(paste0(opt$outdir,"/Kindex"))
  message(date(), " Building Kmer indicies")
  Kmers<-fread(paste0(opt$outdir,"/Kmers.txt"), col.names = "Kmer", header=F)
  Frags<-list.files(paste0(opt$outdir,"/Kmers"))
  Frags<-split(Frags, sample(1:(opt$jobs*opt$threads), length(Frags), replace=T))
  pclust<-makeForkCluster(nnodes = (opt$jobs*opt$threads), outfile=paste0(opt$outdir,"/indexing_log.txt"))
  registerDoParallel(pclust)
  
  jsink<-foreach(Frag=Frags) %do% {
    for(Fragment in Frag){
      Ks<-readDNAStringSet(paste0(opt$outdir,"/Kmers/",Fragment)) %>% as.character(., use.names=FALSE)
      idx<-which(Kmers$Kmer %in% Ks)
      saveRDS(idx,paste0(opt$outdir, "/Kindex/",Fragment,".RDS"))
      rm(idx, Ks)
    }
    return("complete")
  }
    stopCluster(pclust)
  message(date(), " Kmer indices generated")
  gc(verbose=T)
}

#######################################################
# Now Generate Kmer sparse matrix

lookup<-read_tsv(paste0(opt$outdir,"/Frags.md5"))
if(file.exists(paste0(opt$outdir,"/kmat.RDS"))){
  message(date(), "---> Kmer sparse matrix already built... skipping")
  kmat<-readRDS(paste0(opt$outdir,"/kmat.RDS"))
}else {
  message(date(), " Building Kmer sparse matrix")
  idxs<-list.files(paste0(opt$outdir,"/Kindex/"))#[1:10]
  ks<-lapply(idxs, function(x) readRDS(paste0(opt$outdir,"/Kindex/", x)))
  
  j=as.integer()
  for(i in 1:length(ks)){
    j<-c(j, rep(i, length(ks[[i]])))
  }
  
  kmat<-sparseMatrix(i = unlist(ks), j=j ,x=1)
  rm(j)
  colnames(kmat)<-lookup[match(gsub("\\.kmers\\.RDS","", idxs),lookup$md5),]$Fragment
  saveRDS(kmat, paste0(opt$outdir,"/kmat.RDS"))
}


#######################################################
# Now process the matrix

data.frame(Nunique=colSums(kmat[rowSums(kmat)==1,])) %>% # Keep only the singleton kmers
  rownames_to_column("FragmentID") %>%
  as.tibble() %>%
  separate(FragmentID, c("StrainID","ContigID","Start_Stop","Length"),";", remove = FALSE) %>%
  write_tsv(paste0(opt$outdir,"/KmerContent.report"))

kmat<-t(kmat)
kmat<-aggregate(kmat, gsub(";..+","", rownames(kmat)), fun="sum")
data.frame(Nunique=rowSums(kmat[,colSums(kmat)==1])) %>% # Keep only the singleton kmers
  rownames_to_column("StrainID") %>%
  write_tsv(paste0(opt$outdir,"/KmerContent_perStrain.report"))


#######################################################
# Now Generate Kmer distance matrix

if(file.exists(paste0(opt$outdir,"/Kmerclusters.pdf"))){
  message(date(), "---> Distance matrix already built... skipping")
} else{
  message(date(), " Building distance matrix and plot")
  #kdist<-vegdist(as.matrix(kmat), method="jaccard", binary=TRUE)
  kdist<-dist(kmat, method="binary")
  saveRDS(kdist, paste0(opt$outdir,"/KmerDist.RDS"))
  
  kclust<-hclust(kdist, method="average")
  saveRDS(kclust, paste0(opt$outdir,"/KmerClust.RDS"))
  
  pdf(paste0(opt$outdir,"/Kmerclusters.pdf"))
    plot(kclust, xlab = "StrainID", ylab = "Binary Distance", main = "UPGMA clustering of Strain Kmer Profiles")
  dev.off()
}

#######################################################
# Finally, build the BBmap index
message(date(), " Building BBMAP index")
dir.create(paste0(opt$outdir,"/BBindex"))
system(paste0("cat ", opt$outdir, "/Fragments/*frags > ", opt$outdir,"/BBindex/BBIndex.fasta"))

idx<-readDNAStringSet(paste0(opt$outdir,"/BBindex/BBIndex.fasta"))
idx<-idx[!grepl(" RC$", names(idx))]
writeXStringSet(idx, paste0(opt$outdir,"/BBindex/BBIndex.fasta"))

system(paste0(
  opt$bbmap,
  " ref=",opt$outdir,"/BBindex/BBIndex.fasta",
  " path=",opt$outdir,"/BBindex/"
  ))

message(date(), " PreProcessR complete")

