
library("data.table")
library("GenomicRanges")
library("stringr")
library("seqinr")

########################################################### Functions ####################################################################
##########################################################################################################################################
##########################################################################################################################################

read.fasta.1 <- function(file,force2upper = TRUE,rna2dna = TRUE,seqonly = FALSE,subseq=FALSE,start_sub=0)
{
  
  lines <- readLines(file)
  
  ind <- which(substr(lines, 1L, 1L) == ">")
  
  nseq <- length(ind)
  if(nseq == 0){
    stop("no line starting with a > character found")
  }
  
  start <- ind + 1
  end <- ind - 1
  end <- c(end[-1], length(lines))
  
  sequences <- lapply(seq_len(nseq), function(i) paste(lines[start[i]:end[i]], collapse = ""))
  if(seqonly) return(sequences)
  
  nomseq <- lapply(seq_len(nseq), function(i){
    firstword <- strsplit(lines[ind[i]], " ")[[1]][1]
    substr(firstword, 2, nchar(firstword))
  })
  
  if(force2upper){
    sequences <- as.list(toupper(sequences))
  }
  
  
  if(rna2dna){
    sequences <- gsub("U","T",sequences)
  }
  
  names(sequences) <- nomseq
  
  if(subseq==TRUE){
    sequences=substr(sequences,rep(start_sub,length(sequences)),stri_length(sequences))
  }
  return(sequences)
}

exon_features<- function (exons){
  exons[strand == "+"] <- exons[strand == "+"][order(start,transcript_id)]
  exons[strand == "-"] <- exons[strand == "-"][order((-start),`transcript_id`)]
  
  exons$length = exons$end - exons$start +1
  exons <- exons[,region.start := start,by=c("transcript_id")]
  exons <- exons[,region.end := end,by=c("transcript_id")]
  exons <- exons[,region.length := length,by=c("transcript_id")]
  exons <- exons[,region.length := sum(length),by=c("transcript_id")]
  
  exons$exonRelEnd <-  ave(exons$length, exons$`transcript_id`, FUN=cumsum)
  exons$exonRelStart <- exons$exonRelEnd - exons$length + 1
  exons$exon.identifier <- seq(1:NROW(exons))
  exons$seqnames <- paste0(exons$chr,":",exons$transcript_id,"|",exons$strand)
  return(exons)
}

absolute_to_relative <- function(snvs, exons=exons.relative){
  
    
    cols.to.keep <- c("chr","start","end","strand","exon.identifier","exonRelStart","exonRelEnd","transcript_id")
    exons <- exons[,..cols.to.keep]
    
    absolute.table.gr <- makeGRangesFromDataFrame(as.data.frame(snvs),keep.extra.columns = T)
    exons.gr <- makeGRangesFromDataFrame(as.data.frame(exons),keep.extra.columns = T)
    
    overlaps <- findOverlaps(absolute.table.gr,exons.gr)
    absolute.overlaps.gr <-absolute.table.gr[queryHits(overlaps)]
    exons.overlap.gr <-exons.gr[subjectHits(overlaps)]
    
    intersection <- setnames(as.data.table(pintersect(absolute.table.gr[queryHits(overlaps)],exons.gr[subjectHits(overlaps)]))[,c("start","end")],
                             c("element.abs.start.Ex","element.abs.end.Ex"))

    
    element.overlap <- as.data.table(absolute.overlaps.gr)[,c("seqnames","start","end","strand","key")]
    setnames(element.overlap,c("start","end"),c("element.abs.start","element.abs.end"))
    exons.overlap <- as.data.table(exons.overlap.gr); 
    setnames(exons.overlap,c("start","end","width"),c("exon.abs.start","exon.abs.end","exon.width"))
    
    element.exons.merged <- cbind(element.overlap,exons.overlap[,c("exon.abs.start","exon.abs.end","exonRelStart","exonRelEnd","transcript_id")])
    element.exons.merged <- cbind(element.exons.merged,intersection)
    
    
    element.exons.merged$element.abs.coord.ex <- paste0(element.exons.merged$element.abs.start.Ex,":",element.exons.merged$element.abs.end.Ex)
    element.exons.merged$key <- paste0(element.exons.merged$key,"@", element.exons.merged$transcript_id)
    
    element.exons.merged$length <- element.exons.merged$element.abs.end - element.exons.merged$element.abs.start + 1
    
    element.exons.merged <- element.exons.merged[,elementRelStart := exonRelStart + (element.abs.start.Ex - exon.abs.start)
                                                 ][,elementRelEnd := elementRelStart + length - 1
                                                   ][strand == "-", elementRelEnd := exonRelStart + (exon.abs.end - element.abs.start.Ex)
                                                     ][strand == "-", elementRelStart := elementRelEnd - length + 1]
   
    element.relative <- element.exons.merged[,c("key","elementRelStart", "elementRelEnd","transcript_id"), with=F]

    element.relative$key <- as.numeric(gsub("@.*","",element.relative$key))

    element.relative <-as.data.table(merge(snvs,element.relative,by="key"))
    

    return(element.relative)
  
}

extend_snvs <- function(snvs, spliced_transcripts.fa, flank = 25, mode){
  
  snvs <- merge(snvs, spliced_transcripts.fa, by = "transcript_id")
  snvs$length <- 2*as.numeric(flank) + 1
  
  snvs <- snvs[,extended_start := elementRelStart - flank][extended_start < 1, extended_start:=1]
               
  snvs <- snvs[,extended_end := elementRelEnd + flank][extended_end > region.length, extended_end:=region.length]
  
  snvs$extended_length <- snvs$extended_end - snvs$extended_start + 1
  
  snvs <- snvs[,flank := flank][extended_length < length, flank := floor(extended_length/2) - 1
                                  ][extended_length < length,extended_start := elementRelStart - flank
                                    ][extended_length < length,extended_end := elementRelEnd + flank]

  
  snvs$initial_header <- paste0(snvs$chr,":", snvs$start,"-", snvs$end,"(", snvs$strand, ")")
  snvs$header <- paste0(snvs$initial_header,"|",snvs$transcript_id, "|", snvs$extended_start,"-", snvs$extended_end)
  
  snvs <- snvs[, extendedSeq := str_sub(sequence, start = extended_start, end = extended_end)]
  
  if(mode == "mut"){
    
    snvs <- snvs[, extendedSeq_upstream := str_sub(extendedSeq, start = 1, end = flank)]
    snvs <- snvs[, extendedSeq_downstream := str_sub(extendedSeq, start = flank + 1, end = extended_end)]
    
    snvs$extendedSeq <- paste0(snvs$extendedSeq_upstream, snvs$ALLELE, snvs$extendedSeq_downstream)
  }
  
  snvs <- snvs[,c("header","extendedSeq","initial_header", "key", "transcript_id"), with=F]
  
  return(snvs)
   
}

##########################################################################################################################################
##########################################################################################################################################


#Initialize files------

#Tab delimited exon file.
exon_file <- "/data/exons_example.txt"

#Transcript sequences in fasta format
transcript_file_fasta <- "/data/transcripts_example.fa"

#Tab delimited file with SNVs
snvs_file <- "/data/snvs_example.txt"

#Output file
file.out <- "output.fa"

flank_nts <- 25

## Add "ref" for reference allele in the retrieved sequence, "mut" for mutation
ref_allele_mode <- "mut" 

#Read files-----

exons <- setnames(fread(exon_file,header=F,sep="\t"),c("chr","start","end","transcript_id","score","strand","type"))
exons <- exon_features(exons)
snvs <- setnames(fread(snvs_file, header=F, sep="\t"),
                 c("chr","start","end","REF","ALLELE","strand"))

#Read spliced transcript sequences-----

spliced_transcripts.fa<-read.fasta.1(transcript_file_fasta)
file_fasta <- as.data.table(names(spliced_transcripts.fa))
spliced_transcripts.fa <- as.data.table(spliced_transcripts.fa)
spliced_transcripts.fa <- setnames(cbind(file_fasta, spliced_transcripts.fa),c("transcript_id","sequence"))
spliced_transcripts.fa[,region.length := nchar(sequence) ]

#Define relative coordinates-----

exons.relative <- exon_features(exons)
snvs$key <- seq(1:NROW(snvs))
snvs <- absolute_to_relative(snvs,exons.relative)

#Extend relative coordinates-----

snvs.extended <- extend_snvs(snvs, spliced_transcripts.fa, flank = flank_nts, mode = ref_allele_mode)


#Write fasta file-----

write.fasta(as.list(snvs.extended$extendedSeq), snvs.extended$header, file.out, open = "w", nbchar = 60, as.string = FALSE)





