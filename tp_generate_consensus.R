# Adapted for Nextflow Aug 2020

# RSV : This script imports bam files and makes a consensus sequence
# Pavitra Roychoudhury
# Adapted from hsv_generate_consensus.R on 6-Mar-19

# Built to be called from hhv6_wgs_pipeline.sh with input arguments specifying input filename
# Requires wgs_functions.R which contains several utility scripts plus multiple R packages listed below

rm(list=ls()); 
sessionInfo();
library(Rsamtools);
library(GenomicAlignments);
library(ShortRead);
library(Biostrings);
library(RCurl);

#Get args from command line 
args<-(commandArgs(TRUE));
if(length(args)==0){
	print("No arguments supplied.")
}else{
	for(i in 1:length(args)){
		eval(parse(text=args[[i]]))
		print(args[[i]])
	}
}

#For testing (these args should come from command line)
# sampname='2016-01040_S451_L001'
# ref='NC_016842' 

#Files, directories, target site
merged_bam_folder <- './'
mapped_reads_folder <- './'
con_seqs_dir <- './'

#merged_bam_folder<-'./remapped_reads/'; 
#mapped_reads_folder<-'./mapped_reads/';
#con_seqs_dir<-'./consensus_seqs_all';

clean_consensus_tp<-function(sampname,merged_bam_folder,mapped_reads_folder,ref){
  mapping_stats<-data.frame(ref=ref,
                            bamfname_merged=grep(sampname,list.files(merged_bam_folder,'*remapped.sorted.bam$',full.names=T),value=T),
                            bamfname_mapped=grep(sampname,list.files(mapped_reads_folder,'*firstmap.sorted.bam$',full.names=T),value=T),
                            mapped_reads_ref=0,mapped_reads_assemblyref=0,perc_Ns=0,num_Ns=0,width=0,
                            stringsAsFactors=F);
  
  #Import mapped reads + assembly and generate consensus
  con_seqs<-lapply(mapping_stats$bamfname_merged,generate_consensus);
  dummyvar<-lapply(con_seqs,function(x)
    writeXStringSet(x,file=paste('./',names(x),'_consensus.fasta',sep=''),format='fasta'));
  rm(dummyvar)
  
  #Compute #mapped reads and %Ns
  mapping_stats$mapped_reads_ref<-unlist(lapply(mapping_stats$bamfname_mapped,n_mapped_reads));
  mapping_stats$mapped_reads_assemblyref<-unlist(lapply(mapping_stats$bamfname_merged,n_mapped_reads));
  mapping_stats$num_Ns<-unlist(lapply(con_seqs,function(x)sum(letterFrequency(x,c('N','+')))));
  mapping_stats$width<-unlist(lapply(con_seqs,width));
  mapping_stats$perc_Ns<-100*mapping_stats$num_Ns/mapping_stats$width;
  write.csv(mapping_stats,file=paste('./',sampname,'_mappingstats.csv',sep=''),row.names=F);
  
  return(TRUE)
}

#Make consensus sequence--returns TRUE if this worked
conseq<-clean_consensus_tp(sampname,'./','./',ref);

#Prepare seqs for annotation -- will make separate folders for A and B
if(conseq==TRUE){
  #Remove all Ns at the beginning and end of the seq, write to folder
  fname<-grep(sampname,list.files(con_seqs_dir,full.names=T),value=T);
  con_seq<-readDNAStringSet(fname);
  con_seq_trimmed<-DNAStringSet(gsub("N*N$",'',gsub("^N*",'',as.character(con_seq))));
  names(con_seq_trimmed)<-substring(names(con_seq),1,20); #prokka needs contig name to be <=20 chars long
  writeXStringSet(con_seq_trimmed,file=paste(sampname,'_preprokka_consensus.fasta',sep=''),format='fasta');
	
}else{
	print('Failed to generate consensus sequences.')
}

