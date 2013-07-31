#################################################
##  Script to get coordinates of primers and      
##        Generate manifest file for      	
##	 Targeted resequencing on the MiSeq 				
##  Aparicio Lab WSOP 2013-001 developed by			
##	 Dr Damian Yap , Research Associate				
##	 dyap@bccrc.ca  Version 2.0 (Jul 2013) 			
##################################################


# These commands must be specifed in order for this script to work
# source("http://www.bioconductor.org/biocLite.R"); source("http://www.bioconductor.org/biocLite.R"); biocLite("BSgenome"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); library('BSgenome.Hsapiens.UCSC.hg19')
# install.packages("XLConnect")
# install.packages("XML", repos = "http://www.omegahat.org/R")


library(Biostrings)
library("IRanges")
library("GenomicRanges")
library(Rsamtools)
library('BSgenome.Hsapiens.UCSC.hg19')
library(XLConnect)
library(XML)

# Load the latest available dnSNP for the version of bioconductor if not installed 
len <-  length(available.SNPs())
dbSNP <- available.SNPs()[len]
print("Latest SNP database")
dbSNP

SNP <-   installed.SNPs()
print("Installed SNP database")
dbSNP

# Inject the SNPs into the hg19 ref
SNP_Hsapiens <- injectSNPs(Hsapiens, dbSNP)

#################################################
# Directory structure and file names

basedir="/share/lustre/backup/dyap/Projects/MiSeq_Data"

primerdir=paste(basedir,"Primer_Order_Files", sep="/")

manifestdir=paste(basedir,"AmpliconManifests", sep="/")

######################################################

# Input #1 is the primer excel sheet
setwd(primerdir)

file_names = list.files(pattern="*Primers*")
files = paste(primerdir,file_names, sep="/")

#################
#for (fname in files)
#{ 
# For testing
fname<-files[1]
	# Input #2 which is the AmpliconManfest file
	# set input directory
	indir=manifestdir
	# This strips the path and characters after sampel name (check when changed)
	sample=substring(fname,66,70)
	manfile=paste(sample,"AmpliconManifest",sep=".")
	infile2=paste(indir, manfile, sep="/")

# workbook <- system.file(file, package = "XLConnect")
# Load workbook
wb <- loadWorkbook(fname, create = FALSE)
# Query available worksheets
sheets <- getSheets(wb)

# This command reads all the sheets in a workbook into a data_list
primers <- readWorksheet(wb, sheets)

############### Matching the forward and reverse primers ########################
# calculating the number of plates
plates <- length(names(primers))
pairs <- plates/2
if (!isTRUE(all(pairs == floor(pairs)))) stop("Odd number of primer plates cannot be matched")

# If the script continues that means we have even number of plates
# Prepare arrays (data.frames) for primer sets (F and R)

for (pr in seq(plates)) {

	# Fwd and Rev_RC are the primer sequences that were ordered

   d.frame <- data.frame(Plate = rep("", nrow(primers[[plates]])),
		     Well = rep("", nrow(primers[[plates]])),
                     Name = rep("", nrow(primers[[plates]])),
                     Fwd = rep("", nrow(primers[[plates]])),
                     Rev = rep("", nrow(primers[[plates]])),
                     Rev_RC = rep("", nrow(primers[[plates]])),
                     stringsAsFactors = FALSE)

# This only assigns the data frame to the primer pairs
#pair <- pr/2
#print(pair)
#if (isTRUE(all(pair == floor(pair)))) assign(paste("Plate", pair, sep=""), d.frame)
}      

# Sections work till here
# This takes the forward primers and puts them into a dataframe.
# This function takes the last n characters of a string
	substrRight <- function(x, n)	{
  					substr(x, nchar(x)-n+1, nchar(x))
					}


for (rw in seq(names(primers)))  {

			platename <- names(primers)[rw]
			plate <- substr(platename,1,nchar(platename)-1)
			plateid <- substr(platename,nchar(platename)-1,nchar(platename))
			test <- substrRight(platename,1)
			if (test == "R")
				next
				else { 
 			for (rx in seq(nrow(primers[[plates]]))) {

						d.frame$Plate[rx] <- plate
						d.frame$Well[rx] <- primers[[rw]]$Well[rx]
						d.frame$Name[rx] <- substr(primers[[rw]]$Name[rx],1,nchar(primers[[rw]]$Name[rx])-2)


						left <-  primers[[rw]]$Sequence[rx]

				# Get the matching REV plate by primer name (cuts F and R suffixies)

						 matchF <- substr(primers[[rw]]$Name[rx],1,nchar(primers[[rw]]$Name[rx])-2) 
						 matchR <- substr(primers[[rw+1]]$Name[rx],1,nchar(primers[[rw+1]]$Name[rx])-2)
						if ( matchF == matchR ) right <- primers[[rw+1]]$Sequence[rx]

						# Removal of adaptors (Fluidigm)
						leftadapt="ACACTGACGACATGGTTCTACA"
						# (5'->3' of reverse adaptor)
						rightadapt="TACGGTAGCAGAGACTTGGTCT"

				# removal of the adaptor sequences from ordered primers
					lpriseq <- gsub(leftadapt, "", left)
					rpriseq <- gsub(rightadapt, "", right)
		
				# Reverse Complement the Right primer without adaptor
					x <- DNAString(rpriseq)
					rpriseqr <- as.character(reverseComplement(x))
		
					d.frame$Fwd[rx] <- lpriseq
					d.frame$Rev[rx] <- rpriseqr
					d.frame$Rev_RC[rx] <- rpriseq
		
				# for first value 
				if ( rw == 1 ) {
						first <- d.frame }
 
 				# Reverse plates are skipped
				# combining successive data.frames
				if ( rw == 3 ) { 
						sum1 <- rbind(first,d.frame) }

				# combining successive data.frames
				if ( rw > 3 ) { 
						sum1 <- rbind(sum1,d.frame) 
							} else { print("skip") }

                     					}
					}
				}


################### READ IN PRIMER ORDER FILE AND REMOVE BITS #####################
# Need to relabel primers with amplicon IDs 
# Read in the primer order file and remove primer= etc and secondary primers

# This file is all in one column

orderdir="/home/dyap/Projects/Tumour_Evol/positions/SNV"
orderfile=paste(paste(orderdir,sample,sep="/"),"p3_order.txt",sep="_")

file <- read.table(orderfile, sep="\n", skip=2)

pat="PRIMER_LEFT_1_SEQUENCE="
file1 <- subset(file, !grepl(pat, file$V1))

pat="PRIMER_RIGHT_1_SEQUENCE="
file2 <- subset(file1, !grepl(pat, file1$V1))

# This removes the primer_0_ headings
clean <- function(ttt){
	gsub("PRIMER_LEFT_0_SEQUENCE=", "", ttt)	
			}
file2[] <- sapply(file2, clean)

clean <- function(ttt){
	gsub("PRIMER_RIGHT_0_SEQUENCE=", "", ttt)
			}
file2[] <- sapply(file2, clean)

################# EXTRACT WITH CONTEXT ########################
extract.with.context <- function(x, rows, after = 0, before = 0) {

  match.idx  <- which(rownames(x) %in% rows)
  span       <- seq(from = -before, to = after)
  extend.idx <- c(outer(match.idx, span, `+`))
  extend.idx <- Filter(function(i) i > 0 & i <= nrow(x), extend.idx)
  extend.idx <- sort(unique(extend.idx))

  return(x[extend.idx, , drop = FALSE])
}

####### Inject the chr-pos amplicon ID into sum1 ###########

for (ri in seq(nrow(sum1)))
	{
	saname <-  substr(sum1[ri,3],1,5)
	leftpri <- sum1[ri,4]
	# The right primer is reverse complemented
	rightprir <- sum1[ri,5]
	# Reverse Complement the Right primer 
					x <- DNAString(rightprir)
					rightpri <- as.character(reverseComplement(x))

	#This grep command with $V1 specified returns the COLUMN match
	leftmatch <- grep(leftpri,file2$V1)
	rightmatch <- grep(rightpri,file2$V1)

	# The ID chr-pos is 1 before leftmatch and 2 before right match in one column file2

	if (file2$V1[leftmatch-1] == file2$V1[rightmatch-2]) 
				{
				if (length(file2$V1[rightmatch-2]) != 1)
					{
					label <- capture.output(cat(file2$V1[leftmatch-1],sep="@"))
					} else {
						label <- file2$V1[leftmatch-1]
						}
				sum1[ri,3] <- label
				}
	}


######################### Ensure amplicons are unique ####################
# We do this by assigning uniquely assigning the amplicon name to each of the duplicate amplicons

replicates <- sum1$Name[duplicated(sum1$Name)]
for (rj in replicates)
		{
		# This returns the rows of the match
		dup <- grep(rj, sum1$Name)
		splitnames <- strsplit(sum1[dup[1],3], split="@")
		
		for (count in seq(length(dup)))
				{
				sum1[dup[count],3] <- splitnames[[1]][count]		
				}
		}


################## Get hg19 positions from UCSC inSilico PCR web #################################

PCR <- data.frame(myID = rep("", nrow(sum1)),
                     UCSCID = rep("", nrow(sum1)),
                     UCSCchr = rep("", nrow(sum1)),
                     UCSCstart = rep(0, nrow(sum1)),
                     UCSCEnd = rep(0, nrow(sum1)),
                     UCSCAmplen = rep("", nrow(sum1)),
                     UCSCLpri = rep("None", nrow(sum1)),
                     UCSCRpri = rep("NA", nrow(sum1)),
                     UCSCAmp = rep("NA", nrow(sum1)),
                     UCSCOutput = rep("NA", nrow(sum1)),                      
                     stringsAsFactors = FALSE)

for (no in seq(nrow(sum1)))
	{
	ampid <- sum1[no,3]
	fwd <- sum1[no,4]
	rev <- sum1[no,5]
	add1="http://genome.ucsc.edu/cgi-bin/hgPcr?hgsid=342829929&org=Human&db=hg19&wp_target=genome&wp_f="
	add2="&Submit=submit&wp_size=400&wp_perfect=15&wp_good=15&wp_flipReverse=on&boolshad.wp_flipReverse=0"
	sep="&wp_r="
	url=paste(paste(paste(add1,fwd,sep=""),rev,sep=sep),add2,sep="")

	doc2<-readHTMLTable(url)
	parse<-as.character(doc2[[1]][2,1])
	print(ampid)
	
	PCR$UCSCOutput[no] <-parse
	}	

write.csv(file=PCR, header=TRUE)

# started at 7 check time out file to see what time it ended

# Output file
outfile=paste(paste(workdir,sample,sep="/"), "SNP_Checked_all.csv", sep="/")

manifest=as.data.frame(read.table(infile2, header=TRUE, skip=5, sep="\t", stringsAsFactors = FALSE))

			
outdf <- data.frame(ID = rep("", nrow(manifest)),

                     AmpSNP = rep("", nrow(manifest)),
                     LpriSNP = rep("", nrow(manifest)),
                     RpriSNP = rep("", nrow(manifest)),
                     Lpriseq = rep("None", nrow(manifest)),
                     Rpriseq = rep("None", nrow(manifest)),
                     Ampliseq = rep("None", nrow(manifest)),
                     SNPLpriseq = rep("NA", nrow(manifest)),
                     SNPRpriseq = rep("NA", nrow(manifest)),
                     SNPAmpliseq = rep("NA", nrow(manifest)),                      
                     stringsAsFactors = FALSE)

for (ri in seq(nrow(manifest))) {
  
	        	id <- manifest$Name[ri]
    			chr <- manifest$Chromosome[ri]
    			
    			# We have to correct since there was a shift in the manifest
  			start <-  as.numeric(manifest$Amplicon.Start[ri])+2
			end <-  as.numeric(manifest$Amplicon.End[ri])-2
			leftlen <-  as.numeric(manifest$Upstream.Probe.Length[ri])-2
			rgtlen <-  as.numeric(manifest$Downstream.Probe.Length[ri])-2
			
			# find id from manifest in the col 1 (amplicon ID), col 2 (left) and 3 (right)
			leftprimer <-  primerlist[primerlist[,1] %in% c(id),2]
			rightprimer <-  primerlist[primerlist[,1] %in% c(id),3]

# Get the amplicons 1. Reference and 2. SNPmasked
 						ampliseq <- as.character(getSeq(Hsapiens,chr,start,end))
 						snpampliseq <- as.character(getSeq(SNP_Hsapiens,chr,start,end))
 
# Get left and right primers 1. Designed Seq 2. SNPmasked 
						leftend <- start + leftlen
						lpriseq <- leftprimer
						snplpriseq <- as.character(getSeq(SNP_Hsapiens,chr,start,leftend))

						rightstart <- end - rgtlen
						rpriseq <- rightprimer
						snprpriseq <- as.character(getSeq(SNP_Hsapiens,chr,rightstart,end))

						# Testing to see if the sequence are identical, if they are they do not contain SNPs

						if (ampliseq == snpampliseq) ampsnp <- "ok" else { ampsnp <- "SNP" }
						if (lpriseq == snplpriseq) lprisnp <- "ok" else { lprisnp <- "SNP" }
						if (rpriseq == snprpriseq) rprisnp <- "ok" else { rprisnp <- "SNP" }


# writing the output to the dataframe 
			outdf$ID[ri] <- id
                      	outdf$AmpSNP[ri] <- ampsnp
                     	outdf$LpriSNP[ri] <- lprisnp
                     	outdf$RpriSNP[ri] <- rprisnp
                     	outdf$Lpriseq[ri] <- lpriseq
                     	outdf$Rpriseq[ri] <- rpriseq
                     	outdf$Ampliseq[ri] <- ampliseq
                     	outdf$SNPLpriseq[ri] <- snplpriseq
                     	outdf$SNPRpriseq[ri] <- snprpriseq
                     	outdf$SNPAmpliseq[ri] <- snpampliseq

						
  }
  
# output file
write.csv(outdf, file = outfile)

}
