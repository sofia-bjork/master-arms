#!/usr/bin/env Rscript

### COI ###

# define and create output directory

MOTU_dir <- "output/COI/MOTU"
if(!dir.exists(MOTU_dir)) dir.create(MOTU_dir)

# Read the uclust.txt table (from swarm output)

uclust <- read.table("output/COI/swarm/uclust.txt", sep = "\t", stringsAsFactors = F)

# Delete the rows with value "C" in first column (to remove one of the S / C duplicate ASV names)

uclust2 <- subset(uclust, uclust[,1]!="C")

# Rename the 10th column if first column = S to give the most abundant ASv in a cluster a MOTU ID

uclust2[,10] <-ifelse(uclust2[,1]=="S",uclust2[,9],uclust2[,10])

# Subset the columns we need

motu_asv_list<-uclust2[,9:10]

# Remove all characters after _ (incl. _) 

motu_asv_list[] <- lapply(motu_asv_list, function(y) gsub("_.*","", y))

# Rename columns

colnames(motu_asv_list)<-c("ASV", "MOTU")

# Read the ASV count table (output from blank correction)

asv_counts <- read.table("output/COI/blank_corr/asv_no_contaminants_COI.txt", sep = "\t", header = T, stringsAsFactors = F)

# Sort motu_asv_list based on order in asv_counts

motu_asv_list <- motu_asv_list[order(match(motu_asv_list[,1],asv_counts[,1])),]

# Bind corresponding MOTUs to the ASV count table

asv_counts <- cbind(motu_asv_list$MOTU, asv_counts, stringsAsFactors=F)
colnames(asv_counts)[1] <- "MOTU"

# Sum ASV counts based on MOTU
# First, delete the original column containing the ASV names as they cannot be summed

asv_counts <- asv_counts[-2]
motu_table <- aggregate(. ~ MOTU, data=asv_counts, FUN=sum)

# Write MOTU table to file

write.table(motu_table, file = file.path(MOTU_dir, "motu_table_COI.txt"), sep="\t", row.names=F, quote=F)
