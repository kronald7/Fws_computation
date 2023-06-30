# whithin host diversity (Fws) computation with moimix

# set working directory
setwd('/Users/ronky/Projects/WGS_SANRUHRP/updated_vcf')

# load packages

library(moimix)
library(SeqArray)

# upload the vcf and convert it to gds

seqVCF2GDS('wgs.random.309snps.vcf.gz', 
           'wgs.random.309snps.gds')

# to read the gds file in R, use the seqOpen function

isolates <- seqOpen('wgs.random.309snps.gds')
seqSummary(isolates)

# save sample identifiers
sample.id <- seqGetData(isolates, "sample.id")

# get genomic coordinates of all variants
coords <- getCoordinates(isolates)
head(coords, 10)

# filter variant.ids not on apicoplast
seqSetFilter(isolates, 
             variant.id = coords$variant.id[coords$chromosome != "Pf3D7_API_v3"])

#### ESTIMATE THE BAF MATRIX #####

isolate_baf <- bafMatrix(isolates)
class(isolate_baf)

str(isolate_baf)

plot(isolate_baf, "112037")

### Estimating Fws ###
fws_all <- getFws(isolates)
fws.df <- as.data.frame(fws_all)
#hist(fws_all)

### create a final dataframe containing the sample id, the geographic origin
### and fws values

# uploda a panel file containing sample id and origin in a dictionary style

panel_file <- read.table('panel_file.txt', header = F, sep = '\t')
variables = c('sample.id', 'province')
panel_file <- `colnames<-`(panel_file, variables)
head(panel_file)

### build up the dataframe
fws_long.panel <- as.data.frame(panel_file |> dplyr::select(2,1) |> 
  cbind(fws.df))
head(fws_long.panel)

### visualize the fws per province
library(ggplot2)

fws_long.panel %>% 
  ggplot() +
  geom_histogram(aes(x=fws_all, y = (..count../sum(..count..))*100),
                 color = "#000000", fill = "#d9d9d9") +
  facet_wrap(~province) +
  xlab("Fws") + ylab("Frequency (%)") +
  (labs(title = 'Within sample diversity (Fws)'))+
  theme_linedraw()
rownames(fws_long.panel) <- NULL

### Get monoclonal samlples for PCA

# monoclonal infections (Fws >= 0.95)
monoclo_samples <- as.data.frame(fws_long.panel %>% filter(fws_all >= 0.95))
rownames(monoclo_samples) <- NULL
head(monoclo_samples)

### visualize the proportion of polyclonal sample vs monoclo

### write the monoclonal samples id in a text file for PCA
monoclo_samples %>% select(sample.id, province) %>% 
  write.table("/Users/ronky/Projects/WGS_SANRUHRP/updated_vcf/pf.wgs.panel.txt",
              sep = "\t",
              row.names = FALSE)
