### RAN AND DOWNLOADED BY MANISHA A. MUNASINGHE ON 2/26/20
### Visualization updated as original function outputed by REVIGO was deprecated


# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
lib.path()
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0002376","immune system process",4.166,5.5950,0.922,0.000,"immune system process"),
c("GO:0005975","carbohydrate metabolic process",3.683,5.2275,0.817,0.000,"carbohydrate metabolism"),
c("GO:0006457","protein folding",1.099,1.3466,0.899,0.000,"protein folding"),
c("GO:0008152","metabolic process",52.035,4.5925,0.960,0.000,"metabolism"),
c("GO:0019731","antibacterial humoral response",0.369,5.2872,0.626,0.000,"antibacterial humoral response"),
c("GO:0006952","defense response",3.533,3.8617,0.740,0.436,"antibacterial humoral response"),
c("GO:0009113","purine nucleobase biosynthetic process",0.044,2.8314,0.563,0.082,"purine nucleobase biosynthesis"),
c("GO:0006629","lipid metabolic process",4.140,2.1183,0.631,0.255,"purine nucleobase biosynthesis"),
c("GO:0000023","maltose metabolic process",0.114,1.9961,0.743,0.279,"purine nucleobase biosynthesis"),
c("GO:0006030","chitin metabolic process",0.984,1.4001,0.744,0.320,"purine nucleobase biosynthesis"),
c("GO:0006730","one-carbon metabolic process",0.211,1.7809,0.631,0.321,"purine nucleobase biosynthesis"),
c("GO:0035999","tetrahydrofolate interconversion",0.053,1.4379,0.615,0.325,"purine nucleobase biosynthesis"),
c("GO:0016042","lipid catabolic process",1.037,1.7901,0.588,0.359,"purine nucleobase biosynthesis"),
c("GO:0055114","oxidation-reduction process",6.443,9.7027,0.625,0.469,"purine nucleobase biosynthesis"),
c("GO:0033499","galactose catabolic process via UDP-galactose",0.035,1.4379,0.599,0.573,"purine nucleobase biosynthesis"),
c("GO:0006631","fatty acid metabolic process",0.853,1.5030,0.539,0.654,"purine nucleobase biosynthesis"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

#Everything below this line was added to REVIGO output for visualization
levels(stuff$representative)[levels(stuff$representative)=='purine nucleobase biosynthesis'] <- 'purine biosynthesis'
levels(stuff$representative)[levels(stuff$representative)=='primede novoprime protein folding'] <- 'de novo protein folding'



pdf('./Rplots/final/mito_REVIGO_treemap.pdf',width=16,height=9)
treemap(
	stuff,
	index=c("representative","description"),
	vSize = "abslog10pvalue",
	vColor = "representative",
	type = "categorical",
	title = "Gene Ontology TreeMap",
	inflate.labels = FALSE,
	lowerbound.cex.labels = 0,
	bg.labels = "#CCCCCCAA",
	position.legend = "none",
	palette = c("#D0D1E6","#A6BDDB","#74A9CF","#3690C0","#0570B0","#034E7B"),
	fontsize.labels=c(30,20)
)
dev.off()
