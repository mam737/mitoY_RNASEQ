### RAN AND DOWNLOADED BY MANISHA A. MUNASINGHE ON 10/11/19



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
.libPaths('./Rlibs')
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0002376","immune system process",4.166,4.1367,0.922,0.000,"immune system process"),
c("GO:0005975","carbohydrate metabolic process",3.683,5.6198,0.818,0.000,"carbohydrate metabolism"),
c("GO:0006458","'de novo' protein folding",0.158,2.1451,0.899,0.000,"primede novoprime protein folding"),
c("GO:0008152","metabolic process",52.035,3.6882,0.960,0.000,"metabolism"),
c("GO:0019731","antibacterial humoral response",0.369,5.2807,0.672,0.000,"antibacterial humoral response"),
c("GO:0006952","defense response",3.533,3.8125,0.757,0.436,"antibacterial humoral response"),
c("GO:0006457","protein folding",1.099,1.5834,0.897,0.028,"protein folding"),
c("GO:0009113","purine nucleobase biosynthetic process",0.044,2.8013,0.557,0.082,"purine nucleobase biosynthesis"),
c("GO:0006730","one-carbon metabolic process",0.211,1.7773,0.631,0.369,"purine nucleobase biosynthesis"),
c("GO:0033499","galactose catabolic process via UDP-galactose",0.035,1.4342,0.600,0.667,"purine nucleobase biosynthesis"),
c("GO:0006030","chitin metabolic process",0.984,1.3947,0.744,0.320,"purine nucleobase biosynthesis"),
c("GO:0035999","tetrahydrofolate interconversion",0.053,1.4342,0.615,0.497,"purine nucleobase biosynthesis"),
c("GO:0006629","lipid metabolic process",4.140,2.1451,0.633,0.255,"purine nucleobase biosynthesis"),
c("GO:0006635","fatty acid beta-oxidation",0.272,2.1451,0.536,0.327,"purine nucleobase biosynthesis"),
c("GO:0055114","oxidation-reduction process",6.443,9.3605,0.628,0.469,"purine nucleobase biosynthesis"),
c("GO:0016042","lipid catabolic process",1.037,1.7852,0.592,0.685,"purine nucleobase biosynthesis"),
c("GO:0000023","maltose metabolic process",0.114,1.9788,0.742,0.279,"purine nucleobase biosynthesis"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );
levels(stuff$representative)[levels(stuff$representative)=='purine nucleobase biosynthesis'] <- 'purine biosynthesis'
levels(stuff$representative)[levels(stuff$representative)=='primede novoprime protein folding'] <- 'de novo protein folding'



pdf('./mito_REVIGOBP_treemap.pdf',width=16,height=9)
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
	palette = c("#F1EEF6","#D0D1E6","#A6BDDB","#74A9CF","#3690C0","#0570B0","#034E7B"),
	fontsize.labels=c(30,20)
)
dev.off()
