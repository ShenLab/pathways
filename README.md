# DNVs pathway analysis
## Environment
`R>=3.6.0; package: doParallel`
## Usage
Modify your DNV table in `01.prepare.simulation.R`, see `data/06.denovo.annotated.txt` as example.
Modify the pathway file name in `01.prepare.simulation.R`, see `data/genesets/` for ready to use files.

Then Run those scripts **in terminal** step by step:
`Rscript 01.prepare.simulation.R`
`Rscript 02.simulation.R $i` **where i is the random seed to produce one simulation, you need a for loop to produce as many times of simulation as you want**
`Rscript 03.simulation.load.parallel.R` **Note this script require 48 cpu cores, modify the code if your machine does not have that many cpus**
## Output
`pathway.enrichment.results/` contains the csv formated files
`data/genesets/significant.pathways.symbols.gmt` contains the gmt file of significant pathways
`to.cytoscape.edge.tsv` and `to.cytoscape.node.tsv` can be used as input to `cytoscape` software for visualization and clustering analysis.

