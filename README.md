# UniprotToNtFasta

This script takes a list of UniProt protein IDs and creates a fasta file with all nucleotide sequences it can find for those proteins' coding DNA sequence (CDS). It also creates several other files which were relevant to my use case, which was to create gene trees with Spimap (Rasmussen & Kellis et al. 2011)

I created it in the spring of 2025 as part of a rotation in Enzo Carnevale's lab at Temple University.

## Data Sources

The data sources of this version are the European Nucleotide Archive and NCBI's Refseq database (O'Leary NA et al. 2016).

UniProt IDs are matched to ENA and Refseq records via UniProt's ID mapping API (Huang et al. 2011).

The script gets CDS data from ENA first, and then uses Refseq for everything not picked up from ENA.  I chose this because ENA has a data source with cleanly annotated CDS sequences with their own separate accessions. Refseq requires either pasting together the CDS regions in eukaryotes, parsing text from annotations, or making a second API call to get well-formatted species data. I chose to make the second API call.

As of the first push to Github, if the sequence is picked up in ENA but filtered for being full of N values, it won't try for it again in refseq. I intend to fix this, but if this line is still in the ReadMe, that means I didn't get to it. My bad.

## Outputs
- A fasta file with a fasta entry for every CDS
- A list of all species which returned a CDS
- A species mapping file (something specifically for spimap)

## How to use
Currently, you'll have to add your own Uniprot ID string and modify the script to use those files. Input and output names are at the top. As you can tell, this was for personal use. I will change this if I start passing it around the lab or if my shame at having been a professional once catches up to me.

## Known issues
Some sequences (such as the protein M1CEV4 in test_list.txt) do not start with ATG, and are not from species that could be an exception. This is something I'm working on, so I wouldn't recommend using this seriously right now.

The aforementioned mediocre sequence handling for ENA data might mean some seqs get missed.

The API calls take a long time, and currently it doesn't make files halfway when it probably could do.


## Citations
Huang H, McGarvey PB, Suzek BE, et al. A comprehensive protein-centric ID mapping service for molecular data integration. Bioinformatics. 2011;27(8):1190-1191. doi:10.1093/bioinformatics/btr101

Matthew D. Rasmussen, Manolis Kellis, A Bayesian Approach for Fast and Accurate Gene Tree Reconstruction, Molecular Biology and Evolution, Volume 28, Issue 1, January 2011, Pages 273–290, https://doi.org/10.1093/molbev/msq189

O'Leary NA, Wright MW, Brister JR, et al. Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucleic Acids Res. 2016;44(D1):D733-D745. doi:10.1093/nar/gkv1189