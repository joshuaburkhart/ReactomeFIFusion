# Gene Fusion Network Analysis

This directory contains the cleaning of gene fusion data from the [Catalogue of
Somatic Mutations in Cancer (COMSIC)][1] and basic statistical and network analysis
on the resulting data.

[1]: http://cancer.sanger.ac.uk/cosmic

## Written Scripts

- `bin/tsv-to-csv.bash`
    - Run from `bin/` directory to convert the raw data of the COSMIC `.tsv`
      file into a `.csv` file for the analysis.
- `report/cosmic_fusion_extraction.Rmd`
    - Cleans COSMIC data
    - Write file `results/newDescription.txt` with the gene fusion data
    - Takes in Ensembl and functional interaction network data for analysis
    - Performs simple statistics and basic network statistics

## Notes on Downloaded Files

- `raw-data/FIsInGene_121013_with_annotations.txt`
    - Description:
        - This is the functional interaction network file from the Reactome
          website.  The file will be used so that the gene fusions can be
          overlaid on top to calculate metrics that will categorize the
          cancerous gene fusions.
        - “Functional interactions (FIs) derived from Reactome, and other
          pathway and interaction databases.” We downloaded the Version 2013.
    - Source: http://www.reactome.org/pages/download-data/
- `bin/ensembl_GRCh37_BioMart_2014.08.29.pl`
    - Description:
        - This is a Perl script that Ensembl automatically generated, based on
          the parameters I set:
            - Associated Gene name
            - Ensembl Transcript ID
            - 5' UTR Start
            - 5' UTR End
            - Exon Chr Start (bp)
            - Exon Chr End (bp)
            - 3' UTR Start
            - 3' UTR End
            - Strand (directionality)
            - Exon Rank in Transcript (which exon number it is)
    - Downloads: `raw-data/ensembl_GRCh37_BioMart_2014.08.29.csv`
    - Source: http://grch37.ensembl.org/biomart/martview/
- `raw-data/CosmicFusionExport_v69_310514.tsv`
    - Description:
        - "All gene fusion mutation data from the current release in a tab
          separated file."
    - Source: cancer.sanger.ac.uk/cancergenome/projects/cosmic/download

## NOTE: Error in COMIC Data

If you try and compile and run the `cosmic_fusion_extraction.Rmd` analysis, it
will not fail but it will be if you start with the raw data. The reason being is
that the COSMIC gene fusion data set (`CosmicFusionExport_v69_310514.tsv`) is
missing an open bracket on Line 11620.

The converted `.csv` version of the data (`CosmicFusionExport_v69_310514.csv`)
included in this analysis is manually edited so that the `.Rmd` analysis file
will run correctly.

## Analysis Directory Structure

```
.
├── README.md
├── bin
│   ├── ensembl_GRCh37_BioMart_2014.08.29.pl
│   └── tsv-to-csv.bash
├── data
│   └── CosmicFusionExport_v69_310514.csv
├── papers
│   ├── wang2009-an-integrative-approach-to-reveal-driver-gene-fusions.pdf
│   └── wu2013-identification-of-cancer-fusion-drivers.pdf
├── raw-data
│   ├── CosmicFusionExport_v69_310514.tsv
│   ├── FIsInGene_121013_with_annotations.txt
│   └── ensembl_GRCh37_2014.08.29.csv
├── report
│   ├── cosmic_fusion_extraction.Rmd
│   └── cosmic_fusion_extraction.html
└── results
    └── newDescription.txt

6 directories, 12 files
```
