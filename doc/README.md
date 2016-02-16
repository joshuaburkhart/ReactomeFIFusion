# ReactomePPI

## Input

http://reactomews.oicr.on.ca:8080/caBigR3WebApp2014/FIsInGene_121514_with_annotations.txt.zip

Format
```
Gene1    Gene2    Annotation    Direction    Score
A2M    APOA1    inhibit    -|    1.00
A2M    APP    predicted    -    0.61
A2M    BMP1    inhibit    -|    1.00
A2M    CDC42    catalyze    ->    1.00
A2M    CTSB    predicted    -    0.61
```

(200,000+ rows)

## Preprocessing

### MyGene

https://pypi.python.org/pypi/mygene
https://www.biostars.org/p/22/

### Install

```bash
$ pip install mygene
```

### Usage

```python
import mygene
mg = mygene.MyGeneInfo()
geneSymbols = ['ZNF888', 'ZP1', 'ZP2']
mg.querymany(geneSymbols,
             scopes="symbol",
             fields=["pdb","uniprot"],
             species="human",
             as_dataframe=True)
```

## Preprocessing Results

Python script produces 72000+ rows.  
72123 / 217250 = 0.33198158803222094  
33% coverage of FI's

Many genes have several PDB's, how can we prioritize them?

### Missing Data

http://www.sbg.bio.ic.ac.uk/phyre2/webscripts/dbmapper.cgi

If no PDB exists and time allows, get swiss-prot and use phyre2 webserver. PDB File can be downloaded. *Now what do we do with it? Upload to webserver?

## PPI

### Query Interactome3D Database for Known Interactions & Complexes

http://interactome3d.irbbarcelona.org/help.php#restful  
http://interactome3d.irbbarcelona.org/help.php#interactions_dat_file  
https://docs.python.org/3/library/xml.etree.elementtree.html  

### Use a Modelling Program, like IMP to Assess Putative Interactions & Complexes


####IMP

https://integrativemodeling.org/  
review http://www.cgl.ucsf.edu/chimera/
gather data from experimental & theoretical databases

    Experimental techniques, such as:
        X-ray crystallography
        nuclear magnetic resonance (NMR) spectroscopy (CSP, NOE, J-couplings)
        electron microscopy (EM) (2D class averages or 3D maps)
        footprinting
        Immunoprecipitation pull-down
        Cysteine cross-linking
        Chemical cross-linking
        FRET spectroscopy
        small angle X-ray scattering (SAXS)
        proteomics
    Theoretical sources of information, such as:
        template structures used in comparative/homology modeling
        scoring functions used in molecular docking
        statistical preferences
        physics-based energy functions

PDBe (Cryo EM) API
http://www.ebi.ac.uk/pdbe/pdbe-rest-api

Store entry information & rank interactions based on information available
Begin IMP analysis in order of most information -> least


####ZDOCK

Example:
mark_sur receptor.pdb receptor_m.pdb
mark_sur ligand.pdb ligand_m.pdb
zdock -R receptor_m.pdb -L ligand_m.pdb -o zdock.out
create.pl zdock.out

Please note: The file uniCHARMM should be in the current directory when
executing mark_sur. Also, receptor_m.pdb, ligand_m.pdb and create_lig must
be in your current directory when you create all predicted structures
using create.pl.

Standard PDB format files must be processed by mark_sur before being used as
the input to ZDOCK. Formatted PDB files of docking benchmark can be downloaded
at http://zlab.umassmed.edu/benchmark. If you know that some atoms
are not in the binding site, you can block them by changing their ACE type
(column 55-56) to 19. This blocking procedure can improve docking
performance significantly. A blocking script block.pl is included, type
"block.pl" for usage information.

####BUDE

http://www.bris.ac.uk/biochemistry/research/bude

### Other Python Libraries & APIs

http://bmcstructbiol.biomedcentral.com/articles/10.1186/1472-6807-9-27  
http://www.pyrosetta.org/  
https://downloads.ccdc.cam.ac.uk/documentation/API/modules/docking_api.html  
http://structure.bu.edu/content/protein-protein-docking  
http://tools.iedb.org/main/bcell/modeling-docking/

### ClusPro Webserver
http://cluspro.bu.edu/home.php

### HADDOCK Webserver
http://www.bonvinlab.org/software/haddock2.2/pdb.html

### Interpretation
