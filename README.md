# AusARG Phylogenomics Working Group

The [***Aus***tralian ***A***mphibian and ***R***eptile ***G***enomics](https://ausargenomics.com/) initiative is a national collaborative project aiming to facilitate the development of genomics resources for Australia's unique amphibian and reptile fauna. The ***Phylogenomics Working Group*** operates under AusARG with the goal of collecting a consistent set of phylogenomic data for all of Australia's frogs and reptiles, including recognized subspecies and cryptic lineages. We've decided to use the ***Sq***uamate ***C***onserverd ***L***oci (***SqCL***) developed by Sonal Singhal and colleages [see the paper](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.12681), which incorporates roughly 5,000 Ultra Conserved Elements (UCEs), 400 Anchored Hybrid Enrichment targets, and 40 "legacy" exons commonly used in herpetological phylogenetics. 

# What's in this Repository?

This repo holds the scripts and software to go from raw sequence reads provided by [BioPlatforms Australia](https://bioplatforms.com/) to individual gene trees. More information about the background of this workflow and its use can be found in the ***AusARG Phylogenomics Workflow*** document. In a nutshell:

+ The *Generate_Metadata_OZCAM* file provides a tool for generating a BPA consistent metadata file which is necessary for sequence submission. 

+ The *Scripts* directory holds all necessary sequence cleaning, assembly, alignment, and tree building python scripts. This builds heavily on the [SqCL Pipeline](https://github.com/singhal/SqCL) developed by [Sonal Singhal](https://scholar.google.com.au/citations?user=hGRmhQkAAAAJ&hl=en&oi=ao).

+ The *AusARG_Phylogenomics_Workflow* document outlines how to use the tools provided here. 

+ The *Software* text file provides a link to download all the required software files. You can also access the folder by clicking [here](https://drive.google.com/drive/folders/1wb7OgU4nnvpd-RPT7XZHkp4ewABX7IRS?usp=sharing). 


# Where's the Beef?

The moving pieces of this workflow are a series of bioinformatics and phylogenetics software packages that have been containerized using [Singularity](https://sylabs.io/) (similar to *Docker*). This means the packages are free standing (including dependencies), so as long as you have *Singularity v.3+* installed on your machine you can run them without any prior installation.

While the initial steps of this workflow take place on your local machine (desktop, laptop, whatever), you'll want to execute most of the heavy lifting on a more powerful machine or server. That **analysis** machine will also need to have [Python](https://www.python.org/downloads/) installed (including *pandas* and *numpy*), and to take advantage of parallel processing you'll want [GNU Parallel](https://www.gnu.org/software/parallel/) installed as well.

## What about the mtDNA?

If you're looking to put your off-target mitochondrial reads to use, try assembling mitogenomes using the tools in [mitoGenome Assembly](https://github.com/IanGBrennan/mitoGenome_Assembly).

Good luck!