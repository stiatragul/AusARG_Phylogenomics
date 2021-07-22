# AusARG Phylogenomics Working Group

The [**Aus**tralian **A**mphibian and **R**eptile **G**enomics](https://ausargenomics.com/) initiative is a national collaborative project aiming to facilitate the development of genomics resources for Australia's unique amphibian and reptile fauna. The ***Phylogenomics Working Group*** is operating with the goal of collecting a consistent set of phylogenomic data for all of Australia's frogs and reptiles, including recognized subspecies and cryptic lineages. We've decided to use the **Sq**uamate **C**onserverd **L**oci (**SqCL**) developed by Sonal Singhal and colleages [see the paper](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.12681), which incorporates roughly 5,000 Ultra Conserved Elements (UCEs), 400 Anchored Hybrid Enrichment targets, and 40 "legacy" exons commonly used in herpetological phylogenetics. 

# What's in this Repository?

This repo holds the scripts and software to go from raw sequence reads provided by [BioPlatforms Australia](https://bioplatforms.com/) to individual gene trees. More information about the background of this workflow and its use can be found in the ***AusARG Phylogenomics Workflow*** document. 

+ The *Generate_Metadata_OZCAM* file provides a tool for generating a BPA consistent metadata file which is necessary for sequence submission. 

+ The *Scripts* directory holds all necessary sequence cleaning, assembly, alignment, and tree building python scripts. 

+ The *AusARG_Phylogenomics_Workflow* document outlines how to use the tools provided here. 

+ The *Software* text file provides a link to download all the required software files. You can also access the folder by clicking [here](https://drive.google.com/drive/folders/1wb7OgU4nnvpd-RPT7XZHkp4ewABX7IRS?usp=sharing). 

