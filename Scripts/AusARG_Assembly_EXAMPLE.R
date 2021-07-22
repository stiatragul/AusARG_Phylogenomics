


# STEP 1: Dedupe, Clean, Filter
dcf_shell(sample.file.path = "~/Desktop/AusARG_BPA_Download/samples.csv", no.parallel = 2,
          memory = 40, cpu = 40, sample.dir = "Optimization/", source.dir = "/home/ian/SqCL_Pipeline",
          reference = "Scripts/Squamate_AHE_UCE_genes_Pygopodidae.fasta",
          minid = 0.25, project.name = "Pygopodidae")
# ran it with parallel = 2 but I really should have spiked that way up and dropped the mem/cpu down and done 10 at a time or something


# STEP 2: Assemble with Trinity
trinity_shell(sample.file.path = "~/Desktop/AusARG_BPA_Download/samples.csv", no.parallel = 2,
              memory = 40, cpu = 40, sample.dir = "Optimization/", source.dir = "/home/ian/SqCL_Pipeline",
              project.name = "Pygopodidae")


# STEP 3: Match Contigs to Targets
match2targets_shell(sample.file.path = "~/Desktop/AusARG_BPA_Download/samples.csv",
                    no.parallel = 10, sample.dir = "Optimization/", source.dir = "/home/ian/SqCL_Pipeline",
                    project.name = "Pygopodidae", evalue = 1e-30,
                    targets.file.path = "Scripts/squamate_AHE_UCE_genes_unique.flipped.fasta")


# STEP 4: Generate Pseudo-Reference Genomes
prg_shell(sample.file.path = "~/Desktop/AusARG_BPA_Download/samples.csv",
          no.parallel = 8, sample.dir = "Optimization/", source.dir = "/home/ian/SqCL_Pipeline",
          project.name = "Pygopodidae", keep = "easy_recip_match")


# STEP 5: Get Quality Scores
qc2_shell(sample.file.path = "~/Desktop/AusARG_BPA_Download/samples.csv",
          sample.dir = "Optimization/", source.dir = "/home/ian/SqCL_Pipeline",
          project.name = "Pygopodidae", output.dir = "quality_scores",
          sample.or.lineage = "lineage")


# STEP 6: Generate Rough Alignments
# this step can be executed exclusively by the 'phylogeny_make_alignments.py' script

# STEP 7: Trim Alignments and Build Gene Trees
# this step can be executed exclusively by the 'phylogeny_align_genetrees_parallel.py' script