# Software and data references

The versions below match the commands in the analysis runbook. The selected
SRA Toolkit release is 3.4.1; the exact version used for the original download
was not recorded. This release includes long-read support, including reads
longer than 65 kb.

1. **NCBI BioProject PRJNA1133668.** *Carbapenem-resistant Enterobacterales
   from a New York City Hospital sequencing data.* [NCBI BioProject record](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1133668).

2. **NCBI SRA Toolkit 3.4.1** (`prefetch`, `fasterq-dump`). [NCBI SRA Toolkit
   download archive](https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.4.1/);
   [release notes](https://github.com/ncbi/sra-tools/blob/master/README.md).

3. Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ
   preprocessor. *Bioinformatics*. 2018;34(17):i884–i890.
   [doi:10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560).
   Used as fastp 1.4.0 (`staphb/fastp:1.4.0`).

4. Seemann T. Shovill: assemble bacterial isolate genomes from Illumina
   paired-end reads. [Source code](https://github.com/tseemann/shovill).
   Used as Shovill 1.4.2 (`staphb/shovill:1.4.2`).

5. Wick RR, Schultz MB, Zobel J, Holt KE. Bandage: interactive visualization
   of de novo genome assemblies. *Bioinformatics*. 2015;31(20):3350–3352.
   [doi:10.1093/bioinformatics/btv383](https://doi.org/10.1093/bioinformatics/btv383).
   Used as Bandage 0.9.0 (`staphb/bandage:0.9.0`).

6. Wick RR, Judd LM, Gorrie CL, Holt KE. Unicycler: resolving bacterial
   genome assemblies from short and long sequencing reads. *PLoS
   Computational Biology*. 2017;13(6):e1005595.
   [doi:10.1371/journal.pcbi.1005595](https://doi.org/10.1371/journal.pcbi.1005595).
   Used as Unicycler 0.5.1 (`staphb/unicycler:0.5.1`) for hybrid assembly.

7. Schwengers O, Jelonek L, Dieckmann MA, Beyvers S, Blom J, Goesmann A.
   Bakta: rapid and standardized annotation of bacterial genomes via
   alignment-free sequence identification. *Microbial Genomics*.
   2021;7(11):000685.
   [doi:10.1099/mgen.0.000685](https://doi.org/10.1099/mgen.0.000685).
   Used as Bakta 1.12.2 (`staphb/bakta:1.12.2`).

8. Robertson J, Nash JHE. MOB-suite: software tools for clustering,
   reconstruction and typing of plasmids from draft assemblies.
   *Microbial Genomics*. 2018;4(8):e000206.
   [doi:10.1099/mgen.0.000206](https://doi.org/10.1099/mgen.0.000206).
   Used as MOB-suite 3.1.9
   (`quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_1`).

9. Arredondo-Alonso S, Gladstone RA, Pöntinen AK, et al. Mge-cluster: a
   reference-free approach for typing bacterial plasmids. *NAR Genomics and
   Bioinformatics*. 2023;5(3):lqad066.
   [doi:10.1093/nargab/lqad066](https://doi.org/10.1093/nargab/lqad066).
   Used as mge-cluster 1.1.0.

10. **unitig-caller 1.3.2.** [Source code](https://github.com/bacpop/unitig-caller).
   Used by mge-cluster to identify unitigs in plasmid sequences.

11. Tonkin-Hill G, MacAlasdair N, Randon C, et al. Producing polished
    prokaryotic pangenomes with the Panaroo pipeline. *Genome Biology*.
    2020;21:180.
    [doi:10.1186/s13059-020-02090-4](https://doi.org/10.1186/s13059-020-02090-4).
    Used as Panaroo 1.8.1
    (`quay.io/biocontainers/panaroo:1.8.1--pyhdfd78af_0`).

12. Seemann T. Snippy: rapid haploid variant calling and core genome
    alignment. [Source code](https://github.com/tseemann/snippy).
    Used as Snippy 4.6.0.
