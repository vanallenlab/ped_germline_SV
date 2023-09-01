## Change logs for study callset versions

#### Version 2.2 (release date TBD)
- Fixed myriad small formatting/annotation issues  
- Excluded an extra 194 samples that were outliers for rare TSS duplications  
- Flagged variants as non-PASS if they had >20% coverage by loci with alternative haplotypes/contigs in hg38  
- Fixed no-call rate assignment on chrX & chrY to account for sex differences in genotyping  
- Manually reviewed read depth evidence for all 176 rare (AF<1%) CNVs >1Mb, resulting in 27/176 being tagged as non-PASS  
- Manually reviewed read depth evidence for CNV segments >5kb involved in all 47 rare (AF<1%) complex SVs >1Mb, resulting in 16/47 being tagged as non-PASS  

#### Version 2.1 (release date July 26, 2023)
- Added new FILTER tag to mark predicted gene retroduplication splice junctions  
- Masked RD-only GTs in ~5% of samples that were determined to be outliers for rare RD-only events
- Added new non-PASS FILTER tag for RD-only records that had ≥50% of all non-ref GTs contributed by rare RD-only outliers
- Added final layer of outlier exclusion after all site-level filtering targeted at rare exonic deletions (≥30 per sample) & rare copy-gain duplications (≥15 per sample)

#### Version 2.0 (July 20, 2023)
- Major restructuring of study and callset
- Added 2,842 new ancestry- & sex-matched adult controls from TOPMed-BioMe cohort
- Re-batched & joint-called all 10,590 samples in a single GATK-SV cohort
- Forced uniform proportions of cases & controls in each GATK-SV batch
- Improved genotype filtering with GATK GQRecalibrator trained on All of Us matching short- & long-read WGS data daisy-chained to trio-based minGQ model targeting 5% false discovery rate
- Added methods to identify and correct batch- or cohort-specific variants appearing at appreciable (>50%) frequencies
- Dropped 1000 Genomes samples from study for all downstream disease association analyses
- Updated gene reference file from MANE v0.95 to MANE v1.2
- Added variant frequency annotations from gnomAD-SV v2.1
- Incremented version of PedSV R package to v0.0.2
- Inplemented fallback to Firth's bias-reduced logistic regression in `PedSV::pedsv.glm`

#### Version 1.1 (April 7, 2023)  
- Dropped empty records (i.e., those with no non-ref samples remaining)  
- Corrected nomenclature for mCNV frequencies in VCF INFO / BED header  

#### Version 1 (April 6, 2023)  
- Initial VCFs produced by RLC  
