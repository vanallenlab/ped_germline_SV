## Change logs for study callset versions

#### Version 2.5.4 (release date August 8, 2024)
- Excluded 1 sample carrying three suspicious large (500kb-6Mb) mosaic deletions on chr11, chr16, and chr17  
- Upgrade external allele frequencies to the "controls only" subset of gnomAD v4.1  

#### Version 2.5.3 (release date December 21, 2023)
- Excluded 3 samples carrying suspicious translocations resulting in known oncogenic fusions (BCR-ABL, EWSR1-FLI1, NUP98-PRRX1), one of which was unable to be validated by PCR in blood DNA (EWSR1-FLI1)  

#### Version 2.5.2 (release date November 1, 2023)
- Modified outlier sample definition to treat samples of unknown ancestry as European ancestry  
- Increased P-value cutoff for inter-cohort frequency comparisons from P<0.01 to P<0.05  
- Pooled controls from both study arms (discovery & trio) for osteosarcoma, since there are no osteosarcoma trios  

#### Version 2.5 (release date October 20, 2023)
- Re-trained secondary GT filtering method, `MinGQ`, while restricting to pediatric cancer trios (and not including any 1000G trios)  
- Increased global maximum no-call rate to 4%  
- Streamlined/reworked outlier sample exclusion and site-level refinement  
- Increased stringency when defining samples with unreliable RD-only genotypes introduced in version 2.1 (see below)  
- Removed batch- and cohort-specific artifact checks implemented in v2.0  
- Masked all genotypes with `OGQ` = 0 for all variants with AF < 5%  
- Fail sites with SL mean < 0, Manta- and/or Wham-only, <1kb, no BAF evidence, AF<5%, NCR>0.1%, no PE/SR GT overdispersion, SL max < 75  
- Corrected a bug causing 12 extra cases to be included in analyses despite not having ancestry-matched controls  

#### Version 2.4 (release date September 29, 2023)
- Masked individual genotypes with `SL` <= 1 for variants supported by just one algorithm or one evidence type.  
- Dropped `FILTER` corresponding to variants with low `SL_MAX`; it was no longer necessary given sample-level SL-based filtering above.  

#### Version 2.3 (release date September 25, 2023)
- Flagged all variants with `SL_MAX` < 0 as non-PASS  
- Compared population-specific frequencies for all relatively uncommon (AF<5%) SVs between GMKF vs. St. Jude + ICGC. Marked all variants with significant (P<0.01) discrepancies in frequencies between cohorts with `INTERCOHORT_HETEROGENEITY`.  
- Compared population-specific frequencies for all relatively uncommon (AF<5%) SVs between MESA + BioMe. Marked all variants with significant (P<0.01) discrepancies in frequencies between cohorts with `INTERCOHORT_HETEROGENEITY`.  
- Tagged all variants with >20% coverage by hg38 reference fix patches with `HG38_PATCH_LOCUS` in `INFO`; note that `FILTER` for these variants remained unchanged.  
- Increased minimum kinship coefficient cutoff from 1/16 to 0.1 (slightly more lenient than second-degree relatives).  
- Reconfigured study phase assignment for cases such that each phase only ever has cases from St. Jude or GMKF (but not both) for each disease.  
- Switched from using full gnomAD v3.1 frequencies to the non-cancer subset of gnomAD v3.1. Note that the freuqency annotation variable names have not changed, but the underlying data has.  
- Default behavior of `PedSV::load.sv.bed()` now no longer reads all population- and sex-specific frequency information into memory. Old (verbose) behavior can be restored by setting `keep.all.pop.frequencies` and `keep.all.sex.frequencies` to `TRUE`. See `?PedSV::load.sv.bed` for more info.  

#### Version 2.2.1 (release date September 7, 2023)
- Upgraded external AF annotation from gnomAD v2.1 to (pre-release/unpublished) gnomAD v3.1  

#### Version 2.2 (release date September 5, 2023)
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
