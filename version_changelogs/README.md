## Change logs for study callset versions

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
