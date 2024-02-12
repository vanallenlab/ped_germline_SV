# Gameplan

Coming back to devise a formulated gameplan after some discussion with Riaz and Ryan. See my plan that I sent them in the slack below:

## 2023-12-03

_CWAS RNA_
* Need to try to convert the test to a normal (i.e. convert expression to (0, 1) normal and extract the samples with SVs

Going to start breaking out the gameplan into dates that I semi-regularly update. The callset is now v2.5.2. We've started moving away from the targeted SV enrichment in disease relative genes. While this analysis is now accomplished by `targeted-gene-preprocessing-v2.5.2.ipynb`, there wasn't a whole that was substantial or interesting. As a result, I've moved some of the disease specific features (EWSR1-FLI1 joint SV analysis) into `old`. I should probably revisit these before the project is done. We're now moving into digesting Riaz's CWAS results and Ryan's de novo sequencing work.

Before I do move on, I'd like to document lingering questions in the gene-specific space, in case I ever come back to it.

## 2023-02-11

We've advanced these approaches significantly. In total, we've essentially accomplished (or are currently accomplishing) 4 things in the scope of this paper:

1. Individual gene analyses. There is nothing of relevance here (the common results are almost assuredly artifacts), but it was important to do.
2. RNA expression analyses. In this set of analyses, we've discovered that A) coding SVs tend to decrease expression in tumors, and B) that SVs in our categories have a higher effect on expression than random coding SVs. Pretty neat!
3. Gene set enrichment analysis. In this set of analyses, we've identified gene sets that are enriched in our significant categories more than we would expect by chance. Doing these analyses right was a _massive_ headache.
4. Visualizing category hierarchies. This wasn't so much analysis as it was tinkering around with visualizations, but I think I'm approaching some figure designs that communicate what we want.

# Remaining questions

As of our most recent update, there are a few loose ends that I need to tie up:

__Really tidy up the RNA__
* We now have the St Jude RNA, which I think represents a really nice chance to show that the effect on RNA is preserved across contexts. But... the formatting for this stuff is trash. It was generated using HTSeq, which doesn't output TPMs, and I'm really reluctant to calculate TPMs myself.
* I'm finding that several things don't replicate (or are unusually significant, like noncoding analyses).
* Additionally, I was forced to subanalyze GMKF, which threw off some results. Overall, I hate this.
* I want to standardize this whole thing. There are two ways to do that:
    * Use counts alone. Normalized counts should be good enough for rank-based analyses. But this requires me to have the full counts for St Jude's, which we don't yet have. Once I have that, I can try to match up count names and just normalize. I've actually tried this approach quickly, and it's... fine. I'd prefer TPMs because that is what is familar to people.
    * More preferably--process everything from scratch myself. This would allow us to just use RSEM and be done with it. But we'd need the raw bams (which would be a pain to track down).
* Additionally, I could do some more aggressive filtering on the SVs that are included and on the genes themselves.
    * We could exclude coding SVs that don't directly overlap an exon
    * We could exclude genes that are below some expression level, reasoning that low-expressed genes have noise that prevents detection.
    
__Make better visualizations for the gene set enrichment__

* Our current visualizations for this are complete trash. I need to come up with better visualizations (and definitely visualizations that confirm case vs. control enrichment).
* This will require quite a bit of work

__Decide on an approach to visualize categories__

* This shouldn't be too hard, I just need to coordinate with Ryan/Riaz and decide on something.
        

__Other__

* Add in single gene analyses to the paper