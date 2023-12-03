# Gameplan

Coming back to devise a formulated gameplan after some discussion with Riaz and Ryan. See my plan that I sent them in the slack below:

## 2023-12-03

Going to start breaking out the gameplan into dates that I semi-regularly update. The callset is now v2.5.2. We've started moving away from the targeted SV enrichment in disease relative genes. While this analysis is now accomplished by `targeted-gene-preprocessing-v2.5.2.ipynb`, there wasn't a whole that was substantial or interesting. As a result, I've moved some of the disease specific features (EWSR1-FLI1 joint SV analysis) into `old`. I should probably revisit these before the project is done. We're now moving into digesting Riaz's CWAS results and Ryan's de novo sequencing work.

Before I do move on, I'd like to document lingering questions in the gene-specific space, in case I ever come back to it.

*Remaining questions*
* Some SV burden results are pretty consistently significant. Worth thinking about whether that is worth reporting.
* There may also be some significant individual gene hits. We'll need to go back over these with a fine-toothed comb at some point with the most updated SV callset. Most of the hits, however, are common variants with very large effect sizes, which suggests that they are fake.
* More concretely examine MYCN in neuroblastoma and its overlap with somatic SVs. We're seeing that neuroblastoma has the most believable signals, so it's possible something is still here. We also have some expression data.
* Go back and look at the EWSR1-FLI1 joint SV analysis. It's such a cool analysis that we should make sure nothing comes of it.
* Cross reference any potential hits in neuroblastoma using the RNA.