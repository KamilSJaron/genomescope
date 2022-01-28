# GenomeScope assumption exploration

Genomescope is awesome, but lately I run into problems where GenomeScope was unable to do what I wanted to. Why? Well, because not all sequenced genomes follow the assumption of GenomScope.

This repo will deal with very specific problem, if you just want a genome model, I would recommend to look at the [GenomeScope 2.0](https://github.com/tbenavi1/genomescope2.0) or the original [GenomeScope](https://github.com/schatzlab/genomescope) instead.

This fork is of the older GenomScope version for one reason - it's a lot simpler code to tweak, but I will certainly use some of the tricks from GenomeScope 2.0 too (and possibly completely switch to a GenomeScope 2.0 fork). Also, this is a construction zone, I will be pushing commits with my progress, not everything will be necesarily working.

## List of problems

- uneven spacing of k-mer peaks
- the genome has a haploid portion (e.g. genome of the heterogametic sex)

## References

- Vurture, GW, Sedlazeck, FJ, Nattestad, M, Underwood, CJ, Fang, H, Gurtowski, J, Schatz, MC (2017) *Bioinformatics* doi: https://doi.org/10.1093/bioinformatics/btx153
- Ranallo-Benavidez TR, Jaron KS, Schatz MC. GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. Nature Communications. 2020 Mar 18;11(1):1-0.

