

### Using simplistic 2 peak model

```{R}
input_file1 <- 'analysis/real_data/springtails/Afus1_k21_truncated.hist' # high coverage male
input_file4 <- 'analysis/real_data/springtails/WW5-3_k21_truncated.hist' # low coverage female

male_kmer_spectrum <- read.table(input_file1, col.names = c('coverage', 'frequency'))
# this a globular springtail, the is one more biolgical phenomenon that makes it suboptimal for this benchmark - it has two different karyotypes in two tissues
# one is X0, second is monoploid including all chormsomes
# I know the data, so I will "remove" the haploid portion of the sequencing run by subtracting 17x of the coverage, then I will use only 20x and more for fitting (the visual inspection will show you the peak is still there, nothing to worry about)
male_kmer_spectrum[, 'coverage'] <- male_kmer_spectrum[, 'coverage'] - 18

female_kmer_spectrum <- read.table(input_file4, col.names = c('coverage', 'frequency'))

source('R/modeling_functions.R')

x <- female_kmer_spectrum$coverage
y <- female_kmer_spectrum$frequency
k <- 21
estKmercov <- 15
estLength <- 300e6
max_iterations <- 40

female_2peak_explicit <- nls_2peak_kmer_explicit(x, y, k, estKmercov, estLength, max_iterations)
female_2peak_estimates <- coef(female_2peak_explicit)
# 3.703e-03 (full 4 peak heterozygosity was 0.00354, so rather close)
```

So we have all female parameters estimates, heterozygosity and genome size are accurate. So let's use these parameters to estimate the fraction of X / autosomes from male k-mer spectra.

```{R}
# let's subset only the part of the part kmer spectra that is interesting for us
cov_range <- 20:150

# visual inspection
# plot(male_kmer_spectrum[5:150, 'frequency'] ~ male_kmer_spectrum[5:150, 'coverage'])

x <- male_kmer_spectrum[cov_range, 'coverage']
y <- male_kmer_spectrum[cov_range, 'frequency']
estR <- female_2peak_estimates['r'] # FIXED parameter for male
estLength <- female_2peak_estimates['length'] # FIXED parameter for male, we assume X0 sex determination here
estKmercov <- 35 # from visual inspection of k-mer spectra - can be automated in future

male_model_2peaks <- nls_2peak_male_X0(x, y, k, estKmercov, estLength, estR, max_iterations)
# male_regular_model_2peaks <- nls_2peak_kmer_explicit(x, y, k, estKmercov, estLength, max_iterations)

male_parameters <- coef(male_model_2peaks)
```

Wow, the fit worked! It seems I was able to estimate the proportion of diploid tissue. As a matter of fact, that model is simpler compared to the previous one we fir to female, because the three parameters we fir are 1 k-mer coverage 2 "bias" which is just overdispersal factor (I am using genomescope code terminology although this name is a bit unfortunate). Anyway, let's see how the method have performed...

```{R}
barplot <- barplot(y ~ x, col = 'deepskyblue', border = F, xlab = 'Coverage', ylab = 'Frequency')
lines(predict(male_model_2peaks, response = T) ~ barplot, lwd = 3)

disomic_prediction <- predict_disomic_portion_2peak_male_X0(x, estR, k, male_parameters['kmercov'], male_parameters['bias'], estLength, male_parameters['fraction_diploid'])
monosomic_prediction <- predict_monosomic_portion_2peak_male_X0(x, estR, k, male_parameters['kmercov'], male_parameters['bias'], estLength, male_parameters['fraction_diploid'])

lines(disomic_prediction ~ barplot, lwd = 3, col = 'darkgoldenrod1')
lines(monosomic_prediction ~ barplot, lwd = 3, col = 'red')

legend('topright',
       c('kmer histogram','full model', 'autosomes', 'X chromosomes'),
       col = c('deepskyblue','black', 'darkgoldenrod1', 'red'),
       lty = c(NA, 1, 1, 1), lwd = 3,
       pch = c(15, NA, NA, NA), bty = 'n')

total_genome <- round(estLength / 1e6, 2)
X_chrom_size_est <- round(estLength * (1 - male_parameters['fraction_diploid']) / 1e6, 2)
title(paste0('Estimated X chromosome size: ', X_chrom_size_est, ' Mbp out of total ', total_genome, ' Mbp'))
```

![first_model_fit](https://pbs.twimg.com/media/FKStIv7WUAAT86h?format=jpg&name=large)

Oh my goodness, that looks GREAT! I think, I would like to

1. extend to model to XY - there I won't be able to model "proportions" of X and autosomes, because the haploid size of male and female is NOT the same (male is larger by the size of Y).
2. I would like to consider 4 peak model

### 1. XY model

I can use the X0 k-mer spectra and fit there the XY model and see if I get the same estimates of X and Autosomes as I have before. The difference is, that the X-linked size is no longer fit as a proportion of fixed genome size, but rather both X-linked (/ Y-linked) and Autosomal fractions are modelled independently. As a consequence, the model has one more variable to fix (one more degree of freedom), but it still might work. Let's see about that. And the initial values for the two will be 0.9 of genome size estimate and 0.1 of the genome size estimate respectively

```{R}
male_model_2peaks_XY <- nls_2peak_male_XY(x, y, k, estKmercov, estLength, estR, max_iterations)
coef(male_model_2peaks_XY)
```

This looks good, but I am done for today :-)

### 2. using 4 peak model

The first step here is to fit regular GenomeScope to the female k-mer spectra.

```{R}
source('R/regular_genomescope_functions.R')

# in this excercise we WON'T TRY to be smart about estimating initial values. First we need to get the model right, so we use model where we know sane initial values
x <- female_kmer_spectrum$coverage
y <- female_kmer_spectrum$frequency
k <- 21
estKmercov <- 15
estLength <- 300e6
max_iterations <- 40

VERBOSE=TRUE
female_model <- nls_4peak(x, y, k, estKmercov, estLength, max_iterations)
female_estimates <- coef(female_model)
```

now ...