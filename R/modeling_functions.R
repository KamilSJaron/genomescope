# for now I will just try to create a model that will fit 1n and 2n peaks indipendently

nls_2peak <- function(x, y, estKmercov, max_iterations){

    cat("Fitting patemeters: \n")
    cat(paste("Fitting coverage range: ", min(x), max(x), '\n'))
    cat(paste("Fitting frequency range: ", min(y), max(y), '\n'))
    cat(paste0("\testKmercov:\t", estKmercov, "\n"))

    model2 <- NULL

  cat("trying nls_2peak standard algorithm\n")
    try(model2 <- nls(y ~ N * (dnbinom(x, size = kmercov_1n / bias1, mu = kmercov_1n) +
                               dnbinom(x, size = kmercov_2n / bias2, mu = kmercov_2n)),
                    start = list(kmercov_1n=estKmercov, kmercov_2n=(2 * estKmercov), bias1 = 1, bias2 = 2, N = sum(y)),
                    control = list(minFactor=1e-12, maxiter=max_iterations)))

    if(is.null(model2)){
        cat("retrying nls_2peak with port algorithm\n")
        try(model2 <- nls(y ~ N * (dnbinom(x, size = kmercov_1n / bias1, mu = kmercov_1n) +
                                         dnbinom(x, size = kmercov_2n / bias2, mu = kmercov_2n)),
                                start = list(kmercov_1n=estKmercov, kmercov_2n=(2 * estKmercov), bias1 = 1, bias2 = 2, N = sum(y)),
                          algorithm="port", control = list(minFactor=1e-12, maxiter=max_iterations)))
    }

    return(model2)
}


# for now I will just try to create a model that will fit 1n and 2n peaks indipendently
nls_2peak_conditional <- function(x, y, k, estKmercov, estLength, rEst, biasEst, max_iterations){
    model2 = NULL

    cat("Fitting patemeters: \n")
    cat(paste("Fitting coverage range: ", min(x), max(x), '\n'))
    cat(paste("Fitting frequency range: ", min(y), max(y), '\n'))
    cat(paste0("\tk:\t", k, "\n"))
    cat(paste0("\testKmercov:\t", estKmercov, "\n"))
    cat(paste0("\testLength:\t", estLength, "\n"))
    cat(paste0("\tmax_iterations:\t", max_iterations, "\n"))

  cat("trying nls_2peak_conditional standard algorithm\n")

    try(model2 <- nls(y ~ ((2*(1-(1-rEst)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
                              ((1-rEst)^k)        * dnbinom(x, size = kmercov*2 / bias, mu = kmercov2)) * estLength,
                      start = list(kmercov=estKmercov, kmercov2 = (estKmercov * 2), bias = biasEst),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    return(model2)
}


# for now I will just try to create a model that will fit 1n and 2n peaks indipendently

nls_2peak_kmer_explicit <- function(x, y, k, estKmercov, estLength, max_iterations){
    model2 = NULL

    # to do
    cat("Fitting patemeters: \n")
    cat(paste("Fitting coverage range: ", min(x), max(x), '\n'))
    cat(paste("Fitting frequency range: ", min(y), max(y), '\n'))
    cat(paste0("\tk:\t", k, "\n"))
    cat(paste0("\testKmercov:\t", estKmercov, "\n"))
    cat(paste0("\testLength:\t", estLength, "\n"))
    cat(paste0("\tmax_iterations:\t", max_iterations, "\n"))

    cat("trying nls_2peak_kmer_explicit standard algorithm\n")

    try(model2 <- nls(y ~ ((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
                          ((1-r)^k)        * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2)) * length,
                      start = list(r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    return(model2)
}

nls_2peak_male_X0 <- function(x, y, k, estKmercov, estLength, estR, max_iterations){
    model2 = NULL

    cat("trying nls_2peak_kmer_explicit standard algorithm\n")
    # fixed parameters
    r = estR
    length = estLength

    try(model2 <- nls(y ~ ((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
                          ((1-r)^k)        * dnbinom(x, size = kmercov * 2 / bias, mu = kmercov * 2)) * length * fraction_diploid +
                          dnbinom(x, size = kmercov / bias, mu = kmercov) * length * (1 - fraction_diploid),
                      start = list(kmercov=estKmercov, bias = 0.5, fraction_diploid=0.9),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    return(model2)
}

predict_disomic_portion_2peak_male_X0 <- function(x, r, k, kmercov, bias, length, fraction_diploid){
    ((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
    ((1-r)^k)        * dnbinom(x, size = kmercov * 2 / bias, mu = kmercov * 2)) * length * fraction_diploid
}

predict_monosomic_portion_2peak_male_X0 <- function(x, r, k, kmercov, bias, length, fraction_diploid){
    dnbinom(x, size = kmercov / bias, mu = kmercov) * length * (1 - fraction_diploid)
}

nls_2peak_male_XY <- function(x, y, k, estKmercov, estLength, estR, max_iterations){
    model2 = NULL

    cat("trying nls_2peak_kmer_explicit standard algorithm\n")

    # fixed parameter
    r = estR

    try(model2 <- nls(y ~ ((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
                          ((1-r)^k)        * dnbinom(x, size = kmercov * 2 / bias, mu = kmercov * 2)) * disomic_length +
                          dnbinom(x, size = kmercov / bias, mu = kmercov) * monosomic_length,
                      start = list(kmercov=estKmercov, bias = 0.5, monosomic_length = 0.1 * estLength, disomic_length = 0.9 * estLength),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    return(model2)
}