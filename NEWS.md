# multibridge 1.3.0

Removed explicit C++11/C++14 requirement from SystemRequirements to comply with upcoming CRAN policy changes (default C++17 is now used).

This addresses the CRAN NOTE regarding deprecated C++ standard specifications and ensures future compatibility.

# multibridge 1.2.0

The null model now constrains all category proportions to be equal. In previous versions, the null model constrained only the category proportions to be equal that were included in the informed hypothesis, leaving the other category proportions free to vary.

The output of the S3 method summary now included the log marginal likelihoods of the encompassing model, the null model, and the informed model.

# multibridge 1.1.0

In binom_bf_equality, the specification of a predicted value p is supported for all binomial probabilities.

# multibridge 1.0.0

Initial release.
