# The `QR` argument

Details about the `QR` argument to rstanarm's modeling functions.

## Details

The `QR` argument is a logical scalar defaulting to `FALSE`, but if
`TRUE` applies a scaled [`qr`](https://rdrr.io/r/base/qr.html)
decomposition to the design matrix, \\X = Q^\ast R^\ast\\. If
`autoscale = TRUE` (the default) in the call to the function passed to
the `prior` argument, then \\Q^\ast = Q \sqrt{n-1}\\ and \\R^\ast =
\frac{1}{\sqrt{n-1}} R\\. When `autoscale = FALSE`, \\R\\ is scaled such
that the lower-right element of \\R^\ast\\ is \\1\\.

The coefficients relative to \\Q^\ast\\ are obtained and then
premultiplied by the inverse of \\R^{\ast}\\ to obtain coefficients
relative to the original predictors, \\X\\. Thus, when
`autoscale = FALSE`, the coefficient on the last column of \\X\\ is the
same as the coefficient on the last column of \\Q^\ast\\.

These transformations do not change the likelihood of the data but are
recommended for computational reasons when there are multiple
predictors. Importantly, while the columns of \\X\\ are almost generally
correlated, the columns of \\Q^\ast\\ are uncorrelated by design, which
often makes sampling from the posterior easier. However, because when
`QR` is `TRUE` the `prior` argument applies to the coefficients relative
to \\Q^\ast\\ (and those are not very interpretable), setting `QR=TRUE`
is only recommended if you do not have an informative prior for the
regression coefficients or if the only informative prior is on the last
regression coefficient (in which case you should set `autoscale = FALSE`
when specifying such priors).

For more details see the Stan case study *The QR Decomposition For
Regression Models* at
<https://mc-stan.org/users/documentation/case-studies/qr_regression.html>.

## References

Stan Development Team. *Stan Modeling Language Users Guide and Reference
Manual.* <https://mc-stan.org/users/documentation/>.
