# Plot the output of [`mldest()`](https://dcgerard.github.io/ldsep/reference/mldest.md) or [`sldest()`](https://dcgerard.github.io/ldsep/reference/sldest.md) using [`corrplot()`](https://rdrr.io/pkg/corrplot/man/corrplot.html)

Formats the LD estimates in the form of a matrix and creates a heatmap
of these estimates. This heatmap is created using the
[`corrplot`](https://rdrr.io/pkg/corrplot/man/corrplot.html) R package.
I've adjusted a lot of the defaults to suit my visualization
preferences.

## Usage

``` r
# S3 method for class 'lddf'
plot(
  x,
  element = "r2",
  type = c("upper", "full", "lower"),
  method = c("color", "circle", "square", "ellipse", "number", "shade", "pie"),
  diag = FALSE,
  is.corr = NULL,
  tl.pos = "n",
  title = NULL,
  na.label = "square",
  ...
)
```

## Arguments

- x:

  An object of class `lddf`, usually created using either
  [`mldest()`](https://dcgerard.github.io/ldsep/reference/mldest.md) or
  [`sldest()`](https://dcgerard.github.io/ldsep/reference/sldest.md).

- element:

  Which element of `x` should we plot?

- type:

  Character, `"full"`, `"upper"` (default) or `"lower"`, display full
  matrix, lower triangular or upper triangular matrix.

- method:

  See [`corrplot()`](https://rdrr.io/pkg/corrplot/man/corrplot.html) for
  available options. Default value is `"color"`.

- diag:

  Logical, whether display the correlation coefficients on the principal
  diagonal.

- is.corr:

  See [`corrplot()`](https://rdrr.io/pkg/corrplot/man/corrplot.html).
  Default behavior is `TRUE` if an element is constrained between -1 and
  1 and `FALSE` otherwise.

- tl.pos:

  See [`corrplot()`](https://rdrr.io/pkg/corrplot/man/corrplot.html).
  Default value is `"n"`.

- title:

  What should the title be? Defaults to the element name.

- na.label:

  See [`corrplot()`](https://rdrr.io/pkg/corrplot/man/corrplot.html).
  Default value is `"square"`.

- ...:

  Additional arguments to pass to
  [`corrplot()`](https://rdrr.io/pkg/corrplot/man/corrplot.html). See
  the documentation of that function for options.

## Value

(Invisibly) returns a matrix of the selected elements.

## Details

For greater plotting flexibility, see
[`corrplot()`](https://rdrr.io/pkg/corrplot/man/corrplot.html) for the
parameter options.

## Author

David Gerard

## Examples

``` r
set.seed(1)

## Simulate genotypes when true correlation is 0
nloci <- 5
nind  <- 100
K <- 6
nc <- 1
genomat <- matrix(sample(0:K, nind * nloci, TRUE), nrow = nloci)

## Haplotypic LD estimates
lddf <- mldest(geno = genomat,
               K = K,
               nc = nc,
               type = "hap")

## Plot estimates of z
plot(lddf, element = "z")

```
