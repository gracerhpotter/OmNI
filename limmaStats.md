---
editor_options: 
  markdown: 
    wrap: 72
---

## `limma` Statistics for Differential Expression

> This excerpt is pulled directly from the `limma` User Manual, Chapter
> 13 (p 61-64). To see the rest of the manual and associated
> references/bibliography [click
> here](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf).

### Summary Top-Tables

`limma` provides functions `topTable()` and `decideTests()` which
summarize the results of the linear model, perform hypothesis tests and
adjust the p-values for multiple testing. Results include (log2) fold
changes, standard errors, t-statistics and p-values. The basic statistic
used for significance analysis is the *moderated t-statistic*, which is
computed for each probe and for each contrast. This has the same
interpretation as an ordinary t-statistic except that the standard
errors have been moderated across genes, i.e., squeezed towards a common
value, using a simple Bayesian model. This has the effect of borrowing
information from the ensemble of genes to aid with inference about each
individual gene [[35](https://doi.org/10.2202/1544-6115.1027), 
[21](https://doi.org/10.1214/16-aoas920)]. Moderated t-statistics 
lead to p-values in the same way that ordinary t-statistics do except that the degrees of
freedom are increased, reflecting the greater reliability associated
with the smoothed standard errors. The effectiveness of the moderated t
approach has been demonstrated on test data sets for which the
differential expression status of each probe is known 
[[11](https://link.springer.com/book/10.1007/0-387-29362-0)].

A number of summary statistics are presented by `topTable()` for the top
genes and the selected contrast. The `logFC` column gives the value of the
contrast. Usually this represents a log2-fold change between two or more
experimental conditions although sometimes it represents a
log2-expression level. The `AveExpr` column gives the average
log2-expression level for that gene across all the arrays and channels
in the experiment. Column `t` is the moderated t-statistic. Column `P.Value`
is the associated p-value and `adj.P.Value` is the p-value adjusted for
multiple testing. The most popular form of adjustment is "BH" which is
Benjamini and Hochberg’s method to control the false discovery rate [[1](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)].
The adjusted values are often called q-values if the intention is to
control or estimate the false discovery rate. The meaning of "BH"
q-values is as follows. If all genes with q-value below a threshold, say
0.05, are selected as differentially expressed, then the expected
proportion of false discoveries in the selected group is controlled to
be less than the threshold value, in this case 5%. This procedure is
equivalent to the procedure of Benjamini and Hochberg although the
original paper did not formulate the method in terms of adjusted
p-values. 

The B-statistic (`lods` or `B`) is the log-odds that the gene is
differentially expressed [[35](https://doi.org/10.2202/1544-6115.1027), 
Section 5]. Suppose for example that B = 1.5. The odds of differential 
expression is exp(1.5) = 4.48, i.e., about
four and a half to one. The probability that the gene is differentially
expressed is 4.48 / (1 + 4.48) = 0.82, i.e., the probability is about 82% that
this gene is differentially expressed. A B-statistic of zero corresponds
to a 50-50 chance that the gene is differentially expressed. The
B-statistic is automatically adjusted for multiple testing by assuming
that 1% of the genes, or some other percentage specified by the user
in the call to `eBayes()`, are expected to be differentially expressed.
The p-values and B-statistics will normally rank genes in the same
order. In fact, if the data contains no missing values or quality
weights, then the order will be precisely the same. 

As with all model-based methods, the p-values depend on normality and other
mathematical assumptions which are never exactly true for microarray
data. It has been argued that the p-values are useful for ranking genes
even in the presence of large deviations from the assumptions [[34](https://doi.org/10.1385/1-59259-364-x:111), 
[32](https://doi.org/10.1093/bioinformatics/bti270)].
Benjamini and Hochberg’s control of the false discovery rate assumes
independence between genes, although Reiner et al [[24](https://doi.org/10.1093/bioinformatics/btf877)] have argued that
it works for many forms of dependence as well. The B-statistic
probabilities depend on the same assumptions but require in addition a
prior guess for the proportion of differentially expressed probes. The
p-values may be preferred to the B-statistics because they do not
require this prior knowledge. 

The `eBayes()` function computes one more useful statistic. The moderated 
F-statistic (`F`) combines the t-statistics for all the contrasts into 
an overall test of significance for that gene. The F-statistic tests whether 
any of the contrasts are non-zero for that gene, i.e., whether that gene 
is differentially expressed on any contrast. The denominator degrees of 
freedom is the same as that of the moderated-t. Its p-value is stored 
as `fit$F.p.value`.  It is similar to the ordinary F-statistic from analysis 
of variance except that the denominator mean squares are moderated across genes. 

A frequently asked question relates to the occasional occurrence that all 
of the adjusted p-values are equal to 1. This is not an error situation but 
rather an indication that there is no evidence of differential expression in 
the data after adjusting for multiple testing. This can occur even though 
many of the raw p-values may seem highly significant when taken as individual 
values. This situation typically occurs when none of the raw p-values are less 
than 1/G, where G is the number of probes included in the fit. In that case the 
adjusted p-values are typically equal to 1 using any of the adjustment methods 
except for `adjust = "none"`. 

