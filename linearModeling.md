# Linear Modeling Using `limma`

The `limma` package was developed for linear modeling of gene expression data.

The model calculations in this application use `limma` functions. To see a comprehensive use manual for `limma` follow [this link](https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf). The following will be a more concise explanation based on the functionality that is used in Omics Notebook.

### Generating the Fit

In order to start with developing a model there must be data file, annotation file, data type, and data format options selected in the `Data` tab (along with any other optional inputs on that page such as normalization method).

1.  `Model Factors` [REQUIRED]: Input *at least two* groups to include as factors in the model.
2.  `Covariates` [OPTIONAL]: Check the box to add up to *three* covariates to the model. Covariates are added by column, so any covariates you would like to include must be included in the annotation file and have values for every sample. Adding a variable as a covariate helps to account for its effect on the model. It does not make the variable an option for contrasts.
3.  `Time Series` [OPTIONAL]: Check the box to add time as an interaction term in the model. If checked you must select which column contains the numeric time point information, as well as whether to include time as a *discrete* or *continuous* variable. If *discrete*, specific time points will need to be selected. If *continuous*, you will instead need to specify the number of time points to include. Unlike covariates, time as an interaction term *does* impact the contrasts available for selection.
4.  `Contrasts` [OPTIONAL]: Check the box to perform modeling on contrasts. Contrasts are used to compare two groups to one another. Any group selected as a model factor will be included in generating contrast options. If time is included as an interaction term, contrasts will be generated using each group at each time point.

Finally, when the model is generated, by clicking `Generate Fit`, a results table will render which indicates the number of up-regulated and down-regulated rows for each group in the model. The standard cutoff for significance is an adjusted p-value less than or equal to 0.05. The two numeric inputs above the results table allow for customization of the cutoff values for significance, both for adjusted p-value and logFC.

A `limma` top table will also be generated when the model is fit, which contains values such as log fold change and adjusted p-value for each row. Inputs above this table allow for selection of which contrast is being viewed.

> Note: The analysis MUST be run by clicking `Generate Fit` in order for the linear model to be created, which will then be used for the volcano plots and differential heatmaps. If any changes are made, the model must be re-run for those changes to be saved and applied to the plots.


