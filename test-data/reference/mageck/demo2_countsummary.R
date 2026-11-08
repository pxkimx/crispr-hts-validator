Sweave("demo2_countsummary.Rnw");
library(tools);

texi2dvi("demo2_countsummary.tex",pdf=TRUE);

