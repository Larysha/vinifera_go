## GO enrichment for the *Vitis vinifera* T2T v5.1 assembly

This tool runs gene ontology (GO) enrichment analysis on gene lists from the latest
telomere-to-telomere *Vitis vinifera* assembly, PN40024 T2T v5.1. Gene IDs follow the
`Vitvi05_01...` format, for example `Vitvi05_01chr01g00350`.

It exists to fill a gap. [ShinyGO](https://bioinformatics.sdstate.edu/go/), the tool most
people reach for, only carries the older 12X / V1 annotation, so it cannot map current T2T
v5.1 gene IDs at all. This app does the same kind of analysis, but against the current
assembly, so you can paste in your gene list and get results without writing any code.

## What annotations are used, and where they came from

The GO annotations were downloaded from Grapedia, the *Vitis* genome resource, as the
[T2T_5.1_blast2go.zip](https://grapedia.org/wp-content/uploads/2024/07/T2T_5.1_blast2go.zip)
release, and converted into the gene-set format the app needs.

These are blast2go annotations, meaning GO terms were transferred to grapevine genes by
orthology to characterised genes in other species. Coverage sits at roughly 65% of the
genome, which is higher than the MapMan-based alternative but still incomplete: not every
gene carries a GO term. That matters for interpretation. A gene missing from the results is
not necessarily uninvolved, it may simply be unannotated.

The current annotation holds {{N_TERMS}} GO terms ({{N_BP}} biological process, {{N_MF}}
molecular function, and {{N_CC}} cellular component) across {{N_GENES}} annotated genes.

## What the app does

You give it a list of genes of interest. It tests whether any GO terms turn up in that list
more often than you would expect by chance, given a background set of genes. A term that is
over-represented points to a biological process, molecular function, or cellular location
that your genes share.

There are three ways to run the analysis:

Over-representation analysis (ORA) is the default. It takes your gene list as a single set
and tests each GO term for over-representation against the background.

Split by direction becomes available once you supply an expression or log2 fold-change
column. It runs ORA separately on the up-regulated and down-regulated genes, which is useful
when the two directions tell different biological stories.

GSEA (gene set enrichment analysis) also uses the expression column, but instead of splitting
on a threshold it ranks every gene by its value and looks for GO terms whose genes sit
disproportionately towards the top or bottom of that ranking. It picks up coordinated, subtle
shifts that a hard cut-off would miss.

## Getting started

Upload your gene list under Inputs. It can be a CSV, TSV, or plain text file, either one gene
per line or a table with a gene column and, optionally, a numeric expression / log2 fold-change
column. Lines starting with `#` are treated as comments and ignored.

The app tries to detect your gene column and any expression column, but you can override its
choice. Once a file is loaded it shows a quick input check: how many genes you provided, how
many match the expected T2T v5.1 ID format, and how many are present in the annotation. If the
format match is low, your IDs are probably from a different assembly.

Choose a background, the GO categories you care about, and an analysis mode, then run it.

## The options, and what each one does

Background (universe). The background is the set of genes the test compares against. By
default it is every annotated gene in the genome. If you have a more relevant background, for
instance only the genes expressed in your tissue or experiment, upload it instead. This is the
recommended choice when you have one, because it makes "over-represented" mean over-represented
relative to what was actually testable, rather than relative to the whole genome.

GO categories. Biological process, molecular function, and cellular component. You can run any
combination. Each is tested and displayed separately.

Analysis mode. Over-representation, split by direction, or GSEA, as described above. The latter
two appear only when you have selected an expression column.

FDR cutoff. The adjusted p-value below which a term is reported as significant. The default is
0.05, the usual 5% false discovery rate.

FDR method. How p-values are corrected for testing many terms at once. Benjamini-Hochberg is
the sensible default. Bonferroni and Holm are stricter, controlling the chance of any false
positive at all rather than the proportion of them. "None" applies no correction, which is
rarely what you want.

Number of pathways to show. How many of the top terms appear in the plots. The table always
holds the full set of significant terms regardless of this setting.

Minimum and maximum pathway size. The smallest and largest GO terms to test, by gene count.
Very small terms are noisy, and very large terms ("metabolic process" and the like) are too
general to be informative, so both ends are trimmed. Defaults are 10 and 500.

Remove redundant terms. GO terms overlap by design, so a result can be cluttered with parent
and child terms describing nearly the same gene set. Switching this on keeps the most
significant term in each overlapping group and drops the rest. The overlap cut-off controls how
aggressive that is: a lower value removes more.

## What you get out

Results table. Every significant term, with its GO ID, name, ontology, gene count, gene and
background ratios, raw and adjusted p-values, and the genes driving the enrichment. It is
sortable, filterable, and downloadable as a CSV.

Dot plot. The top terms per ontology, with dot size showing how many of your genes fall in each
term and colour showing significance. This is the figure you most often see in papers.

Bar plot. The top terms ranked by significance, as -log10(adjusted p-value), with a dashed line
marking your FDR cut-off.

GSEA plots. When you run GSEA, you also get a dot plot of normalised enrichment scores and a
running-score plot for any term you select, showing where its genes sit along the ranked list.

Both the dot plot and the bar plot can be downloaded as PNG or PDF.

## A note on reading the results

GO enrichment describes the gene list you gave it against a background you chose, so both of
those choices shape the answer. The most significant term is not always the most interesting
one, and broad terms can dominate without saying much. It is worth scanning down the list, not
just reading the top hit, and treating the output as a starting point for interpretation rather
than a conclusion.
