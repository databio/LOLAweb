


- [What is LOLAweb and what is it used for?](#what)
- [How do I use LOLAweb?](#how-to-use)
- [How does LOLAweb evaluate overlap?](#overlap-evaluation)
- [What are the values in the results table?](#results-table)
- [How do I cite LOLAweb?](#how-to-cite)
- [What universe should I choose?](#uni)
- [How could I define my own universe?](#custom-universe)
- [What is in the LOLAweb region set databases?](#database)
- [What are the example query sets?](#examples)
- [Can I use a custom database with LOLAweb?](#customdb)
- [How do I run LOLAweb locally?](#run-locally)
- [Can I share my LOLAweb results with others?](#share-results)
- [How long are LOLAweb results kept?](#results-kept)
- [What can I find more information or support?](#other-resources)

----------------------------------------

<a name="what"></a>
### What is LOLAweb and what is it used for?

LOLAweb is a web server and interactive results viewer for enrichment of overlap
between a user-provided query region set (a bed file) and a database of region
sets. It provides an interactive result explorer to visualize the highest ranked
enrichments from the database. LOLAweb is a web interface to the [LOLA R
package](http://bioconductor.org/packages/LOLA/).

LOLAweb is useful for exploring enrichment of genomic ranges. Frequently, we
uncover sets of genomic regions as relevant for some particular biological
scenario. For example, they may be binding sites of a given transcription
factor, regions that change in methylation after some perturbation, regions that
specify differences between two species or cell states, and so on. From the way
we've found the regions, we know something about them, but we'd now like to know
if anything else is known about these particular regions. Have other experiments
uncovered them as important for something? Have they been highlighted in a
large-scale project in some other related cell type?

This is where LOLAweb can help. We test the query regions for overlap with
thousands of existing database region sets. By ranking the overlap, we can see
which existing datasets are most similar to the new query set. This enables us
to tie previous knowledge to our new genomic regions.


--------------------------------------------------------------------------------

<a name="how-to-use"> </a>
### How do I use LOLAweb?
<img alt="LOLAweb workflow" title="LOLAweb workflow" src="LOLAweb-workflow.png" style="padding:30px; float:right; width:50%; max-width:600px">
Using LOLAweb is as easy as 1-2-3. First, all you really need is a `BED`-like
file (with at least 3 columns: chromosome, start, and end) identifying regions
of interest. These could be the result of a differential ChIP-seq or ATAC-seq
experiment, or anything, really. Upload that file by dragging it into the first
query box. You'll need to indicate which reference assembly your file subscribes
to.

Second, choose a universe. If you have a set of background regions, upload them
here; otherwise, choose one of our built-in universe sets. You can read more
about these choices later on this page.

Finally, choose a database. We offer a few options; just start with the default
Core database at first, but feel free to look around. If you want, you can
always use a custom database by running either LOLAweb or LOLA locally.

Hit "Run LOLA" and watch the magic happen! The results will be displayed as an
interactive interface letting you sort and filter the top hits. The schematic
to the right demonstrates the whole process.


--------------------------------------------------------------------------------

<a name="overlap-evaluation"> </a>
### How does LOLAweb evaluate overlap?

LOLA evaluates overlaps by comparing the query region set to each database region set and calculates the number of overlapping regions for each pairwise comparison. Along with a similar calculation for the universe region set, LOLA uses the number of overlaps and non-overlaps to build a contingency table, and then uses a Fisher's exact test to assess the significance of the overlap. After computing these statistics for each comparison, LOLA ranks each database region set and provides a ranked summary of the top database sets. This procedure effectively pulls out the region sets in the database that are most similar to the query region set. 

LOLA uses three summary statistics to assess the degree of overlap: 1) the p-value and 2) odds ratio from a Fisher's exact test, and 3) the raw number of overlapping regions (called *support*). Each of these statistics emphasizes a different aspect of the comparison; for example, the number of overlapping regions emphasizes sheer volume of overlapping regions but does not correct for significance, while the odds ratio emphasizes relative enrichment but can be dominated by sets with small numbers of regions. To come up with an aggregate score that shares strengths of each of these statistics, LOLA ranks each pairwise comparison for each of these three statistics independently, and then calculates a combined rank for each region set by assigning it the worst (max) rank among these three. To rank highly in the combined rank, then, requires that a comparison do reasonably well on all three measures, because the worst score is taken. In our experience, this process prioritizes biologically relevant associations and eliminates spurious relationships.

--------------------------------------------------------------------------------
<a name="results-table"> </a>
### What are the values in the results table?

The columns in the results table are: 

- **userSet and dbSet**: index into their respective input region sets.
- **pvalueLog**: -log10(pvalue) from the fisher's exact result
- **oddsRatio**: result from the fisher's exact test
- **support**: number of regions in userSet overlapping databaseSet
- **rnkPV, rnkOR, rnkSup:** rank in this table of p-value, oddsRatio, and Support respectively.
- **maxRnk, meanRnk**: max and mean of the 3 previous ranks, providing a combined ranking system. 
- **b, c, d**: 3 other values completing the 2x2 contingency table (with support). 
- The remaining columns describe the dbSet for the row.

--------------------------------------------------------------------------------

<a name="how-to-cite"> </a>
### How do I cite LOLAweb?

LOLAweb is pending publication. In the interim, please cite the LOLA R package:

N. C. Sheffield and C. Bock, “LOLA: enrichment analysis for genomic region sets and regulatory elements in R and Bioconductor,” Bioinformatics, vol. 32, no. 4, pp. 587–589, Oct. 2016.

[http://dx.doi.org/10.1093/bioinformatics/btv612](http://dx.doi.org/10.1093/bioinformatics/btv612)

--------------------------------------------------------------------------------

<a name="uni"> </a>
### What universe should I choose?

One of the key questions when you run LOLA is what your background set, or
"universe" is. The universe is a bit open-ended, and it's reasonable to try a
few different things. Changing the universe isn't right or wrong, it just
changes the question you are asking. The universe set is tested for overlap with
the database, and the counts are used in the contingency tables that determine
significance for your user sets. The reason this is important is that if you
have some regions which were never really possible to end up in a region set of
interest, it's unfair to penalize your regions for not overlapping those regions
in the database, changing the results of the significance test.

LOLAweb provides 3 basic options. These options are briefly introduced here, and
then details on how to decide are provided below:

**Use pre-loaded universe**: Here you can select from a series of universes we
have already created for you. These universes vary based on the reference genome
assembly you have chosen. They include tiling regions, as well as a manually
curated universe set that was derived by merging all of the DNAse
hypersensitivity data from over 100 cell types from the ENCODE project.

**Build universe from user sets**: This option tells LOLAweb to use your query user sets to build a new universe. To use this, you'll have to have uploaded at least 2 files for query. This will simply combine those files into a background set. This will let you test for differential enrichment within
your usersets.

**Upload universe**: Finally, you're free to create and upload any universe you want to create. Details on what you might want to use this for follow.

<img src="universe_selection.svg" style="padding:30px; width:50%; float:right; max-width:800px">

Let's imagine you've done an experiment where you're testing how some epigenomic
signal (say, H3K27ac or DNA methylation) responds to a perturbation. You end up
with a set of genomic regions that were covered by your assay, along with two
subsets of interest: those that increased, and those that decreased in signal
(see figure).

In this case, one approach is to think of this as two user query sets: one for
each subset of regions that change. Now, what should your universe be? Well, you
have a few options. Let's go through 4 common universe choices see what they
mean about the significance test you're performing. The four choices move in
increasing order of specificity:

<b>Choice 1: genome tiles</b>. LOLAweb lets you choose tiling regions of various
sizes. These tiles are little better than generic splits of the genome because
they at least account for reference assembly gaps. But otherwise, this is a
pretty basic background set which doesn't restrict the question at all. This
would be a useful universe to choose if you are not sure what to expect, and
your experiment could have included any genomic region as output. The result
will tell you what database region sets are enriched in your query sets
generally, across the whole genome.

<b>Choice 2: any regulatory elements</b>. A bit more restricted universe set is
the set of active regulatory elements. LOLAweb includes a manually curated
universe set that was derived by merging all of the DNAse hypersensitivity data
from over 100 cell types from the ENCODE project. This set will restrict your
significance test to only active regulatory elements, so it's a reasonable
background set if you're looking at ChIP-seq data or some other active signal,
but may not be relevant for every experiment type. In practice, we find this to
be a useful first pass test for many experiments. Relative to the tiling
approach, this, for example, won't penalize your overlap test for not
overlapping a database region that's hidden away in a heterochromatic region. It
yields which database region sets are enriched in your query set, but only
considering active parts of the genome. An advantage of this approach is that
it's easy: we've already defined that set for you and you can just select it. A
disadvantage is that you're restricted to elements that were covered in the
considered cell types (although, offsetting that is the point that really
specialized stuff is less likely to be in the database, anyway).

<b>Choice 3: all regions covered</b>. Probably a better choice (though requiring
a bit more work) is to derive the universe from your actual data. The natural
choice would be to use the set of regions your experiment actually covered. You
can think of this universe as the set of regions you tested for possible
inclusion in your user sets; or, in other words, it is the restricted background
set of regions that were tested, including everything in your regions of
interest as well as those that were tested but did not get included. If you can
easily come up with such a test for your problem, this is a really good choice
for background universe. The resulting test will give you overlaps with database
regions without penalizing for areas that weren't even covered in your
experiment anyway. This asks: what's special about my regions that increased,
versus all the regions that didn't?

<b>Choice 4: changing regions only</b> Finally, you could make a very
restrictive universe by just combining your two region sets of interest. In this
case, it would be the set of regions that either increased or decreased. The
result of this test will give you the enrichment of the increased regions
<i>relative to changing regions</i>. This is a subtly different question to ask
than using the set of all covered regions (choice 3). The LOLA R package
provides an easy way to do this with the <i>buildRestrictedUniverse</i>
function, and you can accomplish the same thing with the radio button on
LOLAweb. It assumes that we're only interested in regions that change, then we
ask the question: what's special about the ones that increased?

--------------------------------------------------------------------------------

<a name="custom-universe"></a>
### How could I define my own universe?

If you want complete control over the question you're asking, the best thing to
do is define your own universe. How to do this depends on the particular data
you're looking at. Here we outline a few guidelines for common data types to get
you started in how to think about the universe:

**DNA methylation data**: You universe could be all regions that had reasonable
coverage of methylation reads, and your user sets could be those that were either
highly methylated or lowly methylated (or differentially methylated). The universe is the set of CpGs that had
enough reads that they could have been differentially methylated, even if they
weren't. The universe could be quite different for RRBS (reduced representation
bisulfite sequencing) vs. WGBS (whole-genome bisulfite sequencing) protocols.

**Chip-seq data**: The universe could be the set of 1kb tiles that were covered
by an Input experiment. This test would show where ChIP peaks are generally
enriched relative to the whole genome, and will likely show enrichment in active
areas generally. Alternatively, the universe could be the union of all
transcription factor binding sites from large-scale projects like ENCODE, or the
pre- loaded set of active DNase hypersensitive sites. This universe would test
the enrichment of your ChIP peaks relative to other ChIP peaks, which will be
more likely to highlight cell-type differences. Finally, perhaps you have
differential peaks between two conditions. In this case, you could build a
*restricted* universe of just peaks covered by these two experiments. This test
would then show you how your increasing peaks are enriched relative only to
peaks in this one cell type.

--------------------------------------------------------------------------------

<a name="database"></a>

### What is in the LOLAweb region set databases?

LOLAweb provides access to a few curated databases we have put together for you.
You can also [run LOLAweb on a custom database](#customdb), if you want. But the
online version provides a few useful databases that will suffice for the
majority of users.

When you select a genome assembly, LOLAweb will automatically populate the
database dropdown menu with all available databases for the chosen assembly.
Here are the contents of the databases currently available to LOLAweb users:

**LOLA Core**:

* hg19/hg38
  1. Transcription Factor binding sites from [ENCODE](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/)
  2. Tissue clustered DNase hypersensitive sites from [Sheffield et al. (2013)](http://dnase.genome.duke.edu)
  3. [Codex database](http://codex.stemcells.cam.ac.uk/)
  4. A collection of UCSC feature tables (like CpG islands)
  5. Cistrome database from [Cistrome](http://dx.doi.org/10.1186/gb-2011-12-8-r83)
  6. Epigenome databases from [Cistrome](http://dx.doi.org/10.1186/gb-2011-12-8-r83)
* mm10/mm9
  1. [Codex database](http://codex.stemcells.cam.ac.uk/)
  2. Cistrome database
  3. Transcription Factor binding sites from ENCODE
**LOLA Roadmap**:

* hg19/hg38
  All files downloaded from the Roadmap Epigenomics project.

**LOLA Jaspar**:
* hg19/hg38
  We scanned the genome assembly for motif matches to every motif in the JASPAR
  database. Each motif has an entry in the database with all its matches. This
  enables you to do a type of simple motif enrichment analysis to see if your
  regions significantly overlap the set of motif matches for any known motifs.

You can also download these databases and explore them yourself at
[big.databio.org/regiondb](http://big.databio.org/regiondb/).

--------------------------------------------------------------------------------

<a name="examples"></a>
### What are the example user query sets?

**lamina.bed**. This is a set of regions that are identified as being attached
to the nuclear periphery.

**ewing_DHS.bed**. This is a very small example of some DNase hypersensitivity
sites that are specific to Ewing sarcoma, a rare pediatric cancer. This is the
example set used in the LOLA vignettes and the data were from Sheffield et al.
(2017).

**mesenchymal_DHS.bed**. This example set is derived from ENCODE DNase
hypersensitivity data. It was derived by taking DNase hypersensitive sites with
high scores in mesenchymal cell types and low scores in all other cell types.


--------------------------------------------------------------------------------

<a name="customdb"></a>
### Can I use a custom database with LOLAweb?

Yes, you can, but you'll have to run LOLAweb locally because it's not practical
to let users upload an entire database. Luckily, we've containerized LOLAweb, so
it's relatively easy to run it locally. Alternatively, the [LOLA R
package](http://bioconductor.org/packages/LOLA/) can work with custom databases
and may be simpler, if you don't want to deal with a containerized setup or
running your own local shiny server. But if you do want to run LOLAweb locally,
check out *[how to run LOLAweb locally](#run-locally)*.

--------------------------------------------------------------------------------

<a name="run-locally"></a>
### How do I run LOLAweb locally?

You'll need to be familiar with [Docker](http://docker.com). You can find the Docker container for LOLAweb at [Dockerhub (https://hub.docker.com/r/databio/lolaweb/)](https://hub.docker.com/r/databio/lolaweb/). Instructions for run a local LOLAweb instance are here: https://github.com/databio/LOLAweb/tree/master/docker


--------------------------------------------------------------------------------

<a name="share-results"></a>
### Can I share my LOLAweb results with others?
Yes! Once you run your region sets through LOLAweb, you are redirected to a web
page unique to your data. It will look something like this:

<div style="margin-left:20px;padding:10px;"><code>http://lolaweb.databio.org/?key=VWQN3ZC5HFD92EK</code></div>

This link can be sent to other users to share your findings. We do not store
your original data or any identifying information about who submitted it.

--------------------------------------------------------------------------------

<a name="results-kept"></a>
### How long are LOLAweb results kept?

LOLAweb results are stored on our servers for one year, then automatically deleted.

<a name="other-resources"></a>
### Where can I find more information or support?

Here are some links to other resources related to LOLA and LOLAweb:

- [LOLAweb issue tracker](https://github.com/databio/LOLAweb/issues) - please use this for support requests for LOLAweb.
- [LOLAweb GitHub source code](https://github.com/databio/LOLAweb)
- [Dockerhub docker image](https://hub.docker.com/r/databio/lolaweb/)
- Manuscript (pending)
- [LOLA R package at Bioconductor](http://bioconductor.org/packages/LOLA/)
- [LOLA R package source code at GitHub](https://github.com/nsheff/LOLA)
- [LOLA R package documentation](http://code.databio.org/LOLA)
- [Published Bioinformatics paper describing LOLA](http://dx.doi.org/10.1093/bioinformatics/btv612)
- [Region set databases and vignette data](http://cloud.databio.org)




