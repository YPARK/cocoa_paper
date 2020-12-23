#####################################################
# CoCoA-Diff analysis for Mathys et al. (2019) data #
#####################################################

ROW := data/brain_2018-05-03/features.tsv.gz
COL := data/brain_2018-05-03/barcodes.tsv.gz
MTX := data/brain_2018-05-03/filtered_count_matrix.mtx.gz

PHENO_DATA := data/brain_2018-05-03/phenotypes.csv data/brain_2018-05-03/filtered_column_metadata.txt.gz

CT := $(shell cat data/PsychENCODE.marker | awk '{ ct[$$2]++ } END { for(t in ct) print t }')

DATA := $(foreach c, $(CT), $(foreach ext, mtx cols, result/sorted/$(c).$(ext).gz))

all:

################################################################

step1: result/celltype.annot.gz $(DATA) result/phenotyped.txt.gz

result/celltype.annot.gz: data/PsychENCODE.marker $(MTX) $(ROW) $(COL)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_annotate_col --mtx $(MTX) --col $(COL) --row $(ROW) --ann $< --log_scale --em_iter 1000 --batch_size 10000 --em_tol 1e-4 --verbose --out $(shell echo $@ | sed 's/.annot.gz//')

result/sorted/%.mtx.gz: $(MTX) $(COL) result/sorted/%.select.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_select_col $^ $(shell echo $@ | sed 's/.mtx.gz//g')

result/sorted/%.cols.gz: result/sorted/%.mtx.gz
	touch $@

result/sorted/%.select.gz: result/celltype.annot.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $< | awk '$$2 == "$*"{ print $$1 }' | gzip -c > $@

result/phenotyped.txt.gz: R/find_phenotyped.R $(PHENO_DATA)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ $@

################################################################
# note: "result/phenotyped.txt.gz" contains TAG,IID,{Pheno_types}
NPheno := $(shell [ -f result/phenotyped.txt.gz ] && gzip -cd result/phenotyped.txt.gz | head -n1 | awk '{ print (NF -2) }')

pheno_ := 1

step2: $(foreach c, $(CT), $(foreach f, $(pheno_), result/cocoa/$(c)_$(f).resid_mu.gz))

# % = $(celltype)_$(pheno_col)
result/temp/%.pheno_select.gz: result/phenotyped.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $< | awk -F'\t' -vC=$(shell echo $* | awk -F'_' '{ print $$2 }') 'NR > 1 && $$(C + 2) != -9 && $$(C + 2) != "NA" && length($$(C + 2)) > 0 { print $$1 }' | gzip > $@

# % = $(celltype)_$(pheno_col)
result/temp/%.mtx.gz: result/temp/%.pheno_select.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_select_col result/sorted/$(shell echo $* | awk -F'_' '{ print $$1 }').mtx.gz result/sorted/$(shell echo $* | awk -F'_' '{ print $$1 }').cols.gz $< $(shell echo $@ | sed 's/.mtx.gz//g')

result/temp/%.cols.gz: result/temp/%.mtx.gz
	touch $@

# % = $(celltype)_$(pheno_col)
result/temp/%.trt.gz: R/match_pheno.R result/temp/%.cols.gz result/phenotyped.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ $(shell echo $* | awk -F'_' '{ print $$2 }') $@

# % = $(celltype)_$(pheno_col)
result/temp/%.ind.gz: R/match_pheno.R result/temp/%.cols.gz result/phenotyped.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ 0 $@

# % = $(celltype)_$(pheno_col)
result/temp/%.annot.gz: result/temp/%.cols.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $< | awk '{ print $$1 FS "$*" }' | gzip -c > $@

# % = $(celltype)_$(pheno_col)
result/temp/%.lab.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	echo $* | gzip -c > $@

# % = $(celltype)_$(pheno_col)
result/cocoa/%.resid_mu.gz: result/temp/%.mtx.gz result/temp/%.cols.gz result/temp/%.annot.gz result/temp/%.trt.gz result/temp/%.lab.gz result/temp/%.ind.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	OMP_NUM_THREADS=16 mmutil_cfa_col --mtx result/temp/$*.mtx.gz --col result/temp/$*.cols.gz --annot result/temp/$*.annot.gz --trt result/temp/$*.trt.gz --lab result/temp/$*.lab.gz --ind result/temp/$*.ind.gz --verbose --knn 100 --nboot 100 --rank 50 --log_scale --out $(shell echo $@ | sed 's/.resid_mu.gz//g')

################################################################

step3: $(foreach c, $(CT), result/aggregate/$(c).mean.gz)

# % = $(celltype)
result/aggregate/%.mean.gz: result/sorted/%.mtx.gz result/sorted/%.cols.gz result/aggregate/%.ind.gz result/aggregate/%.lab.gz result/aggregate/%.annot.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_aggregate_col --mtx result/sorted/$*.mtx.gz --col result/sorted/$*.cols.gz  --annot result/aggregate/$*.annot.gz --ind result/aggregate/$*.ind.gz --lab result/aggregate/$*.lab.gz --verbose --out $(shell echo $@ | sed 's/.mean.gz//')

result/aggregate/%.annot.gz: result/sorted/%.cols.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $< | awk '{ print $$1 FS "$*" }' | gzip -c > $@

result/aggregate/%.ind.gz: R/match_pheno.R result/sorted/%.cols.gz result/phenotyped.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ 0 $@

result/aggregate/%.lab.gz: result/sorted/%.cols.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	echo $* | gzip -c > $@

#################
# documentation #
#################

RMD := $(wildcard *.rmd)
HTML := $(RMD:.rmd=.html)
PDF := $(RMD:.rmd=.pdf)
DOCX := $(RMD:.rmd=.docx)

doc: $(HTML) $(DOCX) $(PDF)

%.html: %.rmd $(COMMON)
	Rscript -e "library(rmarkdown); render('$<');"

%.pdf: %.rmd $(COMMON)
	Rscript -e "library(rmarkdown); render('$<', 'pdf_document');"

%.docx: %.rmd $(COMMON)
	Rscript -e "library(rmarkdown); render('$<', 'word_document');"
