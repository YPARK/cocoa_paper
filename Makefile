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

step1: result/celltype.annot.gz $(DATA) result/phenotyped.txt.gz $(foreach f, 1, result/combined/aggregate/combined_$(f).mean.gz result/combined/cocoa/combined_$(f).resid_mu.gz result/combined/stat/$(f).stat.gz)

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

#####################################
# What if cell types were combined? #
#####################################

# % = $(pheno_col) in {1, 2}
result/combined/%.pheno_select.gz: result/phenotyped.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $< | awk -F'\t' -vC=$* 'NR > 1 && $$(C + 2) != -9 && $$(C + 2) != "NA" && length($$(C + 2)) > 0 { print $$1 }' | gzip > $@

# % = $(pheno_col) in {1, 2}
result/combined/%.mtx.gz: result/combined/%.pheno_select.gz $(MTX) $(COL)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_select_col $(MTX) $(COL) $< $(shell echo $@ | sed 's/.mtx.gz//g')

# % = $(pheno_col)
result/combined/%.cols.gz: result/combined/%.mtx.gz
	touch $@

# % = $(pheno_col)
result/combined/cocoa/combined_%.trt.gz: R/match_pheno.R result/combined/%.cols.gz result/phenotyped.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ $* $@

# % = $(pheno_col)
result/combined/cocoa/combined_%.ind.gz: R/match_pheno.R result/combined/%.cols.gz result/phenotyped.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ 0 $@

# % = $(pheno_col)
result/combined/cocoa/combined_%.annot.gz: result/combined/%.cols.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $< | awk '{ print $$1 FS "combined_$*" }' | gzip -c > $@

# % = $(pheno_col)
result/combined/cocoa/combined_%.lab.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	echo "combined" | gzip -c > $@

# % = $(pheno_col)
result/combined/cocoa/combined_%.resid_mu.gz: result/combined/%.mtx.gz result/combined/%.cols.gz result/combined/%.annot.gz result/combined/%.trt.gz result/combined/%.lab.gz result/combined/%.ind.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	OMP_NUM_THREADS=16 mmutil_cfa_col --mtx result/combined/$*.mtx.gz --col result/combined/$*.cols.gz --annot result/combined/$*.annot.gz --trt result/combined/$*.trt.gz --lab result/combined/$*.lab.gz --ind result/combined/$*.ind.gz --verbose --knn 100 --nboot 100 --rank 50 --log_scale --out $(shell echo $@ | sed 's/.resid_mu.gz//g')

# % = $(pheno_col)
result/combined/aggregate/combined_%.mean.gz: result/combined/%.mtx.gz result/combined/%.cols.gz result/combined/%.annot.gz result/combined/%.trt.gz result/combined/%.lab.gz result/combined/%.ind.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_aggregate_col --mtx result/combined/$*.mtx.gz --col result/combined/$*.cols.gz --annot result/combined/$*.annot.gz --lab result/combined/$*.lab.gz --ind result/combined/$*.ind.gz --verbose --out $(shell echo $@ | sed 's/.mean.gz//g')

result/combined/aggregate/combined_%.sum.gz: result/combined/aggregate/combined_%.mean.gz
	@touch $@

result/combined/aggregate/combined_%.mu_cols.gz: result/combined/aggregate/combined_%.mean.gz
	@touch $@

result/combined/stat/%.stat.gz: R/calc_glob_stat.R result/combined/cocoa/combined_%.boot_ln_mu.gz result/combined/cocoa/combined_%.mu_cols.gz result/combined/aggregate/combined_%.sum.gz result/combined/aggregate/combined_%.mean.gz result/combined/aggregate/combined_%.mu_cols.gz data/brain_2018-05-03/features.tsv.gz result/phenotyped.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ $@

################################################################
# note: "result/phenotyped.txt.gz" contains TAG,IID,{Pheno_types}
NPheno := $(shell [ -f result/phenotyped.txt.gz ] && gzip -cd result/phenotyped.txt.gz | head -n1 | awk '{ print (NF -2) }')

pheno_ := 1 2

step2: $(foreach c, $(CT), $(foreach f, 1 2, result/cocoa/$(c)_$(f).resid_mu.gz) $(foreach f, 1 2 3 4, result/aggregate/$(c)_$(f).mean.gz))

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

# % = $(celltype)_$(pheno_col)
result/aggregate/%.mean.gz: result/temp/%.mtx.gz result/temp/%.cols.gz result/temp/%.annot.gz result/temp/%.trt.gz result/temp/%.lab.gz result/temp/%.ind.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_aggregate_col --mtx result/temp/$*.mtx.gz --col result/temp/$*.cols.gz --annot result/temp/$*.annot.gz --lab result/temp/$*.lab.gz --ind result/temp/$*.ind.gz --verbose --out $(shell echo $@ | sed 's/.mean.gz//g')

################################################################
# Statistical test after adjustment
GLOB_STAT := $(foreach c, $(CT), $(foreach f, $(pheno_), result/glob_stat/$(c)_$(f).stat.gz result/glob_stat/$(c)_$(f).cf_stat.gz))

step3: $(GLOB_STAT)

result/glob_stat/%.stat.gz: R/calc_glob_stat.R result/cocoa/%.boot_ln_mu.gz result/cocoa/%.mu_cols.gz result/aggregate/%.sum.gz result/aggregate/%.mean.gz result/aggregate/%.mu_cols.gz data/brain_2018-05-03/features.tsv.gz result/phenotyped.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ $@

result/glob_stat/%.cf_stat.gz: R/calc_glob_cf_stat.R result/cocoa/%.cf_mu.gz result/cocoa/%.mu_cols.gz result/aggregate/%.sum.gz result/aggregate/%.mean.gz result/aggregate/%.mu_cols.gz data/brain_2018-05-03/features.tsv.gz result/phenotyped.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ $@

#################
# documentation #
#################

RMD := $(wildcard *.rmd)
HTML := $(RMD:.rmd=.html)
PDF := $(RMD:.rmd=.pdf)

doc: $(HTML) $(PDF)

%.html: %.rmd $(COMMON)
	Rscript -e "library(rmarkdown); render('$<');"

%.pdf: %.rmd $(COMMON)
	Rscript -e "library(rmarkdown); render('$<', 'pdf_document');"

%.docx: %.rmd $(COMMON)
	Rscript -e "library(rmarkdown); render('$<', 'word_document');"
