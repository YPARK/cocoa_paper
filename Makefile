#####################################################
# CoCoA-Diff analysis for Mathys et al. (2019) data #
#####################################################

ROW := data/brain_2018-05-03/features.tsv.gz
COL := data/brain_2018-05-03/barcodes.tsv.gz
MTX := data/brain_2018-05-03/filtered_count_matrix.mtx.gz

PHENO_DATA := data/brain_2018-05-03/phenotypes.csv data/brain_2018-05-03/filtered_column_metadata.txt.gz

CT := Ex In Astro Microglia OPC Oligo

all:

################################################################

BBKNN := $(foreach t, mtx factors svd_D svd_V svd_U cols, result/bbknn/Total.$(t).gz)
DATA := $(foreach c, $(CT), $(foreach ext, mtx cols, result/sorted/$(c).$(ext).gz))

step1: $(BBKNN) result/celltype.annot.gz result/bbknn_celltype.annot.gz $(DATA) result/phenotyped.txt.gz 

result/celltype.annot.gz: data/PsychENCODE.marker $(MTX) $(COL) $(ROW)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_annotate_col --mtx $(MTX) --col $(COL) --row $(ROW) --ann $< --em_iter 1000 --batch_size 10000 --em_tol 1e-6 --verbose --out $(shell echo $@ | sed 's/.annot.gz//')

result/bbknn_celltype.annot.gz: data/PsychENCODE.marker $(BBKNN) $(ROW)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_annotate_col --svd_u result/bbknn/Total.svd_U.gz --svd_v result/bbknn/Total.factors.gz --svd_d result/bbknn/Total.svd_D.gz --col $(COL) --row $(ROW) --ann $< --em_iter 1000 --batch_size 10000 --em_tol 1e-6 --verbose --out $(shell echo $@ | sed 's/.annot.gz//')

result/sorted/%.mtx.gz: $(MTX) $(COL) result/sorted/%.select.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_select_col $^ $(shell echo $@ | sed 's/.mtx.gz//g')

result/sorted/%.cols.gz: result/sorted/%.mtx.gz
	touch $@

result/sorted/%.select.gz: result/celltype.annot.gz result/bbknn_celltype.annot.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $^ | awk '$$2 == "$*" { n[$$1]++ } END { for(j in n) if(n[j] > 1) print j }' | gzip -c > $@

result/bbknn/Total.mtx.gz: $(MTX) $(COL) result/bbknn/Total.batch.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_bbknn --mtx $(MTX) --col $(COL) --batch result/bbknn/Total.batch.gz --knn 100 --rank 50 --out result/bbknn/Total --verbose

result/bbknn/Total.cols.gz: result/bbknn/Total.mtx.gz
	touch $@

result/bbknn/Total.svd_D.gz: result/bbknn/Total.mtx.gz
	touch $@

result/bbknn/Total.svd_U.gz: result/bbknn/Total.mtx.gz
	touch $@

result/bbknn/Total.factors.gz: result/bbknn/Total.mtx.gz
	touch $@

result/bbknn/Total.batch.gz: $(COL)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $< | awk -F'.' '{ print $$2 }' | gzip -c > $@

################################################################

CTS := Ex-L2or3 Ex-L4 Ex-L5or6 In-PV In-SST In-SV2C In-VIP
CTSS := $(foreach s, $(CTS), subtype-$(s))

BBKNN_SUBTYPE := $(foreach ct, Ex In, $(foreach t, mtx factors svd_D svd_V svd_U cols, result/bbknn/$(ct).$(t).gz))
ANNOT_SUBTYPE := $(foreach ct, Ex In, result/subtype/bbknn_$(ct).annot.gz)
DATA_SUBTYPE := $(foreach s, $(CTS), result/sorted/subtype-$(s).mtx.gz)

step1a: $(BBKNN_SUBTYPE) $(ANNOT_SUBTYPE) $(DATA_SUBTYPE)

result/sorted/subtype-Ex-%.select.gz: result/subtype/bbknn_Ex.annot.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $^ | awk '$$2 == "Ex-$*" { print $$1 }' | gzip -c > $@

result/sorted/subtype-In-%.select.gz: result/subtype/bbknn_In.annot.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $^ | awk '$$2 == "In-$*" { print $$1 }' | gzip -c > $@

result/sorted/subtype-%.mtx.gz: $(MTX) $(COL) result/sorted/subtype-%.select.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_select_col $^ $(shell echo $@ | sed 's/.mtx.gz//g')

result/bbknn/%.mtx.gz: result/sorted/%.mtx.gz result/sorted/%.cols.gz result/bbknn/%.batch.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_bbknn --mtx result/sorted/$*.mtx.gz  --col result/sorted/$*.cols.gz --batch result/bbknn/$*.batch.gz --knn 100 --rank 50 --out result/bbknn/$* --verbose

result/bbknn/%.factors.gz: result/bbknn/%.mtx.gz
	touch $@
result/bbknn/%.cols.gz: result/bbknn/%.mtx.gz
	touch $@
result/bbknn/%.svd_D.gz: result/bbknn/%.mtx.gz
	touch $@
result/bbknn/%.svd_V.gz: result/bbknn/%.mtx.gz
	touch $@
result/bbknn/%.svd_U.gz: result/bbknn/%.mtx.gz
	touch $@

result/bbknn/%.batch.gz: result/sorted/%.cols.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $< | awk -F'.' '{ print $$2 }' | gzip -c > $@

result/subtype/In.marker.gz: data/Velmeshev_etal_science_2019_markers.csv
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $^ | awk -F',' '/^IN-/{ gsub("IN","In",$$1);  print $$3 " " $$1 }' | gzip -c > $@

result/subtype/Ex.marker.gz: data/Velmeshev_etal_science_2019_markers.csv
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $^ | awk -F',' '/^L[1-9]+/ { gsub("/","or",$$1); print $$3 " " ("Ex-" $$1) }' | gzip -c > $@

result/subtype/bbknn_%.annot.gz: result/subtype/%.marker.gz $(BBKNN_SUBTYPE)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_annotate_col --svd_u result/bbknn/$*.svd_U.gz --svd_d result/bbknn/$*.svd_D.gz --svd_v result/bbknn/$*.factors.gz --col result/bbknn/$*.cols.gz --row $(ROW) --ann $< --em_iter 1000 --batch_size 5000 --em_tol 1e-8 --verbose --out $(shell echo $@ | sed 's/.annot.gz//')

################################################################
# note: "result/phenotyped.txt.gz" contains TAG,IID,{Pheno_types}

NPheno := $(shell [ -f result/phenotyped.txt.gz ] && gzip -cd result/phenotyped.txt.gz | head -n1 | awk '{ print (NF -2) }')
pheno_ := 1 2 3 4 5
EXT_ := resid_mu

step2: $(foreach f, $(pheno_), $(foreach t, $(EXT_), result/cocoa/$(f).$(t).gz)) $(foreach f, $(pheno_), result/aggregate/$(f).mean.gz)

result/phenotyped.txt.gz: R/find_phenotyped.R $(PHENO_DATA)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ $@

GLIAL := Microglia Astro OPC Oligo
Neun := Ex In

# % = $(pheno_col)
result/temp/%.pheno_select.gz: result/phenotyped.txt.gz result/temp/total_annot.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $< | awk -F'\t' -vA=result/temp/total_annot.gz -vC=$* 'BEGIN { while(("cat " A " | gzip -cd " | getline line) > 0) { split(line,larr," "); annot[larr[1]] = larr[2]; } } NR > 1 && $$(C + 2) != -9 && $$(C + 2) != "NA" && length($$(C + 2)) > 0 && ($$1 in annot) { print $$1 }' | gzip > $@

# % = $(pheno_col)
result/temp/%.mtx.gz: $(MTX) $(COL) result/temp/%.pheno_select.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_select_col $^ $(shell echo $@ | sed 's/.mtx.gz//g')

result/temp/%.cols.gz: result/temp/%.mtx.gz
	touch $@

# % = $(pheno_col)
result/temp/%.trt.gz: R/match_pheno.R result/temp/%.cols.gz result/phenotyped.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ $(shell echo $* | awk -F'_' '{ print $$NF }') $@

# % = $(pheno_col)
result/temp/%.ind.gz: R/match_pheno.R result/temp/%.cols.gz result/phenotyped.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ 0 $@

# % = $(pheno_col)
result/temp/%.annot.gz: result/temp/total_annot.gz result/temp/%.cols.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd result/temp/$*.cols.gz | awk -vA=$< 'BEGIN { while(("cat " A " | gzip -cd " | getline line) > 0) { split(line,larr," "); annot[larr[1]] = larr[2]; } } { a = ($$1 in annot) ? annot[$$1] : "NA"; print $$1 FS a }' | gzip -c > $@

# % = $(pheno_col)
result/temp/%.lab.gz: result/temp/%.annot.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $< | awk '{ n[$$2]++ } END {for(t in n) print t }' | sort | gzip -c > $@

result/temp/total_annot.gz: $(foreach g, $(GLIAL), result/temp/$(g).glial.gz)  $(foreach n, $(Neun), result/temp/$(n).subtype.gz)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $^ > $@

result/temp/%.glial.gz: result/celltype.annot.gz result/bbknn_celltype.annot.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $^ | awk '$$2 == "$*" && $$3 > .99 { n[$$1]++; m[$$1]=$$2 } END { for(j in n) if(n[j] > 1) print j FS m[j] }' | gzip -c > $@

result/temp/%.subtype.gz: result/subtype/bbknn_%.annot.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $^ | awk '$$3 > .99 { print $$1 FS $$2 }' | gzip -c > $@

# % = $(pheno_col)
result/cocoa/%.resid_mu.gz: result/temp/%.mtx.gz result/temp/%.cols.gz result/temp/%.annot.gz result/temp/%.trt.gz result/temp/%.lab.gz result/temp/%.ind.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	OMP_NUM_THREADS=16 mmutil_cfa_col --mtx result/temp/$*.mtx.gz --col result/temp/$*.cols.gz --annot result/temp/$*.annot.gz --trt result/temp/$*.trt.gz --lab result/temp/$*.lab.gz --ind result/temp/$*.ind.gz --verbose --knn 100 --rank 50 --log_scale --out $(shell echo $@ | sed 's/.resid_mu.gz//g')

# % = $(pheno_col)
result/cocoa_knn/%.resid_mu.gz: result/temp/%.mtx.gz result/temp/%.cols.gz result/temp/%.annot.gz result/temp/%.trt.gz result/temp/%.lab.gz result/temp/%.ind.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	OMP_NUM_THREADS=16 mmutil_cfa_col --mtx result/temp/$*.mtx.gz --col result/temp/$*.cols.gz --annot result/temp/$*.annot.gz --trt result/temp/$*.trt.gz --lab result/temp/$*.lab.gz --ind result/temp/$*.ind.gz --verbose --knn 100 --rank 50 --impute_knn --log_scale --out $(shell echo $@ | sed 's/.resid_mu.gz//g')

# % = $(pheno_col)
result/aggregate/%.mean.gz: result/temp/%.mtx.gz result/temp/%.cols.gz result/temp/%.annot.gz result/temp/%.trt.gz result/temp/%.lab.gz result/temp/%.ind.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_aggregate_col --mtx result/temp/$*.mtx.gz --col result/temp/$*.cols.gz --annot result/temp/$*.annot.gz --lab result/temp/$*.lab.gz --ind result/temp/$*.ind.gz --verbose --out $(shell echo $@ | sed 's/.mean.gz//g')

################################################################
# Statistical test after adjustment
GLOB_STAT := $(foreach f, $(pheno_), $(foreach s, stat, result/glob_stat/$(f).$(s).gz))

step3: $(GLOB_STAT)

result/glob_stat/%.stat.gz: R/calc_glob_stat.R result/cocoa/%.ln_resid_mu.gz result/cocoa/%.ln_resid_mu_sd.gz result/cocoa/%.cf_mu.gz result/cocoa/%.cf_mu_sd.gz result/cocoa/%.mu_cols.gz result/aggregate/%.sum.gz result/aggregate/%.mean.gz result/aggregate/%.mu_cols.gz data/brain_2018-05-03/features.tsv.gz result/phenotyped.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ $* $@

#################
# documentation #
#################

RMD := $(wildcard *.rmd)
HTML := $(RMD:.rmd=.html)
PDF := $(RMD:.rmd=.pdf)

doc: $(HTML)

%.html: %.rmd $(COMMON)
	Rscript -e "library(rmarkdown); render('$<');"

%.pdf: %.rmd $(COMMON)
	Rscript -e "library(rmarkdown); render('$<', 'pdf_document');"

%.docx: %.rmd $(COMMON)
	Rscript -e "library(rmarkdown); render('$<', 'word_document');"
