IND := 40
NGENE := 10000
NCAUSAL := 100
NCOVAR_S := 5
NCOVAR_B := 0
NCELL_IND := 50
rhoA := 6.25 # Estimated from Microglia
rhoB := 6.25 # mean = 1, var = 0.16

RSEED := $(shell seq 1 10)
PF := 0

all: $(foreach r, $(RSEED), $(foreach pa, 5, $(foreach p0, 5 0, $(foreach p1, 3 2 1 0, summary/$(p1)_$(p0)_$(pa)_$(PF)_$(r).eval.gz))))

##############
# evaluation #
##############

KNN := 1 10 50 100 200
EXT := ln_resid_mu mu_cols

# % = $(pve_1)_$(pve_0)_$(pve_a)_$(pf)_$(rseed)
summary/%.eval.gz: $(foreach knn, $(KNN), $(foreach ext, $(EXT), $(knn)/%.$(ext).gz)) data/%.causal.gz data/%.label.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	echo $@
	@Rscript --vanilla ../R/evalute_v5.R $* $@
	@rm -f $(foreach k, $(KNN), $(k)/$**.*.gz)
	@rm -f data/$**.*.gz
	@rm -f data/$**.*.gz.index

####################
# model estimation #
####################

200/%.ln_resid_mu.gz: data/%.cols.gz data/%.mtx.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@gzip -cd $< | awk '{ print $$1 FS "$*" }' | gzip -c > 200/$*.annot.gz
	@echo $* | gzip > 200/$*.lab.gz
	@gzip -cd data/$*.ind.gz | awk 'BEGIN { j=0; while(("gzip -cd data/$*.label.gz " | getline line)) { v=sprintf("Ind%02d",(++j)); d[v]=line } } { print d[$$1] }' | gzip -c > 200/$*.trt.gz
	@mmutil_cfa_col --mtx data/$*.mtx.gz --col data/$*.cols.gz --annot 200/$*.annot.gz --trt 200/$*.trt.gz --lab 200/$*.lab.gz --ind data/$*.ind.gz --verbose --knn 200 --rank 50 --nboot 0 --gamma_a0 1 --gamma_b0 1 --log_scale --out $(shell echo $@ | sed 's/.ln_resid_mu.gz//g')

100/%.ln_resid_mu.gz: data/%.cols.gz data/%.mtx.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@gzip -cd $< | awk '{ print $$1 FS "$*" }' | gzip -c > 100/$*.annot.gz
	@echo $* | gzip > 100/$*.lab.gz
	@gzip -cd data/$*.ind.gz | awk 'BEGIN { j=0; while(("gzip -cd data/$*.label.gz " | getline line)) { v=sprintf("Ind%02d",(++j)); d[v]=line } } { print d[$$1] }' | gzip -c > 100/$*.trt.gz
	@mmutil_cfa_col --mtx data/$*.mtx.gz --col data/$*.cols.gz --annot 100/$*.annot.gz --trt 100/$*.trt.gz --lab 100/$*.lab.gz --ind data/$*.ind.gz --verbose --knn 100 --rank 50 --nboot 0 --gamma_a0 1 --gamma_b0 1 --log_scale --out $(shell echo $@ | sed 's/.ln_resid_mu.gz//g')

50/%.ln_resid_mu.gz: data/%.cols.gz data/%.mtx.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@gzip -cd $< | awk '{ print $$1 FS "$*" }' | gzip -c > 50/$*.annot.gz
	@echo $* | gzip > 50/$*.lab.gz
	@gzip -cd data/$*.ind.gz | awk 'BEGIN { j=0; while(("gzip -cd data/$*.label.gz " | getline line)) { v=sprintf("Ind%02d",(++j)); d[v]=line } } { print d[$$1] }' | gzip -c > 50/$*.trt.gz
	@mmutil_cfa_col --mtx data/$*.mtx.gz --col data/$*.cols.gz --annot 50/$*.annot.gz --trt 50/$*.trt.gz --lab 50/$*.lab.gz --ind data/$*.ind.gz --verbose --knn 50 --rank 50 --nboot 0 --gamma_a0 1 --gamma_b0 1 --log_scale --out $(shell echo $@ | sed 's/.ln_resid_mu.gz//g')

10/%.ln_resid_mu.gz: data/%.cols.gz data/%.mtx.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@gzip -cd $< | awk '{ print $$1 FS "$*" }' | gzip -c > 10/$*.annot.gz
	@echo $* | gzip > 10/$*.lab.gz
	@gzip -cd data/$*.ind.gz | awk 'BEGIN { j=0; while(("gzip -cd data/$*.label.gz " | getline line)) { v=sprintf("Ind%02d",(++j)); d[v]=line } } { print d[$$1] }' | gzip -c > 10/$*.trt.gz
	@mmutil_cfa_col --mtx data/$*.mtx.gz --col data/$*.cols.gz --annot 10/$*.annot.gz --trt 10/$*.trt.gz --lab 10/$*.lab.gz --ind data/$*.ind.gz --verbose --knn 10 --rank 50 --nboot 0 --gamma_a0 1 --gamma_b0 1 --log_scale --out $(shell echo $@ | sed 's/.ln_resid_mu.gz//g')

1/%.ln_resid_mu.gz: data/%.cols.gz data/%.mtx.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@gzip -cd $< | awk '{ print $$1 FS "$*" }' | gzip -c > 1/$*.annot.gz
	@echo $* | gzip > 1/$*.lab.gz
	@gzip -cd data/$*.ind.gz | awk 'BEGIN { j=0; while(("gzip -cd data/$*.label.gz " | getline line)) { v=sprintf("Ind%02d",(++j)); d[v]=line } } { print d[$$1] }' | gzip -c > 1/$*.trt.gz
	@mmutil_cfa_col --mtx data/$*.mtx.gz --col data/$*.cols.gz --annot 1/$*.annot.gz --trt 1/$*.trt.gz --lab 1/$*.lab.gz --ind data/$*.ind.gz --verbose --knn 1 --rank 50 --nboot 0 --gamma_a0 1 --gamma_b0 1 --log_scale --out $(shell echo $@ | sed 's/.ln_resid_mu.gz//g')

# % = $(knn)/$(data_flag)
%.mu_cols.gz: %.ln_resid_mu.gz
	@touch $@

###############################
# simulate matrix market data #
###############################

data/%.mtx.gz: data/%.lambda.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mmutil_simulate --mu $< --rho data/$*.rho.gz --col data/$*.cols.gz --out data/$*

data/%.lambda.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@Rscript --vanilla ../R/simulate_glm_v5.R $(IND) $(NGENE) $(NCAUSAL) $(NCOVAR_S) $(NCOVAR_B) $(NCELL_IND) $(shell echo $* | awk -F'_' '{ for(j=1; j<=(NF - 1); ++j) print " 0." $$j; print $$NF }') $(rhoA) $(rhoB) data/$*

data/%.cols.gz: data/%.mtx.gz
	@touch $@

data/%.ind.gz: data/%.mtx.gz
	@touch $@

data/%.label.gz: data/%.lambda.gz
	@touch $@

data/%.causal.gz: data/%.lambda.gz
	@touch $@

reset:
	@rm -rf data result
