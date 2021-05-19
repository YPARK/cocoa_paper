IND := 40
NGENE := 10000
NCAUSAL := 100
NCOVAR_S := 5
NCOVAR_B := 0
NCELL_IND := 20
rhoA := 6.25 # Estimated from Microglia
rhoB := 6.25 # mean = 1, var = 0.16

RSEED := $(shell seq 1 10)

all: $(foreach r, $(RSEED), $(foreach pa, 5, $(foreach p0, 5, $(foreach p1, 3 2 15 1 05, summary/$(p1)_$(p0)_$(pa)_$(r).eval.gz))))

##############
# evaluation #
##############

# % = $(pve_1)_$(pve_0)_$(pve_a)_$(rseed)
summary/%.eval.gz: result/%.mu_cols.gz result/%.sum.gz result/%.mean.gz result/%.ln_resid_mu.gz data/%.causal.gz data/%.label.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@Rscript --vanilla ../R/evalute_v3.R $* $@
	@rm -f result/$**.*.gz
	@rm -f data/$**.*.gz
	@rm -f data/$**.*.gz.index

####################
# model estimation #
####################

result/%.sum.gz: data/%.mtx.gz data/%.cols.gz result/%.annot.gz result/%.lab.gz data/%.ind.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mmutil_aggregate_col --mtx data/$*.mtx.gz --col data/$*.cols.gz --annot result/$*.annot.gz --lab result/$*.lab.gz --ind data/$*.ind.gz --verbose --out $(shell echo $@ | sed 's/.sum.gz//g')

result/%.mean.gz: result/%.sum.gz
	@touch $@

result/%.ln_resid_mu.gz: data/%.mtx.gz data/%.cols.gz result/%.annot.gz result/%.lab.gz result/%.trt.gz data/%.ind.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mmutil_cfa_col --mtx data/$*.mtx.gz --col data/$*.cols.gz --annot result/$*.annot.gz --trt result/$*.trt.gz --lab result/$*.lab.gz --ind data/$*.ind.gz --verbose --knn 100 --rank 50 --nboot 0 --gamma_a0 1 --gamma_b0 1 --log_scale --out $(shell echo $@ | sed 's/.ln_resid_mu.gz//g')

result/%.ln_obs_mu.gz: result/%.mu_cols.gz
	@touch $@

result/%.mu_cols.gz: result/%.ln_resid_mu.gz
	@touch $@

result/%.annot.gz: data/%.cols.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@gzip -cd $< | awk '{ print $$1 FS "$*" }' | gzip -c > $@

result/%.lab.gz:  # just a single cell type
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@echo $* | gzip > $@

result/%.trt.gz: data/%.ind.gz data/%.label.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@gzip -cd $< | awk 'BEGIN { j=0; while(("gzip -cd data/$*.label.gz " | getline line)) { v=sprintf("Ind%02d",(++j)); d[v]=line } } { print d[$$1] }' | gzip -c > $@

###############################
# simulate matrix market data #
###############################

data/%.mtx.gz: data/%.mu.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mmutil_simulate --mu $< --rho data/$*.rho.gz --col data/$*.cols.gz --out data/$*

data/%.mu.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@Rscript --vanilla ../R/simulate_glm_v3.R $(IND) $(NGENE) $(NCAUSAL) $(NCOVAR_S) $(NCOVAR_B) $(NCELL_IND) $(shell echo $* | awk -F'_' '{ for(j=1; j<=(NF - 1); ++j) print " 0." $$j; print $$NF }') $(rhoA) $(rhoB) data/$*

data/%.cols.gz: data/%.mtx.gz
	@touch $@

data/%.ind.gz: data/%.mtx.gz
	@touch $@

data/%.label.gz: data/%.mu.gz
	@touch $@

data/%.causal.gz: data/%.mu.gz
	@touch $@

reset:
	@rm -rf data result
