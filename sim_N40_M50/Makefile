IND := 40
NGENE := 5000
NCAUSAL := 50
NCOVAR := 1
NCELL_IND := 50
rhoA := 6.25 # Estimated from Microglia
rhoB := 6.25 # mean = 1, var = 0.16

RSEED := $(shell seq 1 30)

all: $(foreach r, $(RSEED), $(foreach pa, 1 2 3, $(foreach p0, 5, $(foreach p1, 1 2 225 25 275 3 35 4, result/$(p1)_$(p0)_$(pa)_$(r).eval))))

##############
# evaluation #
##############

# % = $(pve_1)_$(pve_0)_$(pve_a)_$(rseed)
result/%.eval: result/%-clean.mu_cols.gz result/%-clean.mu.gz result/%.mu_cols.gz result/%.ln_obs_mu.gz result/%.ln_resid_mu.gz data/%.causal.gz data/%.label.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@Rscript --vanilla ../R/evalute.R $* $@
	@rm -f result/$**.*.gz
	@rm -f data/$**.*.gz
	@rm -f data/$**.*.gz.index

####################
# model estimation #
####################

result/%-clean.mu.gz: data/%-clean.mtx.gz data/%.cols.gz result/%.annot.gz result/%-clean.lab.gz data/%-clean.ind.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mmutil_aggregate_col --mtx data/$*-clean.mtx.gz --col data/$*.cols.gz --annot result/$*.annot.gz --lab result/$*-clean.lab.gz --ind data/$*-clean.ind.gz --verbose --out $(shell echo $@ | sed 's/.mu.gz//g')

result/%.ln_resid_mu.gz: data/%.mtx.gz data/%.cols.gz result/%.annot.gz result/%.lab.gz result/%.trt.gz data/%.ind.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mmutil_cfa_col --mtx data/$*.mtx.gz --col data/$*.cols.gz --annot result/$*.annot.gz --trt result/$*.trt.gz --lab result/$*.lab.gz --ind data/$*.ind.gz --verbose --knn 100 --nboot 0 --rank 10 --log_scale --out $(shell echo $@ | sed 's/.ln_resid_mu.gz//g')

result/%-clean.mu_cols.gz: result/%-clean.mu.gz
	@touch $@

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

data/%-clean.mtx.gz: data/%.mu-clean.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mmutil_simulate --mu $< --rho data/$*.rho.gz --col data/$*.cols.gz --out data/$*-clean

data/%.mtx.gz: data/%.mu.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mmutil_simulate --mu $< --rho data/$*.rho.gz --col data/$*.cols.gz --out data/$*

data/%.mu.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@Rscript --vanilla ../R/simulate_glm.R $(IND) $(NGENE) $(NCAUSAL) $(NCOVAR) $(NCELL_IND) $(shell echo $* | awk -F'_' '{ for(j=1; j<=(NF - 1); ++j) print " 0." $$j; print $$NF }') $(rhoA) $(rhoB) data/$*

data/%.mu-clean.gz: data/%.mu.gz
	@touch $@

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
