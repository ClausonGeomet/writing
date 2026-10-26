#!/usr/bin/env Rscript
# Reproduce with R_LIBS=/tmp/emc-simcodar-d9/lib Rscript research/emc-part2a/generation.R
suppressPackageStartupMessages(library(simCoDaR))
suppressPackageStartupMessages(library(jsonlite))
RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
root <- 'research/emc-part2a/fixtures'
dir.create(root, recursive=TRUE, showWarnings=FALSE)
base_parts <- c('Quartz','Hematite','Kaolinite','Forsterite')
amb_parts <- c('Quartz','Anatase','Rutile')
lib <- mineral_assay_library()
lib$Anatase <- list(formula='TiO2', atoms=c(Ti=1,O=2))
lib$Rutile <- list(formula='TiO2', atoms=c(Ti=1,O=2))
make_draw <- function(parts, centers, seed, library) {
  colnames(centers) <- parts
  rownames(centers) <- paste0('domain_',seq_len(nrow(centers)))
  spec <- new_assay_scenario(centers, weights=rep(1/nrow(centers),nrow(centers)),logratio_sd=.25)
  simulate_mineral_assay_scenario(spec,60L,seed=seed,library=library,iron_reporting='total_Fe2O3')
}
response <- function(parts, channels, library) {
  identity <- diag(length(parts)); dimnames(identity) <- list(parts,parts)
  assay <- minerals_to_assay(identity,library,basis='fraction',iron_reporting='total_Fe2O3')
  t(as.matrix(assay$oxides_reported[,channels,drop=FALSE]))/100
}
base_c <- response(base_parts,c('SiO2','Fe2O3_total','Al2O3','MgO'),lib)
amb_c <- response(amb_parts,c('SiO2','TiO2'),lib)
write.csv(data.frame(analyte=rownames(base_c),base_c,check.names=FALSE),file.path(root,'C_base.csv'),row.names=FALSE)
write.csv(data.frame(analyte=rownames(amb_c),amb_c,check.names=FALSE),file.path(root,'C_ambiguous.csv'),row.names=FALSE)
meta <- list(schema='emc-feasibility-1',dependency=list(name='simCoDaR',commit='d9a0b27676a0cb3a4acb300a3f6b2ebc845d205d',version=as.character(packageVersion('simCoDaR')),R=R.version.string,RNGkind=RNGkind()),basis='whole-sample mineral mass fraction; oxide wt percent divided by 100 once',cases=list(),C_true=list(base=unname(base_c),ambiguous=unname(amb_c)),channels=list(base=rownames(base_c),ambiguous=rownames(amb_c)),phases=list(base=base_parts,ambiguous=amb_parts))
for (r in 0:2) {
  base <- make_draw(base_parts,rbind(c(.40,.20,.20,.20),c(.20,.40,.15,.25)),20261026+100*r,lib)
  amb <- make_draw(amb_parts,matrix(c(.6,.3,.1),nrow=1),20261027+100*r,lib)
  for (sc in c('easy','ambiguous','noisy','missing_phase','wrong_chemistry','infeasible')) {
    d <- if(sc=='ambiguous') amb else base
    parts <- if(sc=='ambiguous') amb_parts else base_parts
    channels <- if(sc=='ambiguous') rownames(amb_c) else rownames(base_c)
    C <- if(sc=='ambiguous') amb_c else base_c
    M <- as.matrix(d$truth$mineral_fraction[,parts,drop=FALSE])
    A0 <- M %*% t(C)
    ideal <- as.matrix(d$observed$xrf_oxide_wt_percent[,channels,drop=FALSE])/100
    if(max(abs(A0-ideal))>1e-10) stop('converter forward mismatch')
    ids <- as.character(d$sample_id)
    out <- data.frame(sample_id=ids,domain=as.character(d$domains[ids,'domain']),check.names=FALSE)
    for(j in seq_along(parts)) out[[paste0('m_',parts[j])]] <- M[,j]
    for(j in seq_along(channels)) {
      out[[paste0('a0_',channels[j])]] <- A0[,j]
      out[[paste0('ideal_',channels[j])]] <- ideal[,j]
    }
    file <- sprintf('%s/%s_r%d.csv',root,sc,r)
    write.csv(out,file,row.names=FALSE,quote=FALSE)
    d$metadata$scenario_assumptions$spec <- unclass(d$metadata$scenario_assumptions$spec)
    audit <- list(metadata=d$metadata,converter=d$truth$assay)
    write_json(audit,sprintf('%s/%s_r%d_audit.json',root,sc,r),auto_unbox=TRUE,digits=17,pretty=TRUE,dataframe='rows')
    meta$cases[[paste0(sc,'_r',r)]] <- list(file=basename(file),seed=if(sc=='ambiguous') 20261027+100*r else 20261026+100*r,phase_names=parts,analytes=channels,source='simCoDaR public API',audit=paste0(sc,'_r',r,'_audit.json'))
  }
}
write_json(meta,file.path(root,'source_manifest.json'),auto_unbox=TRUE,digits=17,pretty=TRUE)
