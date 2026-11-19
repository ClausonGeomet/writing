args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3L) stop('usage: read_saved.R REPO ROOT OUTPUT')
if (!all(grepl('^/', args))) stop('REPO, ROOT, OUTPUT must be absolute')
REPO <- normalizePath(args[[1L]], mustWork=TRUE)
ROOT <- normalizePath(args[[2L]], mustWork=TRUE)
OUTPUT <- args[[3L]]
if (any(strsplit(OUTPUT, '/', fixed=TRUE)[[1L]] %in% c('.', '..'))) stop('OUTPUT must not contain dot path components')
missing_parts <- character()
ancestor <- OUTPUT
while (!file.exists(ancestor) && !dir.exists(ancestor)) {
  link <- Sys.readlink(ancestor)
  if (!is.na(link) && nzchar(link)) stop('OUTPUT has a dangling symlink parent')
  missing_parts <- c(basename(ancestor), missing_parts)
  ancestor <- dirname(ancestor)
}
OUTPUT <- do.call(file.path, as.list(c(normalizePath(ancestor, mustWork=TRUE), missing_parts)))
if (OUTPUT %in% c(REPO, ROOT) || startsWith(OUTPUT, paste0(REPO, '/')) || startsWith(OUTPUT, paste0(ROOT, '/'))) stop('OUTPUT must be outside inputs')
dir.create(OUTPUT, recursive=TRUE, showWarnings=FALSE)
rlib <- normalizePath(file.path(ROOT, 'Rlib'), mustWork=TRUE)
.libPaths(c(rlib, .Library))
for (pkg in c('ranger','jsonlite')) {
  if (!requireNamespace(pkg, quietly=TRUE)) stop(pkg, ' missing from admitted ROOT/Rlib')
  loaded_dir <- normalizePath(getNamespaceInfo(asNamespace(pkg), 'path'), mustWork=TRUE)
  if (!identical(dirname(loaded_dir), rlib)) stop('package outside admitted ROOT/Rlib: ', pkg)
}
parts <- c('Chalcopyrite','Chalcocite','Pyrite','Quartz','Kaolinite')
E <- file.path(REPO, 'examples/publication_sc020/artifacts/SC020-v1-reproduce-002/root/evidence/SC020-v1')
registry <- jsonlite::fromJSON(file.path(ROOT,'registry/SC020-v1.json'), simplifyVector=FALSE)
id <- 'SC020-v1-LF01-R-y_nonlinear-select_lr_ml_rf'
entry <- Filter(function(z) identical(z$model_id,id), registry$models)[[1L]]
stopifnot(identical(entry$response,'y_nonlinear'), startsWith(entry$root_object_path,'models/SC020-v1/'))
model_file <- file.path(REPO,'examples/publication_sc020',entry$root_object_path)
hash_line <- system2('sha256sum', shQuote(model_file), stdout=TRUE)
stopifnot(identical(strsplit(hash_line,' ')[[1L]][1L],entry$object_sha256))
wrapper <- readRDS(model_file)
fit <- wrapper$fit
stopifnot(identical(unname(lapply(fit$selected_pairs,unname)),unname(lapply(entry$selected_pairs,unlist))))
source(file.path(REPO,'src/R/lr_ml.R'))
fixture <- read.csv(file.path(REPO,'examples/copper/copper.csv'), check.names=FALSE)
stopifnot(nrow(fixture)==240L, identical(names(fixture[,parts,drop=FALSE]),parts), all(fixture[,parts]>0), all(abs(rowSums(fixture[,parts])-1)<=1e-10), sum(fixture$split=='train')==180L, sum(fixture$split=='test')==60L)
ids <- paste0('W07_sample_',c(4,9,10,1,2,3))
rows <- fixture[match(ids,fixture$sample_id),,drop=FALSE]
stopifnot(identical(rows$sample_id,ids))
x <- rows[,parts,drop=FALSE]
levels_pct <- predict_lr_ml(fit,x)
cat('R selected pairs:',paste(vapply(fit$selected_pairs,paste,collapse='/',FUN.VALUE=character(1)),collapse=', '),'\n')
cat('R first complete composition:',paste(sprintf('%s=%.6f',parts,as.numeric(x[1,])),collapse=', '),'\n')
cat(sprintf('R first recovery prediction (%%): %.6f\n',levels_pct[[1L]]))
D <- read.csv(file.path(E,'lf01_effects.csv'),check.names=FALSE,colClasses=c(dose_token='character'))
result_rows <- list(); n <- 0L
for (dose in c(0,.01,.025,.05)) {
  contrast <- paired_contrast(fit,x,delta=dose,donor='Quartz',recipient='Kaolinite')
  frozen <- D[D$language=='R' & D$dose_token==as.character(dose),]
  frozen <- frozen[match(ids,frozen$sample_id),]
  stopifnot(identical(frozen$sample_id,ids),max(abs(contrast$before-frozen$before),abs(contrast$after-frozen$after),abs(contrast$contrast-frozen$effect))<=1e-10)
  stopifnot(max(abs(contrast$contrast-(contrast$after-contrast$before))) <= if (dose==0) 1e-12 else 1e-10)
  if (dose==0) stopifnot(max(abs(contrast$contrast),abs(contrast$after-contrast$before))<=1e-12)
  n <- n+1L
  result_rows[[n]] <- data.frame(sample_id=ids,dose_token=as.character(dose),dose=dose,before_percent=contrast$before,after_percent=contrast$after,effect_pp=contrast$contrast)
}
write.csv(do.call(rbind,result_rows),file.path(OUTPUT,'lf01_predictions.csv'),row.names=FALSE)
cat('R LF01 paired rows:',nrow(do.call(rbind,result_rows)),'\n')
L <- read.csv(file.path(E,'levels.csv'),check.names=FALSE)
F <- read.csv(file.path(E,'lf03_effects.csv'),check.names=FALSE,colClasses=c(dose_token='character'))
result <- list(); n <- 0L
for (language in c('R','Python')) for (method in c('stepr_ols','stepr_rf','select_lr_ml_rf','ilr_rf','raw_rf')) {
  a <- L[L$study=='LF03' & L$language==language & L$method==method,]
  b <- F[F$language==language & F$method==method & F$dose_token=='0.05',]
  stopifnot(nrow(a)==256L,nrow(b)==256L,length(unique(a$sample_id))==256L,length(unique(b$sample_id))==256L,setequal(a$sample_id,b$sample_id),length(unique(c(a$model_id,b$model_id)))==1L)
  good <- b[b$status=='SUCCESS',]; bad <- b[b$status=='OFF_SUPPORT',]
  stopifnot(nrow(good)==235L,nrow(bad)==21L,all(is.na(bad$error)),all(is.na(bad$Delta_pred)),all(good$Delta_oracle < -1e-8))
  pred_sign <- ifelse(abs(good$Delta_pred)<=1e-8,'zero',ifelse(good$Delta_pred<0,'negative','positive'))
  n <- n+1L
  result[[n]] <- data.frame(language=language,method=method,n_level=256,n_effect=235,excluded=21,
    observed_rmse_pp=sqrt(mean((a$pred-a$observed)^2)),noiseless_rmse_pp=sqrt(mean((a$pred-a$noiseless)^2)),
    effect_rmse_pp=sqrt(mean((good$Delta_pred-good$Delta_oracle)^2)),
    negative=sum(pred_sign=='negative'),zero=sum(pred_sign=='zero'),positive=sum(pred_sign=='positive'))
}
table <- do.call(rbind,result)
write.csv(table,file.path(OUTPUT,'lf03_table.csv'),row.names=FALSE)
cat('LF03 frozen-row arithmetic:\n')
for (i in seq_len(nrow(table))) cat(sprintf('%s %s %.6f %.6f %.6f %d/%d/%d\n',table$language[i],table$method[i],table$observed_rmse_pp[i],table$noiseless_rmse_pp[i],table$effect_rmse_pp[i],table$negative[i],table$zero[i],table$positive[i]))
