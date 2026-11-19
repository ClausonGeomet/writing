#' Select grouped out-of-fold pairwise log-ratios
#'
#' @param X Named data frame or numeric matrix of strictly positive fractions.
#' @param y Finite numeric response.
#' @param groups Group identifier, with each group in one fold.
#' @param fold_id Integer fold identifier per row.
#' @param max_ratios Maximum selected ratios, at most three.
#' @param min_improvement Minimum relative out-of-fold MSE improvement.
#' @param n_trees Number of ranger trees, from 1 to 100.
#' @param seed User-controlled forest seed.
#' @return A fitted `lr_ml` list.
#' @export
select_lr_ml <- function(X, y, groups, fold_id, max_ratios = 3,
                         min_improvement = 0.01, n_trees = 50, seed = 104729) {
  v <- .lr_composition(X)
  n <- nrow(v)
  d <- ncol(v)
  if (n < 6L) stop("at least six rows are required")
  if (!.lr_integer(max_ratios) || max_ratios < 0L ||
      max_ratios > min(3L, d * (d - 1L) / 2L)) stop("invalid max_ratios")
  if (!.lr_integer(n_trees) || n_trees < 1L || n_trees > 100L)
    stop("n_trees must be an integer in [1, 100]")
  if (!.lr_integer(seed) || seed < 0L || seed > .Machine$integer.max - 1000L)
    stop("seed must be a nonnegative integer")
  if (!is.numeric(min_improvement) || length(min_improvement) != 1L ||
      !is.finite(min_improvement) || min_improvement < 0 || min_improvement > 1)
    stop("min_improvement must be finite and in [0, 1]")
  if (!is.numeric(y) || length(y) != n || any(!is.finite(y)))
    stop("y must contain one finite numeric value per row")
  if (length(groups) != n || length(fold_id) != n || anyNA(groups) ||
      anyNA(fold_id) || !is.numeric(fold_id) ||
      any(!is.finite(fold_id)) || any(fold_id != floor(fold_id)))
    stop("groups and integer fold_id must have one nonmissing value per row")
  folds_by_group <- split(fold_id, as.character(groups))
  if (any(vapply(folds_by_group, function(z) length(unique(z)) != 1L, logical(1))))
    stop("a group spans multiple folds")
  folds <- unique(fold_id)
  if (any(seed + folds < 0 | seed + folds > .Machine$integer.max))
    stop("seed plus fold_id must be a valid forest seed")
  if (length(folds) < 2L || any(vapply(folds, function(z) {
    sum(fold_id == z) < 2L || sum(fold_id != z) < 2L
  }, logical(1)))) stop("at least two folds with two training and validation rows are required")
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("ranger is required; install it explicitly before fitting")
  parts <- colnames(v)
  pairs <- unlist(lapply(2:d, function(j) lapply(1:(j - 1L),
    function(i) c(parts[i], parts[j]))), recursive = FALSE)
  names(pairs) <- vapply(pairs, paste, character(1), collapse = "/")
  log_v <- log(v)
  ratios <- vapply(pairs, function(p) log_v[, p[1L]] - log_v[, p[2L]], numeric(n))
  colnames(ratios) <- names(pairs)
  fit_forest <- function(x, response, forest_seed) {
    ranger::ranger(x = as.data.frame(x), y = response, num.trees = n_trees,
      mtry = ncol(x), min.node.size = 5L, replace = TRUE,
      sample.fraction = 1, num.threads = 1L, seed = forest_seed)
  }
  oof <- function(chosen) {
    predictions <- numeric(n)
    for (fold in folds) {
      train <- fold_id != fold
      valid <- !train
      if (!length(chosen)) {
        predictions[valid] <- mean(y[train])
      } else {
        model <- fit_forest(ratios[train, chosen, drop = FALSE], y[train], seed + fold)
        predictions[valid] <- predict(model, data = as.data.frame(
          ratios[valid, chosen, drop = FALSE]), num.threads = 1L)$predictions
      }
    }
    mean((y - predictions)^2)
  }
  old_loss <- oof(character())
  chosen <- character()
  path <- list()
  for (step in seq_len(max_ratios)) {
    if (old_loss == 0) break
    best <- NA_character_
    best_loss <- Inf
    for (candidate in names(pairs)) {
      if (candidate %in% chosen) next
      loss <- oof(c(chosen, candidate))
      if (loss < best_loss - 1e-12) {
        best <- candidate
        best_loss <- loss
      }
    }
    gain <- (old_loss - best_loss) / old_loss
    accepted <- gain > 0 && gain >= min_improvement
    path[[length(path) + 1L]] <- list(step = step, pair = unname(pairs[[best]]),
      loss = best_loss, improvement = gain, accepted = accepted)
    if (!accepted) break
    chosen <- c(chosen, best)
    old_loss <- best_loss
  }
  model <- if (length(chosen)) fit_forest(ratios[, chosen, drop = FALSE], y, seed) else NULL
  fit <- list(parts = parts, selected_pairs = unname(pairs[chosen]), path = path,
    model = model, mean_y = mean(y),
    folds = unique(data.frame(group = groups, fold_id = fold_id)),
    settings = list(max_ratios = max_ratios, min_improvement = min_improvement,
      n_trees = n_trees, seed = seed, ranger_version = as.character(utils::packageVersion("ranger")),
      learner = "ranger regression; mtry=p; min.node.size=5; replace=TRUE; sample.fraction=1; num.threads=1"))
  class(fit) <- "lr_ml"
  fit
}

.lr_integer <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x == floor(x)
}

.lr_composition <- function(X, parts = NULL) {
  if (!(is.data.frame(X) || is.matrix(X)) || is.null(dim(X)) || nrow(X) < 1L)
    stop("composition must be a nonempty named data frame or matrix")
  nm <- colnames(X)
  if (is.null(nm) || anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm))
    stop("composition columns must have unique nonempty names")
  if (!is.null(parts)) {
    if (length(nm) != length(parts) || !setequal(nm, parts))
      stop("composition part names differ from training parts")
    X <- X[, parts, drop = FALSE]
  }
  if (!all(vapply(as.data.frame(X), is.numeric, logical(1))))
    stop("composition must be numeric")
  v <- as.matrix(X)
  if (ncol(v) < 2L || any(!is.finite(v)) || any(v <= 0))
    stop("composition needs at least two finite, strictly positive parts")
  if (any(abs(rowSums(v) - 1) > 1e-10)) stop("composition rows must sum to one")
  v
}

.lr_design <- function(v, parts, pairs) {
  log_v <- log(v)
  out <- vapply(pairs, function(p) log_v[, p[1L]] - log_v[, p[2L]], numeric(nrow(v)))
  out <- matrix(out, nrow = nrow(v), ncol = length(pairs))
  colnames(out) <- vapply(pairs, paste, character(1), collapse = "/")
  out
}

#' Predict from whole compositions
#' @param fit Fitted `lr_ml` result.
#' @param new_composition Complete named mineral fractions.
#' @return Numeric predictions in input row order.
#' @export
predict_lr_ml <- function(fit, new_composition) {
  if (!inherits(fit, "lr_ml")) stop("fit must be lr_ml")
  v <- .lr_composition(new_composition, fit$parts)
  if (!length(fit$selected_pairs)) return(rep(fit$mean_y, nrow(v)))
  design <- .lr_design(v, fit$parts, fit$selected_pairs)
  as.numeric(predict(fit$model, data = as.data.frame(design), num.threads = 1L)$predictions)
}

#' Finite donor-to-recipient model prediction contrast
#' @param fit Fitted `lr_ml` result.
#' @param X Complete named mineral fractions.
#' @param delta Nonnegative fraction transferred.
#' @param donor,recipient Distinct training part names.
#' @return Data frame with `before`, `after`, and `contrast`.
#' @export
paired_contrast <- function(fit, X, delta = 0.05, donor = "Quartz",
                            recipient = "Kaolinite") {
  if (!inherits(fit, "lr_ml")) stop("fit must be lr_ml")
  v <- .lr_composition(X, fit$parts)
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta < 0)
    stop("delta must be finite and nonnegative")
  if (!is.character(donor) || !is.character(recipient) ||
      length(donor) != 1L || length(recipient) != 1L ||
      is.na(donor) || is.na(recipient) || donor == recipient ||
      !donor %in% fit$parts || !recipient %in% fit$parts)
    stop("donor and recipient must be distinct training parts")
  if (any(v[, donor] <= delta)) stop("donor must remain strictly positive")
  after <- v
  after[, donor] <- after[, donor] - delta
  after[, recipient] <- after[, recipient] + delta
  before_pred <- predict_lr_ml(fit, as.data.frame(v))
  after_pred <- predict_lr_ml(fit, as.data.frame(after))
  data.frame(before = before_pred, after = after_pred,
             contrast = after_pred - before_pred)
}
