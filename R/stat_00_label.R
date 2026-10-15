# =============================================================================
# Group label access layer
# =============================================================================
# CONTRACT
#   raw.data / clean.data / scale.data  ->  FEATURES ONLY. The label column must
#                                          never appear in any of them.
#   info.data                           ->  metadata; the ONE home of the label.
#   group_col (slot)                    ->  the column NAME, resolved in
#                                          info.data.
#
# WHY A DEDICATED LAYER
#   The previous contract required the label to live inside the feature matrix,
#   so every consumer had to remember to exclude it by name. That convention
#   had already leaked at three sites (PlotPCA, PlotCorrelationHeatmap,
#   PlotFeatureHeatmap), each of which silently fed a numeric 0/1 label into an
#   unsupervised computation -- e.g. a numeric group column made PC1 group
#   separation look ~5x stronger than the features actually justified.
#
#   Here the robustness comes from the data layout instead of from discipline:
#   a feature matrix physically cannot leak the label, so no call site has to
#   remember anything. `.stat_features()` therefore *asserts* the invariant
#   rather than quietly repairing a violation with setdiff().
#
# All functions here are internal; they are not exported.
# =============================================================================


#' Read the group column NAME from an object, tolerating an unset slot
#'
#' Returns NULL (rather than character(0)) when the slot is absent or empty, so
#' callers can use a single `is.null()` test.
#'
#' @param object A Stat / Train_Model / Subtyping / PrognosiX object.
#' @return A single character string, or NULL.
#' @keywords internal
.stat_group_col <- function(object) {
  if (is.null(object)) return(NULL)
  gc <- tryCatch(methods::slot(object, "group_col"), error = function(e) NULL)
  if (is.null(gc)) return(NULL)
  if (length(gc) != 1) return(NULL)
  if (is.na(gc) || !nzchar(as.character(gc))) return(NULL)
  as.character(gc)
}


#' Read `info.data` from any icare object, returning an empty frame if absent
#'
#' @param object An icare object.
#' @return A data.frame (possibly 0 x 0).
#' @keywords internal
.stat_info <- function(object) {
  info <- tryCatch(methods::slot(object, "info.data"), error = function(e) NULL)
  if (is.null(info) || !is.data.frame(info)) return(data.frame())
  info
}


#' Are a frame's row names the default integer sequence?
#'
#' Two frames that both use default row names can only be aligned positionally.
#' That is almost never what the user means when one side carries sample IDs,
#' so the label layer treats it as an error rather than guessing.
#'
#' @param x A data.frame.
#' @return Logical scalar.
#' @keywords internal
.stat_has_default_rownames <- function(x) {
  if (!is.data.frame(x) || nrow(x) == 0) return(TRUE)
  rn <- rownames(x)
  is.null(rn) || identical(rn, as.character(seq_len(nrow(x))))
}


#' Pick the frame that currently defines the analysis sample set
#'
#' `clean.data` is the working frame once cleaning has run; before that it is
#' empty and `raw.data` is authoritative. Every consumer needs the same answer to
#' "which samples are we analysing", so the choice lives here rather than being
#' re-decided at each call site.
#'
#' @param object An icare object.
#' @return A data.frame (possibly 0 x 0).
#' @keywords internal
.stat_primary_frame <- function(object) {
  clean <- .stat_features(object, which = "clean", check = FALSE)
  if (nrow(clean) > 0) return(clean)
  .stat_features(object, which = "raw", check = FALSE)
}


#' Extract the group label column, aligned to a given set of sample rows
#'
#' The single read path for group labels. Resolves `group_col` in `info.data`
#' and reorders by `rows`, so the result is always aligned to whatever row set
#' the caller is working with. Every failure mode raises an error: this function
#' never returns NULL and never degrades to an ungrouped analysis, because a
#' silent fallback is how the previous contract produced wrong-but-plausible
#' figures.
#'
#' @param object An icare object whose `group_col` slot names a column of `info.data`.
#' @param rows Character vector of sample IDs to align to. Defaults to
#'   `rownames(object@clean.data)`.
#' @param group_col Optional explicit column name, overriding the slot.
#' @return A one-column data.frame with the labels, row-named by `rows`.
#' @keywords internal
.stat_labels <- function(object, rows = NULL, group_col = NULL) {
  gc <- if (is.null(group_col)) .stat_group_col(object) else as.character(group_col)

  if (is.null(gc)) {
    stop("This object has no group_col. Create it with ",
         "CreateStatObject(..., group_col = \"<name>\") so the label column is ",
         "recorded in info.data.", call. = FALSE)
  }

  if (is.null(rows)) {
    rows <- rownames(.stat_primary_frame(object))
  }

  info <- .stat_info(object)

  if (nrow(info) == 0) {
    stop(sprintf("info.data is empty, so label column '%s' cannot be found. ", gc),
         "Pass the label through CreateStatObject(info.data = <frame containing '",
         gc, "'>).", call. = FALSE)
  }
  if (!gc %in% colnames(info)) {
    stop(sprintf("Column '%s' is not in info.data (columns present: %s). ", gc,
                 paste(colnames(info), collapse = ", ")),
         sprintf("group_col is resolved against info.data, not against the feature matrix."),
         call. = FALSE)
  }

  if (is.null(rows)) {
    if (.stat_has_default_rownames(info)) {
      stop(sprintf("Cannot align labels for '%s': ", gc),
           "info.data has default row names, so there is no reliable sample key. ",
           "Give both the features and info.data explicit row names (e.g. Case IDs).",
           call. = FALSE)
    }
    rows <- rownames(info)
  }

  missing_rows <- setdiff(rows, rownames(info))
  if (length(missing_rows) > 0) {
    stop(sprintf("info.data is missing %d of the %d sample row names, so labels ",
                 length(missing_rows), length(rows)),
         sprintf("cannot be aligned (e.g. %s). ", paste(head(missing_rows, 3), collapse = ", ")),
         "info.data must retain every sample that clean.data holds.",
         call. = FALSE)
  }

  out <- info[rows, gc, drop = FALSE]
  rownames(out) <- rows
  out
}


#' Extract a feature matrix, asserting the label column is absent
#'
#' @param object An icare object.
#' @param which Which frame: `"clean"` (default), `"raw"`, or `"scale"`.
#' @param check When `FALSE`, skip the label-absence assertion (used by
#'   `.stat_labels()` to resolve its default `rows` without recursing).
#' @return A data.frame of features.
#' @keywords internal
.stat_features <- function(object, which = c("clean", "raw", "scale"), check = TRUE) {
  which <- match.arg(which)

  data <- switch(
    which,
    clean = tryCatch(methods::slot(object, "clean.data"), error = function(e) NULL),
    raw   = tryCatch(methods::slot(object, "raw.data"),   error = function(e) NULL),
    scale = tryCatch(methods::slot(object, "scale.data"), error = function(e) NULL)
  )

  if (is.null(data)) data <- data.frame()
  if (!is.data.frame(data)) data <- as.data.frame(data)

  if (check) {
    gc <- .stat_group_col(object)
    if (!is.null(gc) && gc %in% colnames(data)) {
      stop(sprintf("Internal consistency error: label column '%s' is present in ", gc),
           sprintf("object@%s. Feature matrices must contain features only; ", which),
           "the label lives in info.data.", call. = FALSE)
    }
  }

  data
}


#' Attach the label column to a feature frame, for workers that need (X, y)
#'
#' This is the join point between the Stat-level API and the plain-data.frame
#' workers (`batch_Wilcoxon`, `gaze_analysis`, `prcomp` colouring, caret, ...).
#' Only the ~12 Stat-level wrappers call it; the workers themselves are
#' unchanged, because they still receive a single frame holding both the
#' features and the label.
#'
#' @param object An icare object.
#' @param data Feature frame to annotate. Defaults to `object@clean.data`.
#' @param group_col Optional explicit label column name.
#' @return `data` with the label column appended (unchanged if there is none).
#' @keywords internal
.stat_with_group <- function(object, data = NULL, group_col = NULL) {
  if (is.null(data)) data <- .stat_primary_frame(object)
  if (!is.data.frame(data) || nrow(data) == 0) return(data)

  gc <- if (is.null(group_col)) .stat_group_col(object) else as.character(group_col)
  if (is.null(gc)) return(data)
  if (gc %in% colnames(data)) return(data)   # already joined by the caller

  labels <- .stat_labels(object, rows = rownames(data), group_col = gc)
  out <- cbind(data, labels)
  rownames(out) <- rownames(data)
  out
}


#' Re-index info.data to a feature frame's row set
#'
#' info.data is treated as append-only: operations that drop samples from
#' clean.data leave info.data alone, and alignment happens at read time. This
#' helper exists for the opposite case -- a frame that has just gained rows --
#' and to re-establish the invariant explicitly.
#'
#' @param object An icare object (modified in place for S4, returned otherwise).
#' @param rows Target row names.
#' @return The object, with `info.data` re-indexed to `rows`.
#' @keywords internal
.stat_sync_rows <- function(object, rows) {
  info <- .stat_info(object)
  if (nrow(info) == 0) return(object)

  missing_rows <- setdiff(rows, rownames(info))
  if (length(missing_rows) > 0) {
    stop(sprintf("Cannot synchronise info.data: %d row(s) are missing from it ",
                 length(missing_rows)),
         sprintf("(e.g. %s). ", paste(head(missing_rows, 3), collapse = ", ")),
         "info.data must carry a row for every sample in the feature matrix.",
         call. = FALSE)
  }

  new_info <- info[rows, , drop = FALSE]
  rownames(new_info) <- rows

  if (methods::.hasSlot(object, "info.data")) {
    methods::slot(object, "info.data") <- new_info
  }
  object
}


#' Migrate a pre-1.1 object that still carries the label inside its frames
#'
#' Objects saved before the label moved to info.data are still readable: the
#' label is lifted out of raw/clean/scale and written to info.data. The move is
#' announced rather than performed silently, because a silent compatibility
#' shim is exactly the failure mode this layer exists to remove.
#'
#' @param object An icare object.
#' @return The migrated object.
#' @keywords internal
.stat_adopt_legacy <- function(object) {
  gc <- .stat_group_col(object)
  if (is.null(gc)) return(object)
  if (!methods::.hasSlot(object, "clean.data")) return(object)

  moved <- character()
  for (slotname in c("raw.data", "clean.data", "scale.data")) {
    d <- tryCatch(methods::slot(object, slotname), error = function(e) NULL)
    if (is.null(d) || !is.data.frame(d) || nrow(d) == 0) next
    if (!gc %in% colnames(d)) next

    info <- .stat_info(object)
    if (nrow(info) == 0) {
      info <- data.frame(row.names = rownames(d))
    } else {
      info <- info[rownames(d), , drop = FALSE]
      rownames(info) <- rownames(d)
    }
    if (!gc %in% colnames(info)) {
      info[[gc]] <- d[[gc]]
    }
    methods::slot(object, "info.data") <- info
    methods::slot(object, slotname) <- d[, colnames(d) != gc, drop = FALSE]
    moved <- c(moved, slotname)
  }

  if (length(moved) > 0) {
    message(sprintf(
      "Legacy Stat object: label column '%s' was found in %s and has been moved to info.data. ",
      gc, paste(moved, collapse = ", ")),
      "Re-create the object with CreateStatObject(..., info.data = <frame containing '",
      gc, "'>) to silence this message.")
  }
  object
}


#' Explicitly run an object's validity method
#'
#' S4 `@<-` assignment does NOT re-run a class's validity function (it only
#' checks the slot's own class), so a package that relies on validity alone
#' silently accumulates invalid objects. This package calls `validObject()`
#' explicitly at the few places where an invariant can be broken.
#'
#' @param object An S4 object.
#' @return The object, invisibly, after validation.
#' @keywords internal
.stat_validate <- function(object) {
  tryCatch(methods::validObject(object),
           error = function(e) stop("Object invariant violated: ", conditionMessage(e),
                                    call. = FALSE))
  invisible(object)
}
