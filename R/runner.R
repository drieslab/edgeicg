# Task runners: sequential (workers = 0) or mirai daemons.

# Full garbage collection between windows. Per-type temporaries are large (cells of the type x block genes) and R
# otherwise raises its collection threshold as the heap grows, so dead matrices accumulate across types and blocks.
.edgeicg_gc <- function(st) if (!isFALSE(st$gc)) invisible(gc(verbose = FALSE, full = TRUE))

# Daemons free and allocate many per-type matrices of a few MB. Without tuning, the system allocator keeps the freed
# pages: on the kidney crop a daemon's RSS grew from 3.3 to 5.1 GB over 5 blocks while R's live heap stayed at 1.3 GB
# (macOS medium-size magazine; glibc's sliding mmap threshold does the same on Linux). Malloc reads these variables
# at process start, so they are set only while the daemons are launched (inherited), then restored.
# options(edgeicg.malloc_tune = FALSE) turns this off.
.with_malloc_env <- function(expr) {
  if (isFALSE(getOption("edgeicg.malloc_tune"))) return(expr)
  set <- switch(Sys.info()[["sysname"]],
                Darwin = c(MallocMediumZone = "0"),
                Linux = c(MALLOC_MMAP_THRESHOLD_ = "131072", MALLOC_TRIM_THRESHOLD_ = "131072"),
                character())
  if (!length(set)) return(expr)
  old <- Sys.getenv(names(set), unset = NA)
  do.call(Sys.setenv, as.list(set))
  on.exit({ for (k in names(old)) if (is.na(old[[k]])) Sys.unsetenv(k) else do.call(Sys.setenv, as.list(old[k])) })
  expr
}

# $run(tasks, sweep) -> list of task results; $set(list) adds to the state; $close().
# Daemons get the state once (as a file) and load the installed edgeicg namespace, so the task function is sent
# by reference. What the daemons hold depends on the window plan: with one window, the full per-type state (each
# task is a gene block over all cells); with several, only global settings, since every task carries the row
# slices (context) of its own windows. In-memory counts are never copied to the daemons: tasks carry their block.
.edgeicg_runner <- function(workers, st) {
  if (workers < 1) {
    env <- new.env(); env$st <- st
    return(list(run = function(tasks, sweep) lapply(tasks, function(tk) .edgeicg_task(tk, env$st, sweep)),
                set = function(more) env$st[names(more)] <- more, close = function() invisible()))
  }
  if (!requireNamespace("mirai", quietly = TRUE)) stop("workers > 0 needs the mirai package")
  .with_malloc_env(mirai::daemons(workers))
  closed <- FALSE
  push <- function(obj, name) {
    f <- tempfile(fileext = ".rds"); saveRDS(obj, f)
    mirai::everywhere({                                              # runs on each daemon
      x <- readRDS(f); ge <- globalenv()
      if (nm != "st") { s <- get(".edgeicg_st", envir = ge); s[names(x)] <- x; x <- s }
      assign(".edgeicg_st", x, envir = ge); TRUE
    }, f = f, nm = name)[]
    unlink(f)
  }
  mirai::everywhere({ Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
    loadNamespace("edgeicg"); TRUE })[]
  sp <- st
  if (is.function(st$pv$slice)) sp$pv <- list(genes = st$pv$genes, kind = "shipped")
  if (!st$single) sp[c("dz", "E", "wtypes", "Wlt", "cmu")] <- NULL
  push(sp, "st"); rm(sp)
  call1 <- function(task, sweep, f) { r <- f(task, .edgeicg_st, sweep); gc(verbose = FALSE); r }
  environment(call1) <- globalenv()
  list(run = function(tasks, sweep) {
         r <- mirai::mirai_map(tasks, call1, .args = list(sweep = sweep, f = .edgeicg_task))[]
         bad <- vapply(r, mirai::is_mirai_error, TRUE)
         if (any(bad)) stop("daemon error: ", as.character(r[[which(bad)[1]]]))
         r },
       set = function(more) push(more, "more"),
       close = function() if (!closed) { mirai::daemons(0); closed <<- TRUE })
}
