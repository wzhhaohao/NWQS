# Package load hooks. Conditional S3 method registration for broom /
# generics so NWQS does not need a hard dependency on either.

.onLoad <- function(libname, pkgname) {
  if (requireNamespace("generics", quietly = TRUE)) {
    registerS3method("tidy",   "nwqs",      tidy_nwqs,      envir = asNamespace("generics"))
    registerS3method("tidy",   "nwqs_boot", tidy_nwqs_boot, envir = asNamespace("generics"))
    registerS3method("glance", "nwqs",      glance_nwqs,    envir = asNamespace("generics"))
    registerS3method("glance", "nwqs_boot", glance_nwqs_boot, envir = asNamespace("generics"))
  }
}
