.First.lib <- function(lib,pkg)
{
   library.dynam("ljr",pkg,lib)
   cat("ljr 1.0-1 loaded\n")
}
