.First.lib <- function(lib,pkg)
{
   library.dynam("ljr",pkg,lib)
   cat("ljr 1.1-0 loaded\n")
}
