.First.lib <- function(lib,pkg)
{
   library.dynam("ljr",pkg,lib)
   cat("ljr 1.2-2 loaded\n")
}
