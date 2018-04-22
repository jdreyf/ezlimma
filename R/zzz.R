.onLoad <- function(libname, pkgname){
  tryCatch(rJava::.jpackage(name=pkgname), 
           error="Failed to initialize the Java Virtual Machine (JVM) for rJava.")
}