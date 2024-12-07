.onAttach <- function(libname, pkgname){
  version <- utils::packageVersion("snpR")
  packageStartupMessage("snpR version: ", version, " loaded.\n", .console_hline(), "\n")
  if(grepl("\\.9[0-9]{3}$", version)){
    packageStartupMessage("This is a development build of snpR.\nWhile this version is (probably) passing all tests, there may still be bugs present.\nPlease report any issues to the github issues page (https://github.com/hemstrow/snpR/issues)!\n\n\n")
    packageStartupMessage("Check out the NEWS page on github for update notes (https://github.com/hemstrow/snpR/blob/dev/NEWS.md)\n", .console_hline(), "\n")
  }
  else{
    packageStartupMessage("Please report any issues to the github issues page (https://github.com/hemstrow/snpR/issues)!\n\n\n")
    packageStartupMessage("Check out the NEWS page on github for update notes (https://github.com/hemstrow/snpR/blob/master/NEWS.md)\n", .console_hline(), "\n")
    
  }
  packageStartupMessage("If you use snpR in your work, please cite both it AND the methods it uses! You can always run the citations() function on a snpRdata object to see references for the methods used!\n")
}
