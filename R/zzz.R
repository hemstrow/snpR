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
  
  
  
  
  packageStartupMessage("Support for non-bialleic markers is currently under development.\n\nTo allow for non-biallilic markers, the internal structure of snpRdata objects has been changed slightly. While this version has been written to be backwards compatable, some issues may persist.\n\nIf you encounter any issues relating to a 'bi-allelic' or 'ploidy' slot, please report these and remake the snpRdata object using, for example:\n\t x <- import.snpR.data(genotypes(x), snp.meta(x), sample.meta(x), mDat = 'NN')\n")
  
}
