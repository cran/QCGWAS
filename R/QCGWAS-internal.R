.onAttach <-
function(libname, pkgname) {
  packageStartupMessage("QCGWAS library, version 1.0-9")
  packageStartupMessage("")
  packageStartupMessage("As of 2021, QCGWAS has been replaced by GWASinspector")
  packageStartupMessage(" We strongly recommend to switch to that package,")
  packageStartupMessage(" as it is faster and accepts more variant types.")
  packageStartupMessage(" Visit http://gwasinspector.com/ for more details.")
  packageStartupMessage("")
}
