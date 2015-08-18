args<-commandArgs(TRUE)
str(ip <- installed.packages(priority = "high"))
ip_trans = t(ip)
header = c("Package","LibPath","Version","Priority","Depends","Imports","LinkingTo",
           "Suggests","Enhances","License","License_is_FOSS","License_restricts_use",
           "OS_type","MD5sum","NeedsCompilation","Built")
write(header, file=args[1], append = FALSE, ncolumns=16,sep = "\t")
write(ip_trans, file=args[1], append = TRUE, ncolumns=16,sep = "\t")

