source("scripts/packages.R")
load(file = "/workdir/DOCKER_DUMP_V6/global-Rdata/CFSonly.rlt.RData")

CFSonly.rlt -> sobj
rm(CFSonly.rlt)

sobj$library <- plyr::mapvalues(
  x = sobj$orig.ident,
  from = c("R1744short","R1778long","R1778short","R1798long","R1798short","R1855long","R1855short","R1858long","R1858short","R2002short","R2015short","R2024short","R2033short","R2148short","R2175short","R2183short","R2225short","R2227short","R2258short","R2339short","R4511long","R4531long","R6016long","R6026long","R6088long","R6098long","R7120long","R7120short","R7317long","R7317short","R7327long","R7327short","R7385long","R7385short","R7472long","R7472short","R7513long","R7513short"),
  to = c("R1744","R1778","R1778","R1798","R1798","R1855","R1855","R1858","R1858","R2002","R2015","R2024","R2033","R2148","R2175","R2183","R2225","R2227","R2258","R2339","R4511","R4531","R6016","R6026","R6088","R6098","R7120","R7120","R7317","R7317","R7327","R7327","R7385","R7385","R7472","R7472","R7513","R7513")
)

sobj$readlength <- plyr::mapvalues(
  x = sobj$orig.ident,
  from = c("R1744short","R1778long","R1778short","R1798long","R1798short","R1855long","R1855short","R1858long","R1858short","R2002short","R2015short","R2024short","R2033short","R2148short","R2175short","R2183short","R2225short","R2227short","R2258short","R2339short","R4511long","R4531long","R6016long","R6026long","R6088long","R6098long","R7120long","R7120short","R7317long","R7317short","R7327long","R7327short","R7385long","R7385short","R7472long","R7472short","R7513long","R7513short"),
  to = c("short","long","short","long","short","long","short","long","short","short","short","short","short","short","short","short","short","short","short","short","long","long","long","long","long","long","long","short","long","short","long","short","long","short","long","short","long","short")
)

sobj$batchID <- plyr::mapvalues(
  x = sobj$orig.ident,
  from = c("R1744short","R1778long","R1778short","R1798long","R1798short","R1855long","R1855short","R1858long","R1858short","R2002short","R2015short","R2024short","R2033short","R2148short","R2175short","R2183short","R2225short","R2227short","R2258short","R2339short","R4511long","R4531long","R6016long","R6026long","R6088long","R6098long","R7120long","R7120short","R7317long","R7317short","R7327long","R7327short","R7385long","R7385short","R7472long","R7472short","R7513long","R7513short"),
  to = c("2","18","18","18","18","18","18","18","18","3","3","2","2","2","3","3","3","3","2","2","13","13","13","13","13","13","16","16","16","16","16","16","16","16","18","18","18","18")
)

sobj$sex <- plyr::mapvalues(
  x = sobj$orig.ident,
  from = c("R1744short","R1778long","R1778short","R1798long","R1798short","R1855long","R1855short","R1858long","R1858short","R2002short","R2015short","R2024short","R2033short","R2148short","R2175short","R2183short","R2225short","R2227short","R2258short","R2339short","R4511long","R4531long","R6016long","R6026long","R6088long","R6098long","R7120long","R7120short","R7317long","R7317short","R7327long","R7327short","R7385long","R7385short","R7472long","R7472short","R7513long","R7513short"),
  to = c("F","M","M","M","M","M","M","M","M","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","M","M","M","M")
)

sobj$phenotypeID <- plyr::mapvalues(
  x = sobj$orig.ident,
  from = c("R1744short","R1778long","R1778short","R1798long","R1798short","R1855long","R1855short","R1858long","R1858short","R2002short","R2015short","R2024short","R2033short","R2148short","R2175short","R2183short","R2225short","R2227short","R2258short","R2339short","R4511long","R4531long","R6016long","R6026long","R6088long","R6098long","R7120long","R7120short","R7317long","R7317short","R7327long","R7327short","R7385long","R7385short","R7472long","R7472short","R7513long","R7513short"),
  to = c("1","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","1","1","1","2","2","2","1","1","1","1","1","1","2","2","2","2","1","1","1","1","1","1")
)

sobj$day <- plyr::mapvalues(
  x = sobj$orig.ident,
  from = c("R1744short","R1778long","R1778short","R1798long","R1798short","R1855long","R1855short","R1858long","R1858short","R2002short","R2015short","R2024short","R2033short","R2148short","R2175short","R2183short","R2225short","R2227short","R2258short","R2339short","R4511long","R4531long","R6016long","R6026long","R6088long","R6098long","R7120long","R7120short","R7317long","R7317short","R7327long","R7327short","R7385long","R7385short","R7472long","R7472short","R7513long","R7513short"),
  to = c("D2","D1","D1","D2","D2","D2","D2","D1","D1","D1","D2","D1","D2","D1","D1","D2","D2","D1","D1","D2","D1","D2","D1","D1","D2","D2","D1","D1","D1","D1","D2","D2","D2","D2","D1","D1","D2","D2")
)

sobj$batchID <- plyr::mapvalues(
  x = sobj$orig.ident,
  from = c("R1744short","R1778long","R1778short","R1798long","R1798short","R1855long","R1855short","R1858long","R1858short","R2002short","R2015short","R2024short","R2033short","R2148short","R2175short","R2183short","R2225short","R2227short","R2258short","R2339short","R4511long","R4531long","R6016long","R6026long","R6088long","R6098long","R7120long","R7120short","R7317long","R7317short","R7327long","R7327short","R7385long","R7385short","R7472long","R7472short","R7513long","R7513short"),
  to = c("02","18","18","18","18","18","18","18","18","03","03","02","02","02","03","03","03","03","02","02","13","13","13","13","13","13","16","16","16","16","16","16","16","16","18","18","18","18")
)

#useful for grouping/splitting
sobj$batch.readlen <- plyr::mapvalues(
  x = sobj$orig.ident,
  from = c("R1744short","R1778long","R1778short","R1798long","R1798short","R1855long","R1855short","R1858long","R1858short","R2002short","R2015short","R2024short","R2033short","R2148short","R2175short","R2183short","R2225short","R2227short","R2258short","R2339short","R4511long","R4531long","R6016long","R6026long","R6088long","R6098long","R7120long","R7120short","R7317long","R7317short","R7327long","R7327short","R7385long","R7385short","R7472long","R7472short","R7513long","R7513short"),
  to = c("early","late.long","late.short","late.long","late.short","late.long","late.short","late.long","late.short","early","early","early","early","early","early","early","early","early","early","early","early","early","early","early","early","early","late.long","late.short","late.long","late.short","late.long","late.short","late.long","late.short","late.long","late.short","late.long","late.short")
)

sobj.meta = sobj@meta.data
saveRDS(sobj, "data/filtered_seurat_object_with_metaData_for_RLT.Rds")
