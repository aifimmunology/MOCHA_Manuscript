### Runtime analyses
require(readxl)
require(data.table)

mocha = as.data.table(read_xlsx('~/Downloads/MOCHA_runtime_comparison.xlsx',sheet = 1))
homer = as.data.table(read_xlsx('~/Downloads/MOCHA_runtime_comparison.xlsx',sheet = 2))
macs2 = as.data.table(read_xlsx('~/Downloads/MOCHA_runtime_comparison.xlsx',sheet = 3))


mocha[Ncells %in% c(100,1000, 10000), list(mean(Time), sd(Time)), by=Ncells]
homer[Ncells %in% c(100,1000, 10000), list(mean(Time), sd(Time)), by=Ncells]
macs2[Ncells %in% c(100,1000, 10000), list(mean(Time), sd(Time)), by=Ncells]

ks.test(mocha$Time, macs2$Time)
ks.test(mocha$Time, homer$Time)

summary(macs2$Time/mocha$Time)
summary(homer$Time/mocha$Time)
