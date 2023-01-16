# czyszczenie obszaru roboczego

gc(full = T)
rm(list = ls())
.rs.restartR()

# ustawienie ziarna losowania

set.seed(1)

# wczytanie pliku z potrzebny funkcjami do obliczen

source('user_functions_rps.R')

# istalowanie, ladowanie lub wylodowanie pakietow

install.packages('data.table')
require(data.table)
detach('package:data.table', unload = T)

install.packages('PerformanceAnalytics')
require(PerformanceAnalytics)
detach('package:PerformanceAnalytics', unload = T)

install.packages('MonteCarlo')
require(MonteCarlo)
detach('package:MonteCarlo', unload = T)

install.packages('ggplot2')
require(ggplot2)
detach('package:ggplo2', unload = T)

# zapisanie obliczen

save.image('obliczenia_24082019.RData')

# wczytanie poprzednich obliczen

load('obliczenia_24082019.RData')

# zapisanie danych w plikach RDS

saveRDS(file1, file = 'file1.RDS')
saveRDS(file2, file = 'file2.RDS')
saveRDS(file3, file = 'file3.RDS')

# wczytanie danych z plików RDS

readRDS(file = 'file1.RDS')
readRDS(file = 'file2.RDS')
readRDS(file = 'file3.RDS')