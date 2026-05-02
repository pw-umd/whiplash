library(raster);library(exactextractr);library(terra);library(sf);library(trend);library(zoo)

setwd("G:/Whiplash US/")

# read gridded SPEI data 
spei03 <- rast("data/spei/nclimgrid-spei-pearson-03.nc")
spei12 <- rast("data/spei/nclimgrid-spei-pearson-12.nc")
spei_base03 <- spei03[[1:660]] # 1895-01 through 1949-12 as baseline period used to calculate the threshold
spei03 <- spei03[[361:1561]] # from 1925-01 through 2025-01 (1mo fewer when calculating the difference)
spei_base12 <- spei12[[1:660]]
spei12 <- spei12[[361:1561]]

# define the start and end year
start_years <- seq(1925, 2020, by=5)
end_years   <- start_years + 4

# calculate the total & 5-year numbers of whiplash events 
## perform the layer-wise analysis first and grid-wise analysis second
thres_all <- thres_d2w <- thres_w2d <- list()
whip.national <- data.frame(matrix(nrow=1800,ncol=5))
names(whip.national) <- c("timescale","interval","type","year","total")
whip.national$timescale <- rep(c("SPEI-3","SPEI-12"),each=900)
whip.national$interval <- rep(rep(c("5-year","10-year","20-year"),each=300),times=2)
whip.national$type <- rep(rep(c("Overall","Dry-to-wet","Wet-to-dry"),each=100),times=6)
whip.national$year <- rep(1925:2024,times=18)
for (c in c("03","12")) { # 2 SPEI timescales
  print(paste0("SPEI timescale: ",c))
  spei_base.c <- get(paste0("spei_base",c))
  spei.c <- get(paste0("spei",c))
  ## calculate difference in monthly SPEI
  spei_base_diff.c <- spei_base.c[[2:nlyr(spei_base.c)]] - spei_base.c[[1:(nlyr(spei_base.c)-1)]]
  spei_diff.c <- spei.c[[2:nlyr(spei.c)]] - spei.c[[1:(nlyr(spei.c)-1)]]
  spei_base_diff_abs.c <- abs(spei_base_diff.c) # use absolute values to identify thresholds of OVERALL whiplash
  spei_diff_abs.c <- abs(spei_diff.c)
  
  for (cc in c("05","10","20")) { # 3 recurrence intervals
    print(paste0("Recurrence interval: ",cc," years"))
    cutoff1.cc <- ifelse(cc=="05",0.9833,                # cutoff for overall and dry-to-wet whiplash
                         ifelse(cc=="10",0.9917,0.9958)) # 1-1/60 for 5-year recurrence interval
                                                         # 1-1/120 for 10-year recurrence interval
                                                         # 1-1/240 for 20-year recurrence interval
    cutoff2.cc <- ifelse(cc=="05",0.0167,                # cutoff for wet-to-dry whiplash
                         ifelse(cc=="10",0.0083,0.0042))
    # A “10-year recurrence” means a probability of 1/120 per month (≈0.0083). 
    # In other words, the threshold is the value that only about 0.83% of month-to-month changes exceed in magnitude.
    # https://www.nature.com/articles/s43017-024-00624-z
    
    ## calculate thresholds based on years 1895-1949
    thres_all[[c]][[cc]] <- app(spei_base_diff_abs.c, fun = function(x) quantile(x, cutoff1.cc, na.rm = TRUE)) # overall
    thres_d2w[[c]][[cc]] <- app(spei_base_diff.c, fun = function(x) quantile(x, cutoff1.cc, na.rm = TRUE)) # dry-to-wet
    thres_w2d[[c]][[cc]] <- app(spei_base_diff.c, fun = function(x) quantile(x, cutoff2.cc, na.rm = TRUE)) # wet-to-dry
    
    ## compare the difference during 1925-2024 with the thresholds
    whip_all.cc <- spei_diff_abs.c >= thres_all[[c]][[cc]] # SpatRaster with 1200 binary layers
    whip_d2w.cc <- spei_diff.c >= thres_d2w[[c]][[cc]]
    whip_w2d.cc <- spei_diff.c <= thres_w2d[[c]][[cc]]
    
    ## calculate total number of whiplash events during 1925-2024, 1925-1974, and 1975-2024
    whip_all_total.cc <- app(whip_all.cc, sum, na.rm = TRUE)
    whip_d2w_total.cc <- app(whip_d2w.cc, sum, na.rm = TRUE)
    whip_w2d_total.cc <- app(whip_w2d.cc, sum, na.rm = TRUE)
    whip_all_1st.cc <- app(whip_all.cc[[1:600]], sum, na.rm = TRUE)
    whip_d2w_1st.cc <- app(whip_d2w.cc[[1:600]], sum, na.rm = TRUE)
    whip_w2d_1st.cc <- app(whip_w2d.cc[[1:600]], sum, na.rm = TRUE)
    whip_all_2nd.cc <- app(whip_all.cc[[601:1200]], sum, na.rm = TRUE)
    whip_d2w_2nd.cc <- app(whip_d2w.cc[[601:1200]], sum, na.rm = TRUE)
    whip_w2d_2nd.cc <- app(whip_w2d.cc[[601:1200]], sum, na.rm = TRUE)

    ## calculate number of whiplash events every 5 years
    whip_all_5y.cc <- whip_d2w_5y.cc <- whip_w2d_5y.cc <- list()
    for (i in 1:20) {
      print(i)
      start <- (i-1)*60+1
      end <- i*60
      whip_all_5y.cc[[i]] <- app(whip_all.cc[[start:end]], sum, na.rm = TRUE)
      whip_d2w_5y.cc[[i]] <- app(whip_d2w.cc[[start:end]], sum, na.rm = TRUE)
      whip_w2d_5y.cc[[i]] <- app(whip_w2d.cc[[start:end]], sum, na.rm = TRUE)
    }
    whip_all_5y.cc <- rast(whip_all_5y.cc)
    whip_d2w_5y.cc <- rast(whip_d2w_5y.cc)
    whip_w2d_5y.cc <- rast(whip_w2d_5y.cc)
    
    ## stack for easy use - total and every 5 years during 1925-2024 & single years during 1973-2022 (for population exposure calculation)
    whip_all_all.cc <- c(whip_all_total.cc, whip_all_1st.cc, whip_all_2nd.cc, whip_all_5y.cc, whip_all.cc[[577:1176]]) # SpatRaster with 623 layers
    whip_d2w_all.cc <- c(whip_d2w_total.cc, whip_d2w_1st.cc, whip_d2w_2nd.cc, whip_d2w_5y.cc, whip_d2w.cc[[577:1176]])
    whip_w2d_all.cc <- c(whip_w2d_total.cc, whip_w2d_1st.cc, whip_w2d_2nd.cc, whip_w2d_5y.cc, whip_w2d.cc[[577:1176]])
    
    names(whip_all_all.cc) <- names(whip_d2w_all.cc) <- names(whip_w2d_all.cc) <- 
      c("total","total_1st","total_2nd",paste0("y", substr(start_years, 3, 4), "_", substr(end_years, 3, 4)),
        format(seq(as.yearmon("1973-01"), as.yearmon("2022-12"), by = 1/12), "%Y-%m"))
    
    # sum every 12 layers → gives 100-layer raster
    index <- rep(1:100, each=12)
    
    whip_all_year <- tapp(whip_all.cc, index, sum, na.rm=TRUE)
    whip_d2w_year <- tapp(whip_d2w.cc, index, sum, na.rm=TRUE)
    whip_w2d_year <- tapp(whip_w2d.cc, index, sum, na.rm=TRUE)

    tot_all <- global(whip_all_year, sum, na.rm=TRUE)[,1]
    tot_d2w <- global(whip_d2w_year, sum, na.rm=TRUE)[,1]
    tot_w2d <- global(whip_w2d_year, sum, na.rm=TRUE)[,1]
    
    years <- 1925:2024
    for (j in 1:100) {
      idx <- whip.national$timescale == c &
             whip.national$interval  == cc &
             whip.national$year      == years[j]
      
      whip.national$total[idx & whip.national$type=="all"] <- tot_all[j]
      whip.national$total[idx & whip.national$type=="d2w"] <- tot_d2w[j]
      whip.national$total[idx & whip.national$type=="w2d"] <- tot_w2d[j]
    }
    
    ## save data for next steps
    writeCDF(whip_all_all.cc, paste0("study1/results/rasters/1925~2024_base_1895~1949/count/whip_all_spei",c,"_",cc,"y.nc"), overwrite = TRUE, compression = 9)
    writeCDF(whip_d2w_all.cc, paste0("study1/results/rasters/1925~2024_base_1895~1949/count/whip_d2w_spei",c,"_",cc,"y.nc"), overwrite = TRUE, compression = 9)
    writeCDF(whip_w2d_all.cc, paste0("study1/results/rasters/1925~2024_base_1895~1949/count/whip_w2d_spei",c,"_",cc,"y.nc"), overwrite = TRUE, compression = 9)
  }
}
write.csv(whip.national, "study1/results/whip.national.csv", row.names = FALSE)


# test the trend at the grid level
## whiplash trend (spatial)
trend_fun <- function(x) {
  if (all(is.na(x))) return(c(NA, NA))
  res <- try(sens.slope(x), silent = TRUE)
  c(res$estimates, res$p.value)
}

whip.grid.trend <- list()
for (i in 1:2) { # 1st or 2nd 50 years
  idx <- ifelse(i==1, "1925_1974", "1975_2024")
  print(paste0("Time period: ",idx))
  for (c in c("03","12")) {
    print(paste0("SPEI timescale: ",c))
    for (cc in c("05","10","20")) {
      print(paste0("Recurrence interval: ",cc," years"))
      for (ccc in c("all","d2w","w2d")) {
        print(paste0("Whiplash type: ",ccc))
        whip.ccc <- rast(paste0("study1/results/rasters/1925~2024_base_1895~1949/count/whip_",ccc,"_spei",c,"_",cc,"y.nc"))[[(10*i-6):(10*i+3)]] # 1st-3rd layers are numbers for 1925-2024, 1925-1974, and 1975-2024
        
        whip.grid.trend[[c]][[cc]][[ccc]] <- app(whip.ccc, trend_fun)
        names(whip.grid.trend[[c]][[cc]][[ccc]]) <- c("slope", "pval")
        
        r.trend <- whip.grid.trend[[c]][[cc]][[ccc]][["slope"]]
        r.sig <- whip.grid.trend[[c]][[cc]][[ccc]][["pval"]]
        outdir <- paste0("study1/results/rasters/1925~2024_base_1895~1949/trend_", idx)
        writeRaster(r.trend,paste0(outdir,"/trend_",ccc,"_spei",c,"_",cc,"y.tif"),overwrite = TRUE)
        writeRaster(r.sig,paste0(outdir,"/sig_",ccc,"_spei",c,"_",cc,"y.tif"),overwrite = TRUE)
      }
    }
  }
}


## population trend
library(geodata)
usa <- geodata::gadm("USA", level=1, path=".")
conus <- usa[!usa$NAME_1 %in% c("Alaska", "Hawaii"), ]
files <- list.files("G:/Gridded population/GHSL/", full.names = TRUE)
pop2020 <- rast(files[10])
conus <- project(conus, crs(pop2020))

pop.grid <- pop.grid.trend <- list()
for (i in 1:10) {
  print((i-1)*5+1975)
  pop.grid[[i]] <- rast(files[i])
  pop.grid[[i]] <- crop(mask(pop.grid[[i]], conus), conus)
}
pop.grid <- rast(pop.grid)
pop.grid.trend <- app(pop.grid, trend_fun)
names(pop.grid.trend) <- c("slope", "pval")
pop.trend <- pop.grid.trend[["slope"]]
pop.sig <- pop.grid.trend[["pval"]]
outdir <- "study1/results/rasters/1925~2024_base_1895~1949/trend_50y/"
writeRaster(pop.trend, paste0(outdir,"trend_pop.tif"),overwrite = TRUE)
writeRaster(pop.sig, paste0(outdir,"sig_pop.tif"),overwrite = TRUE)


# population exposure
files <- list.files("G:/Gridded population/GHSL/", full.names = TRUE)

r.03 <- rast("study1/results/rasters/1925~2024_base_1895~1949/count/whip_all_spei03_10y.nc")[[1]]
usa <- geodata::gadm("USA", level=1, path=".")
conus <- usa[!usa$NAME_1 %in% c("Alaska", "Hawaii"), ]
conus <- project(conus, crs(r.03))

national.df <- data.frame(matrix(nrow=180,ncol=8))
names(national.df) <- c("timescale","interval","type","year","exp_total","pop_total","ever_exposed","months_per_person")
national.df$timescale <- rep(c("03","12"),each=90)
national.df$interval <- rep(rep(c("05","10","20"),each=30),times=2)
national.df$type <- rep(rep(c("all","d2w","w2d"),each=10),times=6)
national.df$year <- rep(1975+(1:10-1)*5,times=18)

whip.grid.5y <- permo.grid.5y <- list()
for (c in c("03","12")) {
  print(paste0("SPEI timescale: ",c))
  for (cc in c("05","10","20")) {
    print(paste0("Recurrence interval: ",cc," years"))
    for (ccc in c("all","d2w","w2d")) {
      print(paste0("Whiplash type: ",ccc))
      whip.ccc <- rast(paste0("study1/results/rasters/1925~2024_base_1895~1949/count/whip_",ccc,"_spei",c,"_",cc,"y.nc"))[[24:623]] # single years during 1973-2022
      names(whip.ccc) <- seq(as.yearmon("1973-01"), as.yearmon("2022-12"), by = 1/12)
      
      for (i in 1:10) {
        pop.i <- rast(files[i])
        pop.proj.i <- project(pop.i, r.03, method = "sum")
        pop.proj.i <- mask(pop.proj.i, conus)
        
        year.i <- 1975+(i-1)*5
        start.i <- as.yearmon(paste0(year.i-2,"-01"))
        end.i <- as.yearmon(paste0(year.i+2,"-12"))
        
        whip.ccc.i <- whip.ccc[[names(whip.ccc) >= start.i & names(whip.ccc) <= end.i]]
        
        ## FREQUENCY - national fraction of population exposed during the 5-year period
        exp.ever.i <- as.int(app(whip.ccc.i, sum, na.rm = TRUE) > 0) # binary raster indicating any whiplash in the 5-year window
        exp.pop.total.i <- global(exp.ever.i * pop.proj.i, sum, na.rm = TRUE)[1, 1] # national exposed population
        pop.total.i <- global(pop.proj.i, sum, na.rm = TRUE)[1, 1] # total population
        
        national.df[national.df$timescale==c &
                      national.df$interval==cc &
                      national.df$type==ccc &
                      national.df$year==year.i,"exp_total"] <- exp.pop.total.i
        national.df[national.df$timescale==c &
                      national.df$interval==cc &
                      national.df$type==ccc &
                      national.df$year==year.i,"pop_total"] <- pop.total.i
        national.df[national.df$timescale==c &
                      national.df$interval==cc &
                      national.df$type==ccc &
                      national.df$year==year.i,"ever_exposed"] <- exp.pop.total.i/pop.total.i
        
        ## total number of events during 5-year period
        whip.grid.5y[[c]][[cc]][[ccc]][[i]] <- sum(whip.ccc.i, na.rm = TRUE)
        
        ## INTENSITY - gridded annual person-months exposed during 1975-2020
        permo.grid.5y[[c]][[cc]][[ccc]][[i]] <- whip.grid.5y[[c]][[cc]][[ccc]][[i]] * pop.proj.i / 5
        
        ## INTENSITY - national annual average number of months per person exposed
        exp.total.i <- global(permo.grid.5y[[c]][[cc]][[ccc]][[i]], sum, na.rm = TRUE)[1, 1]
        national.df[national.df$timescale==c &
                      national.df$interval==cc &
                      national.df$type==ccc &
                      national.df$year==year.i,"months_per_person"] <- exp.total.i/pop.total.i
        
      }
      r.whip <- rast(whip.grid.5y[[c]][[cc]][[ccc]])
      r.permo <- rast(permo.grid.5y[[c]][[cc]][[ccc]])
      outdir <- "study1/results/rasters/1925~2024_base_1895~1949/exposure/"
      writeCDF(r.whip,paste0(outdir,"whip_5y_",ccc,"_spei",c,"_",cc,"y.nc"),overwrite = TRUE,compression = 9)
      writeCDF(r.permo,paste0(outdir,"permo_",ccc,"_spei",c,"_",cc,"y.nc"),overwrite = TRUE,compression = 9)
    }
  }
}
national.df$timescale <- rep(c("SPEI-3","SPEI-12"),each=90)
national.df$interval <- rep(rep(c("5-year","10-year","20-year"),each=30),times=2)
national.df$type <- rep(rep(c("Overall","Dry-to-wet","Wet-to-dry"),each=10),times=6)
write.csv(national.df, "study1/results/national.df.csv", row.names = FALSE)


# population exposure decomposition
dE_pop_list <- dE_whip_list <- dE_int_list <- list()
contri <- data.frame(matrix(ncol=6,nrow=18))
names(contri) <- c("timescale","interval","type","pct_pop","pct_whip","pct_int")
contri$timescale <- rep(c("03","12"),each=9)
contri$interval <- rep(rep(c("05","10","20"),each=3),times=2)
contri$type <- rep(c("all","d2w","w2d"),times=6)

files <- list.files("G:/Gridded population/GHSL/", full.names = TRUE)
r.03 <- rast("study1/results/rasters/1925~2024_base_1895~1949/count/whip_all_spei03_10y.nc")[[1]]
usa <- geodata::gadm("USA", level=1, path=".")
conus <- usa[!usa$NAME_1 %in% c("Alaska", "Hawaii"), ]
conus <- project(conus, crs(r.03))

for (c in c("03","12")) {
  print(paste0("SPEI timescale: ",c))
  for (cc in c("05","10","20")) {
    print(paste0("Recurrence interval: ",cc," years"))
    for (ccc in c("all","d2w","w2d")) {
      print(paste0("Whiplash type: ",ccc))
      permo.grid <- rast(paste0("study1/results/rasters/1925~2024_base_1895~1949/exposure/permo_", ccc, "_spei", c, "_", cc, "y.nc"))
      E_1975 <- as.numeric(global(permo.grid[[1]],  "sum", na.rm=TRUE))
      E_2020 <- as.numeric(global(permo.grid[[10]], "sum", na.rm=TRUE))
      delta_E <- E_2020 - E_1975
      
      for (i in 1:9) {
        print(i)
        pop.t  <- rast(files[i])
        pop.t1 <- rast(files[i+1])
        pop.t.proj <- mask(project(pop.t, r.03, method = "sum"), conus)
        pop.t1.proj <- mask(project(pop.t1, r.03, method = "sum"), conus)
        
        exp.t  <- permo.grid[[i]]
        exp.t1 <- permo.grid[[i+1]]
        
        # --- recover M (whiplash months per year) ---
        mo.t  <- exp.t / pop.t.proj
        mo.t1 <- exp.t1 / pop.t1.proj
        
        # avoid division issues
        mo.t[!is.finite(mo.t)]   <- 0
        mo.t1[!is.finite(mo.t1)] <- 0
        
        # --- decomposition ---
        dE_pop  <- (pop.t1.proj - pop.t.proj) * mo.t
        dE_whip <- pop.t.proj * (mo.t1 - mo.t)
        dE_int  <- (pop.t1.proj - pop.t.proj) * (mo.t1 - mo.t)
        
        dE_pop_list[[c]][[cc]][[ccc]][[i]]  <- dE_pop
        dE_whip_list[[c]][[cc]][[ccc]][[i]] <- dE_whip
        dE_int_list[[c]][[cc]][[ccc]][[i]]  <- dE_int
      }
      
      ts_pop <- sapply(dE_pop_list[[c]][[cc]][[ccc]], function(x) {
        global(x, "sum", na.rm = TRUE)[1,1]
      })
      ts_whip <- sapply(dE_whip_list[[c]][[cc]][[ccc]], function(x) {
        global(x, "sum", na.rm = TRUE)[1,1]
      })
      ts_int <- sapply(dE_int_list[[c]][[cc]][[ccc]], function(x) {
        global(x, "sum", na.rm = TRUE)[1,1]
      })
      
      C_pop <- sum(ts_pop)
      C_whip <- sum(ts_whip)
      C_int <- sum(ts_int)

      # should be ~ equal to delta_E:
      # C_pop + C_whip + C_int
      
      # % contributions
      contri[contri$timescale==c &
               contri$interval==cc &
               contri$type==ccc, "pct_pop"] <- C_pop/delta_E*100
      contri[contri$timescale==c &
               contri$interval==cc &
               contri$type==ccc, "pct_whip"] <- C_whip/delta_E*100
      contri[contri$timescale==c &
               contri$interval==cc &
               contri$type==ccc, "pct_int"] <- C_int/delta_E*100
    }
  }
} 
contri$timescale <- rep(c("SPEI-3","SPEI-12"),each=9)
contri$interval <- rep(rep(c("5-year","10-year","20-year"),each=3),times=2)
contri$type <- rep(c("Overall","Dry-to-wet","Wet-to-dry"),times=6)
write.csv(contri, "study1/results/contributions.csv", row.names = FALSE)


# detect the hotspots for whiplash & population trend
library(spdep);library(spatialreg)
dir <- "study1/results/rasters/1925~2024_base_1895~1949/trend_50y/"
outdir <- "study1/results/rasters/1925~2024_base_1895~1949/lisa/"
slope.whip <- rast(paste0(dir, "trend_all_spei03_05y.tif"))
slope.pop <- rast(paste0(dir, "trend_pop.tif"))
slope.pop <- project(slope.pop, slope.whip, method = "sum")

moran <- data.frame(matrix(nrow=18,ncol=9))
names(moran) <- c("timescale","interval","type","global_whip","p_whip","global_pop","p_pop","global_bv","p_bv")
moran$timescale <- rep(c("03","12"),each=9)
moran$interval <- rep(rep(c("05","10","20"),each=3),times=2)
moran$type <- rep(c("all","d2w","w2d"),times=6)

lisa <- list()
nsim <- 1000
for (c in c("03","12")) {
  print(paste0("SPEI timescale: ",c))
  for (cc in c("05","10","20")) {
    print(paste0("Recurrence interval: ",cc," years"))
    for (ccc in c("all","d2w","w2d")) {
      print(paste0("Whiplash type: ",ccc))
      slope.whip.ccc <- rast(paste0(dir, "trend_",ccc,"_spei",c,"_",cc,"y.tif"))
      slope.all <- c(slope.whip.ccc, slope.pop)

      names(slope.all) <- c("whiplash", "population")
      slope.df <- as.data.frame(slope.all, xy = TRUE, na.rm = TRUE)
      coords <- as.matrix(slope.df[, c("x", "y")])
      # Create spatial neighbors (k-nearest; stable for raster grids)
      knn <- knearneigh(coords, k = 8) # at ~4 km resolution across US, k=4–8 is stable and avoids edge problems.
      nb  <- knn2nb(knn)
      # Spatial weights (row-standardized)
      lw <- nb2listw(nb, style = "W")
      
      ## Global Moran’s I (Each Raster)
      moran[moran$timescale==c & moran$interval==cc & moran$type==ccc, "global_whip"] <- moran.test(slope.df$whiplash, lw)$estimate[1]
      moran[moran$timescale==c & moran$interval==cc & moran$type==ccc, "p_whip"] <- moran.test(slope.df$whiplash, lw)$p.value
      moran[moran$timescale==c & moran$interval==cc & moran$type==ccc, "global_pop"] <- moran.test(slope.df$population, lw)$estimate[1]
      moran[moran$timescale==c & moran$interval==cc & moran$type==ccc, "p_pop"] <- moran.test(slope.df$population, lw)$p.value
      
      ## Global Bivariate Moran’s I
      x <- scale(slope.df$whiplash)[,1]
      y <- scale(slope.df$population)[,1]
      Wy <- lag.listw(lw, y)
      I_bv <- sum(x * Wy) / sum(x^2)
      
      perm_I <- replicate(nsim, {
        y_perm <- sample(y)
        Wy_perm <- lag.listw(lw, y_perm)
        sum(x * Wy_perm) / sum(x^2)
      })
      r <- sum(abs(perm_I) >= abs(I_bv))
      p_value <- (r + 1) / (nsim + 1)
      
      moran[moran$timescale==c & moran$interval==cc & moran$type==ccc, "global_bv"] <- I_bv
      moran[moran$timescale==c & moran$interval==cc & moran$type==ccc, "p_bv"] <- p_value
      
      ## Local Bivariate Moran’s I (Hotspots)
      local_I <- x * Wy
      local_perm <- matrix(NA, nrow = nrow(slope.df), ncol = nsim)
      
      for (i in 1:nsim) {
        y_perm <- sample(y)
        Wy_perm <- lag.listw(lw, y_perm)
        local_perm[, i] <- x * Wy_perm
      }
      r_local <- rowSums(abs(local_perm) >= abs(local_I))
      p_local <- (r_local + 1) / (nsim + 1)
      
      mean_x <- mean(x)
      mean_y <- mean(y)
      
      cluster <- rep("Not significant", length(x))
      sig <- p_local < 0.05
      cluster[sig & x > 0 & Wy > 0] <- "High-High"
      cluster[sig & x < 0 & Wy < 0] <- "Low-Low"
      cluster[sig & x > 0 & Wy < 0] <- "High-Low"
      cluster[sig & x < 0 & Wy > 0] <- "Low-High"
      cluster[!sig] <- NA
      
      slope.df$cluster <- cluster
      
      lisa[[c]][[cc]][[ccc]] <- rast(slope.whip)
      values(lisa[[c]][[cc]][[ccc]]) <- NA
      
      # match order
      cells <- cellFromXY(slope.whip.ccc, coords)
      cluster_factor <- factor(
        slope.df$cluster,
        levels = c("High-High", "Low-Low", "High-Low", "Low-High")
      )
      
      lisa[[c]][[cc]][[ccc]][cells] <- cluster_factor
      levels(lisa[[c]][[cc]][[ccc]]) <- data.frame(
        id = 1:4, 
        label = c("High-High", "Low-Low", "High-Low", "Low-High")
      )
      writeCDF(lisa[[c]][[cc]][[ccc]], paste0(outdir,"lisa_",ccc,"_spei",c,"_",cc,"y.nc"),overwrite = TRUE,compression = 9)
    }
  }
}

