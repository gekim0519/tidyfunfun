data preparation
================
Gaeun Kim
2/27/2019

handwrting data

``` r
# cursive handwriting coordinates
#load data
load("../data/handwrit.rda")

handwrit_time <- handwritTime %>% as.tibble() %>%
  janitor::clean_names()
rm(handwritTime)

handwrit <- handwrit %>%
  as.tibble() %>%
  janitor::clean_names()
```

``` r
plot(handwrit$rep01_x, handwrit$rep01_y, type="l")
```

![](preparation_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
plot(handwrit$rep02_x, handwrit$rep02_y, type="l")
```

![](preparation_files/figure-markdown_github/unnamed-chunk-2-2.png)

Making handwriting data to tfd format!

``` r
hw_x = handwrit %>% 
  as.matrix() %>%
  .[,1:20] %>%
  t()
hw_y = handwrit %>% 
  as.matrix() %>%
  .[,21:40] %>%
  t()
  
handw_tfd = data.frame(id = 1:20) %>%
  as.tbl %>% 
  mutate(x = tfd(hw_x),
         y = tfd(hw_y))

handw_tfd %>% 
  filter(id == 1) %>% 
  ggplot(aes(y = y)) + 
  geom_spaghetti()
```

![](preparation_files/figure-markdown_github/unnamed-chunk-3-1.png)

we can't combine x and y coordinates in tdf format :( (x and y are not all unique values)

CanadianWeather
---------------

Consists of the following components:

dailyAv a three dimensional array c(365, 35, 3) summarizing data collected at 35 different weather stations in Canada on the following:

\[,,1\] = \[,, 'Temperature.C'\]: average daily temperature for each day of the year

\[,,2\] = \[,, 'Precipitation.mm'\]: average daily rainfall for each day of the year rounded to 0.1 mm.

\[,,3\] = \[,, 'log10precip'\]: base 10 logarithm of Precipitation.mm after first replacing 27 zeros by 0.05 mm (Ramsay and Silverman 2006, p. 248).

place Names of the 35 different weather stations in Canada whose data are summarized in 'dailyAv'. These names vary between 6 and 11 characters in length. By contrast, daily\[\["place"\]\] which are all 11 characters, with names having fewer characters being extended with trailing blanks.

province names of the Canadian province containing each place

coordinates a numeric matrix giving 'N.latitude' and 'W.longitude' for each place.

region Which of 4 climate zones contain each place: Atlantic, Pacific, Continental, Arctic.

monthlyTemp A matrix of dimensions (12, 35) giving the average temperature in degrees celcius for each month of the year.

monthlyPrecip A matrix of dimensions (12, 35) giving the average daily precipitation in milimeters for each month of the year.

geogindex Order the weather stations from East to West to North

``` r
load("../data/CanadianWeather.rda")

cwtemp = CanadianWeather[[1]] %>% as.tibble() %>%
  t() %>% as.matrix() %>%
  .[1:35,]

place = c(CanadianWeather$place)
temp_df = data.frame(place = CanadianWeather$place)
temp_df$temp = cwtemp

temp_tfd = with(temp_df, 
  data.frame(place = place)) %>% as.tbl %>% 
  mutate(temp = tfd(temp_df$temp))
```

    ## Column names not suitable as arg. Using 1:ncol(data).

``` r
# tfd(temp_df$temp)) (1, -4) to jan1, -4? 
```

Made just the temp data to tdf format!

``` r
temp_tfd %>% 
  filter(place == "Halifax") %>% 
  ggplot(aes(y = temp)) + 
  geom_spaghetti() +
  ggtitle(label = "Average daily temperature for each day of the year in Halifax")
```

![](preparation_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
temp_tfd %>%  
  filter(place %in% c("Quebec", "Vancouver")) %>% 
  ggplot(aes(y = temp, color = place)) + 
  geom_spaghetti(alpha = .3) + 
  geom_meatballs(aes(alpha = .5)) +
  facet_grid(~ place)
```

![](preparation_files/figure-markdown_github/unnamed-chunk-5-2.png)

saving datasets to r format

``` r
save(temp_tfd, file = "../data/temp_tdf.RData")
save(handw_tfd, file = "../data/handw_tdf.RData")
```

Example at Tidyfun page

``` r
# tidyfun example
DTI = refund::DTI

dti = with(refund::DTI, 
  data.frame(id = ID, sex = sex, 
    case = factor(ifelse(case, "MS", "control")))) %>% as.tbl %>% 
        mutate(cca = tfd(DTI$cca, seq(0,1, l = 93), signif = 2) %>%
                     tfd(arg = seq(0,1,l = 93)),
               rcst = tfd(DTI$rcst, seq(0, 1, l = 55), signif = 3))

dti %>%
  filter(id == 2031)
```

    ## # A tibble: 2 x 5
    ##      id sex    case  cca                         rcst                      
    ##   <dbl> <fct>  <fct> <S3: tfd_reg>               <S3: tfd_irreg>           
    ## 1  2031 female MS    [1]: (0.000,0.51);(0.011,0… [1]: (0.000,0.51);(0.019,…
    ## 2  2031 female MS    [2]: (0.000,0.55);(0.011,0… [2]: (0.000,0.50);(0.019,…

Refund fpca
-----------

refund cd4 example

``` r
data(cd4)
# CD4 cell counts for 366 subjects between months -18 and 42 since seroconversion. Each subject's observations are contained in a single row.
# subject * weeks
Fit.MM = fpca.sc(cd4, var = TRUE, simul = TRUE)

Fit.mu = data.frame(mu = Fit.MM$mu,
                    d = as.numeric(colnames(cd4)))
Fit.basis = data.frame(phi = Fit.MM$efunctions, # d × npc matrix of estimated eigenfunctions of the functional covariance, i.e., the FPC basis functions.
                       d = as.numeric(colnames(cd4)))

## for one subject, examine curve estimate, pointwise and simultaneous itervals
EX = 1
EX.MM = data.frame(fitted = Fit.MM$Yhat[EX,],
           ptwise.UB = Fit.MM$Yhat[EX,] + 1.96 * sqrt(Fit.MM$diag.var[EX,]),
           ptwise.LB = Fit.MM$Yhat[EX,] - 1.96 * sqrt(Fit.MM$diag.var[EX,]),
           simul.UB = Fit.MM$Yhat[EX,] + Fit.MM$crit.val[EX] * sqrt(Fit.MM$diag.var[EX,]),
           simul.LB = Fit.MM$Yhat[EX,] - Fit.MM$crit.val[EX] * sqrt(Fit.MM$diag.var[EX,]),
           d = as.numeric(colnames(cd4)))

## plot data for one subject, with curve and interval estimates
EX.MM.m = melt(EX.MM, id = 'd')
ggplot(EX.MM.m, aes(x = d, y = value, group = variable, color = variable, linetype = variable)) +
  geom_path() +
  scale_linetype_manual(values = c(fitted = 1, ptwise.UB = 2,
                        ptwise.LB = 2, simul.UB = 3, simul.LB = 3)) +
  scale_color_manual(values = c(fitted = 1, ptwise.UB = 2,
                     ptwise.LB = 2, simul.UB = 3, simul.LB = 3)) +
  labs(x = 'Months since seroconversion', y = 'Total CD4 Cell Count')
```

![](preparation_files/figure-markdown_github/unnamed-chunk-8-1.png)

``` r
## plot estimated mean function
ggplot(Fit.mu, aes(x = d, y = mu)) + geom_path() +
  labs(x = 'Months since seroconversion', y = 'Total CD4 Cell Count')
```

![](preparation_files/figure-markdown_github/unnamed-chunk-8-2.png)

``` r
## plot the first two estimated basis functions
Fit.basis.m = melt(Fit.basis, id = 'd')
ggplot(subset(Fit.basis.m, variable %in% c('phi.1', 'phi.2')), aes(x = d,
y = value, group = variable, color = variable)) + geom_path()
```

![](preparation_files/figure-markdown_github/unnamed-chunk-8-3.png)

``` r
## input a dataframe instead of a matrix
nid <- 20
nobs <- sample(10:20, nid, rep=TRUE)
ydata <- data.frame(
    .id = rep(1:nid, nobs),
    .index = round(runif(sum(nobs), 0, 1), 3)) # long format
ydata$.value <- unlist(tapply(ydata$.index,
                              ydata$.id,
                              function(x)
                                  runif(1, -.5, .5) +
                                  dbeta(x, runif(1, 6, 8), runif(1, 3, 5))
                              )
                       )

Fit.MM = fpca.sc(ydata=ydata, var = TRUE, simul = FALSE)
```

``` r
fit.cw = fpca.sc(cwtemp, var = TRUE, simul = TRUE)

fit.mu = data.frame(mu = fit.cw$mu,
                    days = 1:365)
fit.basis = data.frame(phi = fit.cw$efunctions, #the FPC basis functions.
                       days = 1:365)

temp_tfd %>%  
  ggplot(aes(y = temp, color = place)) + 
  geom_spaghetti(alpha = .3) + xlab("days") + theme(legend.position="none")
```

![](preparation_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
## plot the first two estimated basis functions
fit.basis.m = melt(fit.basis, id = 'days')
ggplot(subset(fit.basis.m, variable %in% c('phi.1', 'phi.2')), aes(x = days,
y = value, group = variable, color = variable)) + geom_path() + ggtitle("first two estimated basis functions")
```

![](preparation_files/figure-markdown_github/unnamed-chunk-9-2.png)

``` r
## plot estimated mean function
ggplot(fit.mu, aes(x = days, y = mu)) + geom_path() +
  labs(x = 'Day', y = 'Average daily temperature')
```

![](preparation_files/figure-markdown_github/unnamed-chunk-9-3.png)

``` r
# estimated eigenvalues of the covariance operator, i.e., variances of FPC scores.
fit.cw$evalues
```

    ## [1] 41.4545221  3.9435175  0.8774614

``` r
fit.cw$evalues*100/sum(fit.cw$evalues)
```

    ## [1] 89.582006  8.521826  1.896168

``` r
# rcst A 382 x 55 matrix of fractional anisotropy tract profiles from the right corticospinal tract;
rcst <- dti$rcst %>% 
  as.data.frame() %>%
  spread(key = arg, value = value) %>%
  select(-id) %>%
  as.matrix() 

fit.rcst = fpca.sc(rcst, var = TRUE, simul = TRUE)

# df
rcst_df <- dti$rcst %>% 
  as.data.frame() %>%
  select(.id = id, .index = arg, .value = value) %>%
  mutate(.id = as.integer(.id))
```

``` r
#testing diff tfd
cca.fit = dti$cca %>% 
  as.data.frame() %>%
  spread(key = arg, value = value) %>%
  select(-id) %>%
  as.matrix() %>%
  fpca.sc()
  

temp.fit = temp_tfd$temp %>% 
  as.data.frame() %>%
  spread(key = arg, value = value) %>%
  select(-id) %>%
  as.matrix() %>%
  fpca.sc(Y =.)
```