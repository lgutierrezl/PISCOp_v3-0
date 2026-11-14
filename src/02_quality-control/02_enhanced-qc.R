rm(list = ls())

suppressMessages(library(xts))
source("R/utils/pipe.R")
source("R/02_quality-control/rounding-and-precision.R")
source("R/02_quality-control/weekly-cycles.R")
source("R/02_quality-control/small-gaps.R")
source("R/02_quality-control/truncation.R")

ncores <- 6

raw_data <- readRDS(
  file.path(
    "data",
    "processed",
    "point",
    "datos_PE_1981-2025_v3.0_02.RDS"
  )
)

raw_data$data[raw_data$data < 0] <- NA

raw_data_xyz <- raw_data$xyz[raw_data$xyz$SIZE == 1, ]
raw_data_data <- raw_data$data[,
  raw_data$xyz[raw_data$xyz$SIZE == 1, ]$ID
]


enhanced_qc_res <-  
  raw_data_data %>% 
  parallel::mclapply(function(z) {

    data.frame(truncation = get_level_trunc(xts_obj = z),
               small_gaps = get_levels_sg(xts_obj = z),
               weekly_cycles = get_levels_wc(xts_obj = z),
               precision_rounding = get_levels_rpp(xts_obj = z),
               lenght_decimals_patterns = get_dec_patterns_mode(xts_obj = z),
               wet_day = wd_value(xts_obj = z))

  }, mc.cores = ncores) %>%
  do.call(rbind, .)

# adding more parameters due to that SA
# has a diverse range of wet days and 
# problems with the precision 

enhanced_qc_res <- transform(
  enhanced_qc_res,
  truncation2 = ifelse(wet_day < 15, 0, truncation),
  small_gaps2 = ifelse(wet_day < 15, 0, small_gaps),
  weekly_cycles2 = ifelse(wet_day < 5, 0, weekly_cycles),
  precision_rounding2 = ifelse(wet_day < 15, 0, precision_rounding)
)

enhanced_qc_res <- transform(
  enhanced_qc_res,
  small_gaps3 = ifelse(lenght_decimals_patterns > 25, 0, small_gaps2)
)

enhanced_qc_res <- enhanced_qc_res[, c("truncation2", "small_gaps3", "weekly_cycles2", "precision_rounding2")]
colnames(enhanced_qc_res) <- c("truncation", "small_gaps", "weekly_cycles", "precision_rounding")

all(raw_data_xyz$ID == rownames(enhanced_qc_res))
raw_data_xyz <- cbind(raw_data_xyz, enhanced_qc_res)

# best stations
## truncation == 2 | small_gaps == 2 are flagged for the next steps
raw_data_xyz <- transform(
  raw_data_xyz,
  EQC = ifelse(truncation == 2 | small_gaps == 2, 0, 1)
)

saveRDS(
  raw_data_xyz,
  file.path(
    "output",
    "02_quality-control",
    "eqc-01-summary.RDS"
  )
)

# add qc visual
raw_data_data['1985/1990','PE105040'] <- NA#
raw_data_data['1981/1998','PE107002'] <- NA#
raw_data_data['1981/1991-06','PE106061'] <- NA#
raw_data_data[c('1988-02/1988-06','1988-11/1989-05','1990/1993','1994-01/1994-05'),'PE106077'] <- NA#
raw_data_data['1988/1992','PE108045'] <- NA
raw_data_data['1984/1997','PE109014'] <- NA
raw_data_data[c('1996/1999', '2003/2011'),'PE111023'] <- NA
raw_data_data['2007/2013','PE104065'] <- NA
raw_data_data['/1985','PE111025'] <- NA
raw_data_data['1994/2010','PE110017'] <- NA
raw_data_data['/1984','CO21010100'] <- NA
raw_data_data['1997','CO44015070'] <- NA
raw_data_data['2013/2015','CO44040020'] <- NA
raw_data_data[c('2007-09/2009', '2013/2015'),'CO44050010'] <- NA
raw_data_data[c('1986-07/1998-03','2021/2025'),'CO44187020'] <- NA
raw_data_data['1997-03/1999-07','CO47010020'] <- NA
raw_data_data['1982-08/1986-03','CO47010050'] <- NA
raw_data_data[c('1986/1987','1995-11/1996-11','2000-12/2002-03','2005-04/2005-10','2011-07/2013-06'),'CO47087010'] <- NA
raw_data_data[c('1995-06/1996-03','2013-09/2014-10'),'CO47107010'] <- NA
raw_data_data['1993-10/1999-02','CO52070010'] <- NA
raw_data_data[c('1994-03/1994-09','1995-02/1997'),'BR000964001'] <- NA
raw_data_data['2000/2004','EC0150'] <- NA
raw_data_data['2010/2012','EC0414'] <- NA
raw_data_data['2004/2006-08','EC0421'] <- NA
raw_data_data['1992/1993','EC0429'] <- NA
#

all(colnames(raw_data_data) == raw_data_xyz$ID)
raw_data <- list(
  data = raw_data_data,
  xyz = raw_data_xyz[,
   c("ID", "NAM", "LON", "LAT", "ALT", "COUNTRY", "SOURCE","ECOREGIONS", "SIZE", "EQC")
  ]
)                  

# evl eqc
raw_data_xyz <- raw_data$xyz
raw_data_dat <- raw_data$data

raw_data_dat['2022/2024','PE103044'] <- NA
raw_data_dat[c('/1995','2005'),'PE105095'] <- NA
raw_data_dat['1983/1993','PE105023'] <- NA
raw_data_dat['/1996','PE105042'] <- NA
raw_data_dat['1992/1994','PE100132'] <- NA
raw_data_dat[c('1981/1986', '1993'),'PE100092'] <- NA
raw_data_dat['1981/1986','PE106019'] <- NA
raw_data_dat['1981','PE107012'] <- NA
raw_data_dat['1995/1997','PE107068'] <- NA
raw_data_dat['1990/1994','PE107054'] <- NA
raw_data_dat['1984-07/1985-05','PE107092'] <- NA
raw_data_dat[c('1981/1994','2014-07/2015-05'),'PE108044'] <- NA
raw_data_dat['1998-07/2003-05','PE109018'] <- NA
raw_data_dat['1984-07/1986-05','PE109045'] <- NA
raw_data_dat['2016/2025','PE110041'] <- NA
raw_data_dat[c('1987-06/1992','2006/2009'),'PE102006'] <- NA
raw_data_dat[c('/2002'),'PE105092'] <- NA
raw_data_dat['1998/2000','PE100105'] <- NA
raw_data_dat['1987-07/1989-05','PE108028'] <- NA
raw_data_dat['/1993','PE109022'] <- NA
raw_data_dat['/1992','PE109049'] <- NA
raw_data_dat['/1998','PE111018'] <- NA
raw_data_dat['/1984','PE111028'] <- NA
raw_data_dat['1997-01/1997-03','PE110063'] <- NA
raw_data_dat['1990','PE111070'] <- NA
raw_data_dat['/1990','PE111087'] <- NA
raw_data_dat['1987/1993','PE111026'] <- NA
raw_data_dat['1987/1993','PE111085'] <- NA
raw_data_dat[c('/1992','2015-06/2020'),'PE112130'] <- NA
raw_data_dat[c('1981/1983', '1986/1992-05'),'PE112127'] <- NA
raw_data_dat['/1993-05','PE112012'] <- NA
raw_data_dat['1983/1992-05','PE112057'] <- NA
raw_data_dat[c('1982-07/1983-05', '1986-10/1993-05', '2004-07/2005-05'),'PE113086'] <- NA
raw_data_dat['1988-10/1989-03','PE114023'] <- NA
raw_data_dat['1981/1983','PE115087'] <- NA
raw_data_dat['1984-06/1987-03','CO47010030'] <- NA
raw_data_dat['1997/1999-01','EC0127'] <- NA
raw_data_dat[c('1991','1992-12/1993-03'),'EC0479'] <- NA
raw_data_dat['1981/1995','EC0485'] <- NA
raw_data_dat['1988/1993','PE115021'] <- NA
raw_data_dat[c('1994-01-31','2002-02-05','2020-02-14/2020-02-15', '2020-02-22'),'PE111060'] <- NA

# new qc
raw_data_dat[c('1984/1991','1995-07/1997-01'),'PE102009'] <- NA
raw_data_dat[c('1988/1996-06'),'PE106083'] <- NA
raw_data_dat['1981/1990','PE106088'] <- NA
raw_data_dat['1981/1992','PE106112'] <- NA
raw_data_dat['1981/1990','PE100104'] <- NA
raw_data_dat['1981/1985','PE107018'] <- NA
raw_data_dat['1981/1985','PE107069'] <- NA
raw_data_dat['1981/1984','PE106092'] <- NA
raw_data_dat['1981/2008-07','PE111084'] <- NA
raw_data_dat['1981/1992','PE111095'] <- NA
raw_data_dat['1981/1994','PE112055'] <- NA
raw_data_dat['1981/1998','PE113122'] <- NA
raw_data_dat[c('1981/1984','2021/2025'),'PE113044'] <- NA
raw_data_dat[c('1981/1996-06','1999-01/1999-03'),'PE113101'] <- NA
raw_data_dat[c('1981/1992'),'PE113035'] <- NA
raw_data_dat['1981/1996','PE107068'] <- NA
raw_data_dat['1993','PE108045'] <- NA
raw_data_dat['1987/1995','PE108021'] <- NA
raw_data_dat['1981/1985','PE108047'] <- NA
raw_data_dat['1981/1992','PE108017'] <- NA
raw_data_dat['2012-02-06/2012-02-07','PE109014'] <- NA
raw_data_dat['2005/2014','PE110018'] <- NA
raw_data_dat['1993-02/1996-03','PE110007'] <- NA
raw_data_dat['1981/1984','PE103024'] <- NA
raw_data_dat['1981/1984','PE106014'] <- NA
raw_data_dat[c('1998-01-07', '1998-03-19', '1998-08-23'),'PE110015'] <- NA
raw_data_dat['2018-07/2022','PE110037'] <- NA
raw_data_dat['1985/1996-05','PE111097'] <- NA
raw_data_dat[c('2011-02','2018-06/2018-07','2019-06/2019-08'),'PE112171'] <- NA
raw_data_dat['2020-06/2025','PE113067'] <- NA
raw_data_dat[c('1988/1991','1995/1997'),'PE113081'] <- NA
raw_data_dat['1981/1984','PE113034'] <- NA
raw_data_dat['1985/1998','PE113116'] <- NA
raw_data_dat['1981/1985','PE113025'] <- NA
raw_data_dat['2014/2019','PE113038'] <- NA
raw_data_dat['1981/1985','PE114033'] <- NA
raw_data_dat['1981/1994','PE114072'] <- NA
raw_data_dat['1981/1985','PE114034'] <- NA
raw_data_dat['1981/1993-06','PE115101'] <- NA
raw_data_dat[c('1981/1983','1989/1994'),'PE115088'] <- NA
raw_data_dat['1981/1986','PE116052'] <- NA
raw_data_dat['1985-11/1986-04','PE117048'] <- NA
raw_data_dat['1985-11/1986-04','PE112028'] <- NA
raw_data_dat['1981/1985','PE112037'] <- NA
raw_data_dat['1981/1990','PE113119'] <- NA
raw_data_dat['2023/2025','PE116022'] <- NA
raw_data_dat[c('2007/2009','2019/2021'),'PE104071'] <- NA
raw_data_dat[c('1985','1993/1998'),'PE103030'] <- NA
raw_data_dat[c('1984-08/1985-01'),'CO47030010'] <- NA
raw_data_dat[c('2020-12/2021-09'),'CH1310019'] <- NA
raw_data_dat[c('1991/1995','2005/2007'),'BR82610'] <- NA
raw_data_dat[c('1982/1989-11','2008','2016/2018-11','2023-10/'),'COM00080398'] <- NA
raw_data_dat[c('2021/2022'),'PEM00084515'] <- NA
raw_data_dat[c('1981/1990','2016/'),'EC0122'] <- NA
raw_data_dat[c('2012'),'EC0490'] <- NA
raw_data_dat[c('2009/2011-06'),'CH1200001'] <- NA
raw_data_dat[c('1992/1996-04'),'PE103046'] <- NA
raw_data_dat[c('1986/1989'),'PE105090'] <- NA
raw_data_dat[c('1981/1985'),'PE106016'] <- NA
raw_data_dat[c('1984-11/1984-12-15','1985-01-15/1985-03-14','2020-08'),'PE106017'] <- NA
raw_data_dat[c('1981/1992'),'PE110019'] <- NA
raw_data_dat[c('2003/2004'),'PE110027'] <- NA
raw_data_dat[c('/1997'),'PE111022'] <- NA
raw_data_dat[c('2000/2003'),'PE100028'] <- NA

raw_data_dat[c('1988/1989','2016/2020'),'CO21010100'] <- NA
raw_data_dat[c('2002/2004'),'CO21010190'] <- NA
raw_data_dat[c('2014-03-02','2014-03-08','2014-04-25','2016/2019'),'CO21010200'] <- NA
raw_data_dat[c('1992'),'CO21010250'] <- NA
raw_data_dat[c('1983/1984','2018/2019'),'CO21030060'] <- NA
raw_data_dat[c('1986/1989','2023'),'CO42050020'] <- NA
raw_data_dat[c('1996/2004'),'CO42070010'] <- NA
raw_data_dat[c('2011-07/2015'),'CO42070020'] <- NA
raw_data_dat[c('1985/1987','1994','2020/2022'),'CO42080020'] <- NA
raw_data_dat[c('1994-11/1995-03'),'CO44055010'] <- NA
raw_data_dat[c('2001-10'),'CO44100010'] <- NA
raw_data_dat[c('1989-05/1989-06'),'CO44127020'] <- NA
raw_data_dat[c('1993-08/1994-08'),'CO47010090'] <- NA
raw_data_dat[c('1991/1994'),'CO47010170'] <- NA
raw_data_dat[c('1989','1992'),'CO47017030'] <- NA
raw_data_dat[c('1994-08/1995'),'CO47030030'] <- NA
raw_data_dat[c('1992-08/1994-10','2002-04/2003'),'CO47040040'] <- NA
raw_data_dat[c('2014-11/2015-06'),'CO47107030'] <- NA
raw_data_dat[c('1992-07/1993-08'),'CO52050060'] <- NA
raw_data_dat[c('2003/2015'),'EC0101'] <- NA
raw_data_dat[c('2004-12/2005'),'EC0167'] <- NA
raw_data_dat[c('1993-07/1994-06'),'CO47010060'] <- NA

raw_data_xyz[raw_data_xyz$ID %in% c('PE106083','PE106112','PE100104','PE106092','PE105040',"PE111095",'PE111027','PE112055',
                                    'PE113122','PE113100','PE113044',"PE113101","PE114042",'PE100007','PE100072','PE113035',
                                    'PE115047','PE117041','BR82610','CO44115020','CO47010180','BO85201','CH1020016','CH1020017',
                                    'COM00080398','PEM00084515','CO44015070','CO52010140','EC0122','EC0351','EC0794','CH1200001',
                                    "PE102009","PE103016","PE103012","PE103044","PE103011","PE103015","PE103046","PE103016",
                                    "PE107054","PE103049","PE105095","PE105078","PE105055","PE105090","PE105020","PE105023",
                                    "PE105042","PE100132","PE106115","PE106053","PE100092","PE106019","PE107091",
                                    "PE107012","PE107013","PE107035","PE107038","PE107069","PE107068","PE107054",
                                    "PE107092","PE108044","PE100137","PE108064","PE108021","PE108047","PE109018",
                                    "PE108025","PE109045","PE109019","PE109016",#"id_109021",
                                    "PE109038","PE100020",
                                    "PE110018","PE110040","PE110059","PE110041","PE102006","PE104097","PE100131",
                                    "PE105092","PE100105","PE106054","PE100113","PE107010","PE108033","PE108028",
                                    "PE109022","PE109049","PE111583","PE111018","PE111028","PE110020","PE110057",
                                    "PE110021","PE110063","PE111070","PE111081","PE111087","PE111067","PE111026",
                                    "PE111085","PE111088","PE111114","PE111091","PE111093","PE112130","PE112126",
                                    "PE112080","PE112127","PE112012","PE112060","PE112057","PE113067","PE113034",
                                    "PE100101","PE113086","PE113041","PE113025","PE113022","PE114023","PE114035",
                                    "PE115082","PE115084","PE115023","PE115085","PE115089","PE115052","PE116046",
                                    "PE111089","PE111144","PE111061","PE111175","PE112051","PE100109","PE100055",
                                    "PE114050","PE115076","PE115043","PE115087","PE116022","CO47010030","CO51020010",
                                    "CO52040040","CO52080010","PE153312","PE152134", "PE140406", "BO2090", "EC0127",
                                    "EC0143","EC0147", "EC0479","EC0485"),'EQC'] <- 1

saveRDS(
  raw_data <- list(
    data = raw_data_dat,
    xyz = raw_data_xyz),
  file.path(
    "data",
    "processed",
    "point",
    "datos_PE_1981-2025_v3.0_03.RDS"
  )
)

