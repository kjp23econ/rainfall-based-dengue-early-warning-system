# ============================================================
# Dengue Early Warning System — Full Reproducible Script
# Author: Keanu John Pelitro
# ============================================================

# ---------------------------
# 0. Housekeeping
# ---------------------------
rm(list = ls())
set.seed(123)

# ---------------------------
# 1. Load required libraries
# ---------------------------
packages <- c(
  "here","readxl","dplyr","tidyr","lubridate","ISOweek","zoo",
  "ggplot2","scales","cowplot","gridExtra","grid",
  "RColorBrewer","dlnm","splines","reshape2","viridis"
)

need <- setdiff(packages, rownames(installed.packages()))
if(length(need) > 0) install.packages(need, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

# ---------------------------
# 2. Paths and folders
# ---------------------------
DATA_PATH <- here::here("Analysis 2.xlsx")
FIG_DIR <- here::here("figures")
if(!dir.exists(FIG_DIR)) dir.create(FIG_DIR)

# ---------------------------
# 3. Load data
# ---------------------------
df <- read_excel(DATA_PATH, sheet = "ExtRF_Lancet") %>%
  mutate(
    YR = as.integer(YR),
    WN = as.integer(WN),
    DC = as.numeric(DC),
    RF = as.numeric(RF),
    ISOweek = sprintf("%d-W%02d", YR, WN),
    Date = ISOweek2date(paste0(ISOweek, "-1"))
  ) %>%
  arrange(Date)

df_raw <- read_excel(DATA_PATH, sheet = "ExtRF_Lancet") # Must include PedAdmit, PedDead

# ---------------------------
# 4. Global parameters
# ---------------------------
wet_weeks <- 23:39
pandemic_years <- c(2020, 2021)
target_years <- 2019:2025
cal_years <- c(2019, 2022, 2023)
eval_years <- c(2024, 2025)

# ============================================================
# SHARED FUNCTIONS
# ============================================================
calc_threshold <- function(data){
  mean(data$DC, na.rm = TRUE) + 2 * sd(data$DC, na.rm = TRUE)
}

rolling_threshold_from_map <- function(df, rolling_map){
  bind_rows(lapply(names(rolling_map), function(y){
    yrs <- rolling_map[[y]]
    data.frame(
      YR = as.integer(y),
      Threshold = calc_threshold(df %>% filter(YR %in% yrs)),
      PriorYears = paste(yrs, collapse = ",")
    )
  }))
}

compute_youden_curve <- function(df_input, lag_weeks = 3, rf_grid = NULL){
  if(is.null(rf_grid)){
    RF_lower <- median(df_input$RF[df_input$RF > 0], na.rm = TRUE)
    RF_cap   <- quantile(df_input$RF, 0.99, na.rm = TRUE)
    rf_grid  <- seq(RF_lower, RF_cap, length.out = 200)
  }
  bind_rows(lapply(rf_grid, function(t){
    df_tmp <- df_input %>% mutate(pred_flag = RF >= t)
    metrics <- bind_rows(lapply(unique(df_input$YR), function(y){
      df_y <- df_tmp %>% filter(YR == y)
      surge_wk <- df_y$WN[which.max(df_y$DC)]
      df_y <- df_y %>% mutate(obs_flag = WN >= surge_wk - lag_weeks + 1 & WN <= surge_wk)
      TP <- sum(df_y$pred_flag & df_y$obs_flag, na.rm=TRUE)
      FP <- sum(df_y$pred_flag & !df_y$obs_flag, na.rm=TRUE)
      TN <- sum(!df_y$pred_flag & !df_y$obs_flag, na.rm=TRUE)
      FN <- sum(!df_y$pred_flag & df_y$obs_flag, na.rm=TRUE)
      data.frame(
        sensitivity = ifelse(TP+FN==0, NA, TP/(TP+FN)),
        specificity = ifelse(TN+FP==0, NA, TN/(TN+FP))
      )
    }))
    sens <- mean(metrics$sensitivity, na.rm=TRUE)
    spec <- mean(metrics$specificity, na.rm=TRUE)
    data.frame(RF_threshold = t, sensitivity = sens, specificity = spec, youden = sens + spec -1)
  }))
}

# ============================================================
# FIGURE 1 — Weekly Dengue Cases With vs Without Pandemic
# ============================================================
rolling_map_with_pandemic <- list(
  "2019"=2016:2018, "2020"=2016:2019, "2021"=2016:2020,
  "2022"=2017:2021, "2023"=2018:2022, "2024"=2019:2023, "2025"=2020:2024
)
rolling_map_no_pandemic <- list(
  "2019"=2016:2018, "2022"=2016:2019, "2023"=c(2017:2019,2022),
  "2024"=c(2017:2019,2022:2023), "2025"=c(2018:2019,2022:2024)
)

df_plot_full <- df %>% filter(YR %in% target_years) %>% mutate(WeekSeq=row_number())
df_plot_no <- df %>% filter(YR %in% target_years & !YR %in% pandemic_years) %>% mutate(WeekSeq=row_number())

df_thresh_full <- rolling_threshold_from_map(df, rolling_map_with_pandemic)
df_thresh_no <- rolling_threshold_from_map(df %>% filter(!YR %in% pandemic_years), rolling_map_no_pandemic)

plot_dengue <- function(df_data, df_thresh, title_text){
  ggplot(df_data, aes(WeekSeq, DC)) +
    geom_line(color="black") +
    geom_rect(
      data=df_data %>% filter(WN %in% wet_weeks),
      aes(xmin=WeekSeq-0.5, xmax=WeekSeq+0.5, ymin=-Inf, ymax=Inf),
      inherit.aes=FALSE, fill="lightblue1", alpha=0.2
    ) +
    geom_line(
      data=left_join(df_data, df_thresh, by="YR"),
      aes(y=Threshold),
      color="orange", linetype="dashed", linewidth=1
    ) +
    labs(title=title_text, x="Year", y="Weekly Dengue Cases") +
    theme_minimal()
}

p1 <- plot_dengue(df_plot_full, df_thresh_full, "Figure 1: Weekly Dengue Cases With Pandemic Years")
p2 <- plot_dengue(df_plot_no, df_thresh_no, "Figure 1: Weekly Dengue Cases Without Pandemic Years")
plot_grid(p1,p2,ncol=2)
ggsave(file.path(FIG_DIR,"Figure1_WeeklyDengue.png"), width=12, height=6)
message("✓ Figure 1 generated.")

# ============================================================
# FIGURE 2 — Weekly Rainfall & Dengue
# ============================================================
df_fig2 <- df %>% filter(YR %in% target_years)
scale_factor <- max(df_fig2$RF, na.rm = TRUE)/max(df_fig2$DC, na.rm = TRUE)
p_fig2 <- ggplot(df_fig2, aes(Date)) +
  geom_col(aes(y=RF), alpha=0.6) +
  geom_line(aes(y=DC*scale_factor), linewidth=0.8) +
  scale_y_continuous(
    name="Weekly Rainfall (mm)",
    sec.axis = sec_axis(~./scale_factor, name="Weekly Dengue Cases")
  ) +
  labs(title="Figure 2: Weekly Rainfall and Dengue Cases",
       subtitle="Dual-axis visualization for temporal alignment",
       x="Date") +
  theme_minimal()
print(p_fig2)
ggsave(file.path(FIG_DIR,"Figure2_Rainfall_Dengue.png"), width=12, height=6)
message("✓ Figure 2 generated.")

# ============================================================
# FIGURE 3 — Per-Year Youden Curves
# ============================================================
plot_youden_year <- function(df_year, year_label, side="right"){
  yc <- compute_youden_curve(df_year)
  opt <- yc[which.max(yc$youden),]
  hjust_val <- ifelse(side=="left",1.1,-0.1)
  ggplot(yc, aes(RF_threshold, youden)) +
    geom_line(linewidth=1) +
    geom_vline(xintercept=opt$RF_threshold, linetype="dashed") +
    annotate("text", x=opt$RF_threshold, y=max(yc$youden,na.rm=TRUE),
             label=paste0("Opt=",round(opt$RF_threshold,1)," mm"),
             hjust=hjust_val, size=3) +
    labs(title=paste("Figure 3: Youden Curve —",year_label),
         x="Rainfall Threshold (mm)", y="Youden Index") +
    theme_minimal()
}
p_y2019 <- plot_youden_year(df %>% filter(YR==2019), 2019,"left")
p_y2022 <- plot_youden_year(df %>% filter(YR==2022), 2022,"right")
p_y2023 <- plot_youden_year(df %>% filter(YR==2023), 2023,"right")
plot_grid(p_y2019,p_y2022,p_y2023,ncol=1)
ggsave(file.path(FIG_DIR,"Figure3_Youden_PerYear.png"), width=6, height=10)
message("✓ Figure 3 generated.")

# ============================================================
# FIGURE 4 — Pooled Youden Curve
# ============================================================
df_cal <- df %>% filter(YR %in% cal_years)
yc_pooled <- compute_youden_curve(df_cal)
opt_pooled <- yc_pooled[which.max(yc_pooled$youden),]
per_year_opts <- sapply(cal_years, function(y){
  yc <- compute_youden_curve(df %>% filter(YR==y))
  yc$RF_threshold[which.max(yc$youden)]
})
case_thr_range <- range(per_year_opts)
p_fig4 <- ggplot(yc_pooled, aes(RF_threshold,youden)) +
  geom_line(linewidth=1) +
  geom_vline(xintercept=case_thr_range, linetype="dashed") +
  geom_vline(xintercept=opt_pooled$RF_threshold, linewidth=1) +
  annotate("text", x=opt_pooled$RF_threshold, y=max(yc_pooled$youden,na.rm=TRUE),
           label=paste0("Pooled Opt=",round(opt_pooled$RF_threshold,1)," mm"), hjust=-0.1, size=3) +
  labs(title="Figure 4: Pooled Youden Curve",
       subtitle=paste0("Per-year optimal rainfall threshold range: ", round(case_thr_range[1],1),"-",round(case_thr_range[2],1)," mm"),
       x="Rainfall Threshold (mm)", y="Youden Index") +
  theme_minimal()
print(p_fig4)
ggsave(file.path(FIG_DIR,"Figure4_Youden_Pooled.png"), width=8, height=5)
message("✓ Figure 4 generated.")

# ============================================================
# FIGURE 5 — Early Warning Classification
# ============================================================
df_eval <- df %>% filter(YR %in% eval_years) %>% mutate(Alert=RF>=opt_pooled$RF_threshold)
p_fig5 <- ggplot(df_eval, aes(Date, RF)) +
  geom_col(aes(fill=Alert), alpha=0.7) +
  geom_hline(yintercept=opt_pooled$RF_threshold, linetype="dashed", linewidth=1) +
  scale_fill_manual(values=c("FALSE"="grey70","TRUE"="red"), labels=c("No Alert","Alert")) +
  labs(title="Figure 5: Rainfall-Based Early Warning Classification",
       subtitle=paste0("Evaluation period using pooled threshold = ", round(opt_pooled$RF_threshold,1)," mm"),
       x="Date", y="Weekly Rainfall (mm)", fill="Warning Status") +
  theme_minimal()
print(p_fig5)
ggsave(file.path(FIG_DIR,"Figure5_EarlyWarning.png"), width=12, height=6)
message("✓ Figure 5 generated.")

# ============================================================
# FIGURE 6 — DLNM: Lagged Effect of Rainfall
# ============================================================
lag_max <- 8
cb_rf <- crossbasis(df$RF, lag=lag_max, argvar=list(fun="lin"), arglag=list(fun="ns", df=3))
mdl_dlnm <- glm(DC ~ cb_rf + factor(YR), family=quasipoisson(), data=df)
pred <- crosspred(cb_rf, mdl_dlnm, at=quantile(df$RF, probs=seq(0.1,0.9,0.1), na.rm=TRUE), bylag=1)
plot(pred,"contour", xlab="Rainfall (mm)", ylab="Lag (weeks)", main="Figure 6: Lagged Effect of Rainfall on Dengue Incidence")
ggsave(file.path(FIG_DIR,"Figure6_DLNMLagEffect.png"), width=8, height=6)
message("✓ Figure 6 generated.")

# ============================================================
# FIGURE 7 — Rolling Dengue Cases & Rainfall
# ============================================================
plot_years <- target_years[target_years %in% df$YR]
roll_thr <- rolling_threshold_from_map(df, rolling_map_no_pandemic)
df_plot <- df %>% filter(YR %in% plot_years, WN %in% wet_weeks) %>%
  mutate(WN_cont = as.numeric(factor(paste0(YR,"-",WN),
                                     levels=paste0(rep(plot_years, each=length(wet_weeks)),"-",wet_weeks)))) %>%
  left_join(roll_thr %>% select(YR, Threshold), by="YR")

rf_min <- 175; rf_max <- 200
heavy_rows <- which(df_plot$RF >= rf_min)
lags <- 1:4
rain_detect <- bind_rows(lapply(lags, function(l){
  idx <- heavy_rows + l; idx <- idx[idx<=nrow(df_plot)]
  if(length(idx)==0) return(NULL)
  data.frame(WN_cont=df_plot$WN_cont[idx], DC=df_plot$DC[idx], Lag_type=ifelse(l==4,"4 weeks","1-3 weeks"))
}))
case_detect <- df_plot %>% filter(DC>=Threshold) %>% select(WN_cont,DC)
thr_annot <- df_plot %>% group_by(YR) %>% summarise(x=median(WN_cont), y=max(Threshold,na.rm=TRUE),
                                                    label=paste0(round(unique(Threshold),0)," cases")) %>% ungroup()

pC_base <- ggplot(df_plot, aes(x=WN_cont)) +
  geom_col(aes(y=RF, fill=ifelse(RF>=rf_min,paste0("≥",rf_min," mm"),paste0("<",rf_min," mm"))), width=0.8) +
  geom_line(aes(y=DC, color="Dengue cases"), linewidth=1.2) +
  geom_line(aes(y=Threshold, color="Rolling case threshold"), linetype="dashed", linewidth=1) +
  geom_point(data=case_detect, aes(WN_cont,DC,color="Crossing case threshold"), shape=4, size=3, stroke=1.2) +
  geom_point(data=rain_detect, aes(WN_cont,DC,color=Lag_type,shape=Lag_type), size=3) +
  geom_hline(yintercept=c(rf_min,rf_max), linetype="dashed", color="blue") +
  geom_text(data=thr_annot, aes(x=x,y=y,label=label), color="red", fontface="bold", size=3.6, vjust=-0.8) +
  labs(title="Figure 7: Dengue Cases and Weekly Rainfall with Rolling Thresholds",
       x="Year (Wet Season Weeks)", y="Cases / Weekly Rainfall (mm)", fill="Weekly Rainfall") +
  scale_fill_manual(values=c("lightblue","darkblue")) +
  scale_color_manual(values=c("Dengue cases"="red","Rolling case threshold"="darkgreen",
                              "Crossing case threshold"="red","1-3 weeks"="orange","4 weeks"="#004d00")) +
  scale_shape_manual(values=c("1-3 weeks"=16,"4 weeks"=24)) +
  scale_x_continuous(breaks=df_plot %>% group_by(YR) %>% summarise(mid=median(WN_cont)) %>% pull(mid),
                     labels=df_plot %>% group_by(YR) %>% summarise(mid=median(WN_cont)) %>% pull(YR)) +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold", hjust=0.5))
print(pC_base)
ggsave(file.path(FIG_DIR,"Figure7_RollingDengue.png"), width=12, height=6)
message("✓ Figure 7 generated.")

# ============================================================
# FIGURE 8 — Weekly Pediatric Admissions
# ============================================================
years_to_plot <- c(2019,2022:2025)
weekly_summary <- df_raw %>%
  filter(YR %in% years_to_plot, WN %in% wet_weeks) %>%
  group_by(YR, WN) %>%
  summarise(Cases=sum(PedAdmit,na.rm=TRUE), .groups="drop") %>%
  mutate(Week_cont=rep(1:17,times=length(years_to_plot)))
week_labels <- rep(23:39, times=length(years_to_plot))
weekly_plot <- ggplot(weekly_summary, aes(x=Week_cont,y=Cases,color=as.factor(YR),group=YR)) +
  geom_line(linewidth=1.2)+geom_point(size=2) +
  scale_x_continuous(breaks=1:(17*length(years_to_plot)), labels=week_labels) +
  labs(title="Figure 8: Weekly Pediatric Admissions (Weeks 23–39, 2019, 2022–2025)",
       x="Week Number", y="Number of Pediatric Admissions", color="Year") +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
print(weekly_plot)
ggsave(file.path(FIG_DIR,"Figure8_PedAdmissions.png"), width=12, height=6)
message("✓ Figure 8 generated.")

# ============================================================
# FIGURE 9 — Weekly Pediatric Case Fatality Rate (CFR)
# ============================================================
weekly_summary <- weekly_summary %>%
  left_join(df_raw %>% select(YR,WN,PedDead), by=c("YR","WN")) %>%
  mutate(CFR=ifelse(Cases==0,0,(PedDead/Cases)*100))
cfr_threshold <- 0.1; annotation_y <- cfr_threshold*2
cfr_plot <- ggplot(weekly_summary, aes(x=Week_cont, y=CFR)) +
  geom_line(size=1.2,color="black") +
  geom_point(size=2,color="red") +
  geom_hline(yintercept=cfr_threshold, linetype="dashed", color="red", size=0.8) +
  annotate("text", x=min(weekly_summary$Week_cont)+0.5, y=annotation_y,
           label=paste0("CFR Target = ", cfr_threshold, "%"), color="red", fontface="bold", size=4, hjust=0) +
  scale_x_continuous(breaks=weekly_summary$Week_cont,
                     labels=sapply(strsplit(paste0(weekly_summary$WN,"-",weekly_summary$YR), "-"),
                                   function(z) paste0("W", z[1], "\n", substring(z[2],3,4)))) +
  labs(title="Figure 9: Weekly Pediatric Case Fatality Rate (CFR) with Threshold",
       x="Week Number (top) / Year (bottom)", y="CFR (%)") +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        axis.title=element_text(face="bold"),
        axis.text.x=element_text(size=9),
        plot.margin=margin(5,5,5,5))
print(cfr_plot)
ggsave(file.path(FIG_DIR,"Figure9_PedCFR.png"), width=12, height=6)
message("✓ Figure 9 generated.")

# ============================================================
# FIGURE 10 — Regional Dengue & Rainfall (Rolling Thresholds & Pooled Youden)
# ============================================================
# NOTE: To replicate for other regions, just replace "BARMM" in the sheet name below
#       with the corresponding region name in your Excel file.
PATH_rainfall <- here::here("Rainfall-Dengue Cases.xlsx")
region_name <- "BARMM"  # <--- Change this to other region for replication

df_region <- read_excel(PATH_rainfall, sheet=region_name) %>%
  mutate(YR = as.integer(YR),
         WN = as.integer(WN),
         DC = as.numeric(DC),
         RF = as.numeric(RF),
         ISOweek = sprintf("%d-W%02d", YR, WN),
         Date = ISOweek2date(paste0(ISOweek,"-1"))) %>%
  arrange(Date)

all_years <- unique(df_region$YR)
plot_years <- all_years[all_years >= min(all_years)+1]  # e.g., skip first year for evaluation
wet_weeks <- 23:39
df_region <- df_region %>% filter(WN %in% wet_weeks)

# Continuous week index
df_region <- df_region %>% mutate(WN_cont = as.numeric(factor(paste(YR,WN,sep="-"),
                                                              levels = expand.grid(YR = all_years, WN = wet_weeks) %>%
                                                                arrange(YR,WN) %>% transmute(lvl = paste(YR,WN,sep="-")) %>% pull(lvl))))

# Rolling case threshold
case_thresholds <- rolling_threshold_from_map(df_region, rolling_map_no_pandemic)
df_region <- df_region %>% left_join(case_thresholds %>% select(YR, Threshold), by="YR")

# Detection points & RF threshold
rf_min_fixed <- 175
rf_max_fixed <- 200
case_detect <- df_region %>% filter(DC >= Threshold, YR %in% plot_years) %>% select(WN_cont, DC)
if(nrow(case_detect) == 0) case_detect <- data.frame(WN_cont = numeric(), DC = numeric())

heavy_rows <- which(df_region$RF >= rf_min_fixed)
lags <- 1:4
rain_detect_fixed <- bind_rows(lapply(lags, function(l){
  idx <- heavy_rows + l
  idx <- idx[idx <= nrow(df_region)]
  if(length(idx) == 0) return(NULL)
  data.frame(WN_cont = df_region$WN_cont[idx], DC = df_region$DC[idx], Lag = l)
}))
if(nrow(rain_detect_fixed) > 0){
  rain_detect_fixed <- rain_detect_fixed %>%
    filter(WN_cont %in% df_region$WN_cont[df_region$YR %in% plot_years]) %>%
    mutate(Lag_type = ifelse(Lag == 4, "4 weeks after RF threshold", "1–3 weeks after RF threshold"))
} else {
  rain_detect_fixed <- data.frame(WN_cont = numeric(), DC = numeric(), Lag = numeric(), Lag_type = character())
}

# Prepare plot
df_plot <- df_region %>% filter(YR %in% plot_years)
subtitle_fixed <- paste0(
  "QC-Calibrated RF Threshold: ", sprintf("%.2f", rf_min_fixed), "-", sprintf("%.2f", rf_max_fixed),
  " mm; Pooled Rolling Case Threshold: ", round(min(case_thresholds$Threshold),1), "-", round(max(case_thresholds$Threshold),1), " cases"
)

rf_min_scaled_fixed <- rf_min_fixed * max(df_plot$DC, na.rm=TRUE) / max(df_plot$RF, na.rm=TRUE)
rf_max_scaled_fixed <- rf_max_fixed * max(df_plot$DC, na.rm=TRUE) / max(df_plot$RF, na.rm=TRUE)

year_breaks <- df_plot %>% group_by(YR) %>% summarise(x = min(WN_cont)) %>% pull(x)

pFixed <- ggplot(df_plot, aes(WN_cont)) +
  geom_col(aes(y = RF * max(DC, na.rm=TRUE)/max(RF, na.rm=TRUE), fill = RF >= rf_min_fixed), width = 0.8, alpha = 0.7) +
  geom_line(aes(y = DC, color = "Dengue cases"), linewidth = 1.2) +
  geom_line(aes(y = Threshold, color = "Rolling case threshold"), linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = c(rf_min_scaled_fixed, rf_max_scaled_fixed), linetype = "dashed", color = "blue") +
  geom_point(data = case_detect, aes(WN_cont, DC, color = "Crossing case threshold"), shape = 4, size = 3, stroke = 1.2) +
  geom_point(data = rain_detect_fixed, aes(WN_cont, DC, color = Lag_type, shape = Lag_type), size = 3) +
  scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "darkblue"), name = "Weekly Rainfall") +
  scale_color_manual(values = c(
    "Dengue cases" = "red",
    "Rolling case threshold" = "darkgreen",
    "Crossing case threshold" = "red",
    "1–3 weeks after RF threshold" = "orange",
    "4 weeks after RF threshold" = "#004d00"
  )) +
  scale_shape_manual(values = c("1–3 weeks after RF threshold" = 16, "4 weeks after RF threshold" = 24)) +
  scale_x_continuous(breaks = year_breaks, labels = plot_years, expand = c(0,0)) +
  scale_y_continuous(
    name = "Weekly dengue cases",
    sec.axis = sec_axis(~ . * max(df_plot$RF, na.rm=TRUE)/max(df_plot$DC, na.rm=TRUE),
                        name = paste0("Weekly rainfall (mm) — Fixed RF: ", sprintf("%.2f", rf_min_fixed), "-", sprintf("%.2f", rf_max_fixed), " mm"))
  ) +
  labs(
    title = paste0(region_name, " (", min(plot_years), "-", max(plot_years), ") — Dengue Increase Detection (QC-Calibrated RF Threshold)"),
    subtitle = subtitle_fixed,
    x = "Year"
  ) +
  theme_minimal(base_size=12)

print(pFixed)
message(paste0("✓ Figure 10 generated for ", region_name))