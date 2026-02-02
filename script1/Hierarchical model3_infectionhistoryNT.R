library(readr)
library(dplyr)
library(rstan)
library(ggplot2)
library(cowplot)

supdata<-read_csv('supplemental_data_NT.csv')
Hmodel3 <- stan_model('hierarchical M3.stan')

d3_stan<- supdata

#data curation
d3_stan$infection_status0=case_when(d3_stan$infectionstatus3=="uninfected"~1,
                                    d3_stan$infectionstatus3=="Omicron"~2,
                                    d3_stan$infectionstatus3=="Pre-Omicron"~NA) #exclude 9 Pre-Omicron infected

d3_stan<-d3_stan %>% mutate(
  vax_infection=as.numeric(case_when(
    vax==2&infection_status0==1~1, 
    vax==3&infection_status0==1~2,
    vax==4&infection_status0==1~2,
    vax==5&infection_status0==1~2,
    vax==2&infection_status0==2~3,
    vax==3&infection_status0==2~4,
    vax==4&infection_status0==2~4,
    vax==5&infection_status0==2~4 
  )),
  
  vax_mixed=factor(
    case_when(
      vax == 2     ~ "Twice",
      vax >= 3     ~ "3-5 times"),
    levels = c("Twice", "3-5 times")
  ))
d3_stan<-d3_stan %>% filter(is.na(infection_status0)==FALSE,is.na(vax_infection)==FALSE)

infection_status<-unique(d3_stan[ , c('vax_infection','infection_status0')])$infection_status0

data_infectionstNT <-list(N=nrow(d3_stan),Y=d3_stan$NT,X=d3_stan$general_intervals,K=4,S=2,
                          VAX=as.numeric(d3_stan$vax_infection),
                          AGE=infection_status,
                          AGE0=d3_stan$infection_status0)

#fitting
fit_infectionstexpNT2 <- sampling(Hmodel3,
                                  data=data_infectionstNT,
                                  iter = 9000,
                                  warmup = 1000,
                                  thin = 5,
                                  chains = 4, core=4,
                                  seed=1234,
                                  control = list(adapt_delta = 0.99, max_treedepth=15)
                                  )

#visualise
ms <- rstan::extract(fit_infectionstexpNT2)

x_days  <- 0:730
q_probs <- c(0.025, 0.25, 0.50, 0.75, 0.975)

# 点プロット用
d3_stan <- d3_stan %>% mutate(Ylin = NT)

# ---- quantile summary (列名を扱いやすく) ----
ci_df <- function(mat, x, probs = q_probs) {
  mat <- as.matrix(mat)
  qs  <- apply(mat, 2, quantile, probs = probs)
  
  out <- data.frame(X = x, t(qs), check.names = FALSE)
  names(out) <- c("X", "q025", "q25", "q50", "q75", "q975")
  out
}

# ---- exp-decay curve draws: base = c + a*exp(-b*x) ----
draw_exp_curve <- function(c, a, b, s_y, x, seed = 1234, probs = q_probs) {
  c   <- drop(c)
  a   <- drop(a)
  b   <- drop(b)
  s_y <- drop(s_y)
  
  n  <- length(a)
  nx <- length(x)
  
  if (length(c) == 1) c <- rep(c, n)
  stopifnot(length(c) == n, length(b) == n)
  
  # sd は (scalar or length n) を許容
  sd_row <- if (length(s_y) == 1) rep(s_y, n) else rep(s_y, length.out = n)
  
  eta  <- exp(-outer(b, x, `*`))      # n x nx
  base <- sweep(eta, 1, a, `*`)       # a*exp(-b*x)
  base <- sweep(base, 1, c, `+`)      # + c
  
  set.seed(seed)
  eps <- rnorm(n * nx, mean = 0, sd = rep(sd_row, times = nx))
  y   <- base + matrix(eps, nrow = n, ncol = nx)
  
  list(
    base_draws = as.data.frame(base),
    pred_draws = as.data.frame(y),
    base_ci    = ci_df(base, x, probs),
    pred_ci    = ci_df(y, x, probs)
  )
}

# ---- 4 groups ----
vax_levels <- c("Twice", "3–5 times")
inf_levels <- c("Uninfected", "Omicron")

group_map <- expand.grid(
  infection_ix = 1:2,
  vax_ix       = 1:2
) %>%
  as_tibble() %>%
  mutate(
    group       = (infection_ix - 1) * 2 + vax_ix,  # 1..4
    vaccination = vax_levels[vax_ix],
    infection   = inf_levels[infection_ix],
    
    sd_col      = c(1, 1, 1, 2)[group]
    
  ) %>%
  arrange(group)

grp_res <- setNames(
  lapply(seq_len(nrow(group_map)), function(i) {
    g <- group_map$group[i]
    w <- group_map$sd_col[i]
    draw_exp_curve(
      c   = ms$c[, g],
      a   = ms$a[, g],
      b   = ms$b[, g],
      s_y = ms$s_y[, w],
      x   = x_days,
      seed = 1234
    )
  }),
  paste0("g", group_map$group)
)

base_ci <- lapply(grp_res, `[[`, "base_ci")   
pred_ci <- lapply(grp_res, `[[`, "pred_ci")   

# ---- base = a1*exp(-b1*x)（c=0） ----
age_res <- setNames(
  lapply(1:2, function(z) {
    draw_exp_curve(
      c   = 0,
      a   = ms$a1[, z],
      b   = ms$b1[, z],
      s_y = ms$s_y[, z],
      x   = x_days,
      seed = 1234
    )
  }),
  paste0("age", 1:2)
)
d_infection1.2 <- age_res$age1$pred_ci
d_infection2.2 <- age_res$age2$pred_ci

# ---- total: base = c0 + a0*exp(-b0*x) ----
total_res <- draw_exp_curve(
  c   = ms$c0,
  a   = ms$a0,
  b   = ms$b0,
  s_y = ms$s_y,
  x   = x_days,
  seed = 1234
)
d_estTotal <- total_res$base_ci
d_est      <- total_res$pred_ci


colp_inf <- c("#cccccc50", "#7744a4")

plot_nt_infection_pair <- function(ci_uninf, ci_omi, g_ids,
                                   ribbon_cols = c("#99999950", "#7744a4"),
                                   title = "",
                                   points_df = d3_stan,
                                   ylim = c(0, 4),
                                   hline = log10(5)) {
  pts <- points_df %>% filter(vax_infection %in% g_ids)
  
  ggplot() +
    theme_bw(base_size = 12) +
    geom_point(
      data = pts,
      aes(x = general_intervals, y = Ylin,
          color = factor(infection_status0),
          fill  = factor(infection_status0)),
      shape = 1, alpha = .6, size = 1
    ) +
    geom_ribbon(data = ci_uninf, aes(x = X, ymin = q025, ymax = q975),
                fill = ribbon_cols[1], alpha = 3/6) +
    geom_line(  data = ci_uninf, aes(x = X, y = q50),
                colour = ribbon_cols[1], size = 1, alpha = .8) +
    geom_ribbon(data = ci_omi, aes(x = X, ymin = q025, ymax = q975),
                fill = ribbon_cols[2], alpha = 3/6) +
    geom_line(  data = ci_omi, aes(x = X, y = q50),
                colour = ribbon_cols[2], size = 1, alpha = .8) +
    coord_cartesian(ylim = ylim) +
    scale_color_manual(values = colp_inf) +
    scale_fill_manual(values  = colp_inf) +
    labs(
      y = expression(paste("Log"["10"], " BA.5 NT")),
      x = "Days post exposure",
      color = "Infection Status",
      title = title
    ) +
    theme(
      axis.ticks = element_line(colour = "black"),
      axis.ticks.length = grid::unit(-3, "mm"),
      axis.text = element_text(colour = "black", size = 12),
      axis.title.x = element_text(colour = "black"),
      axis.title.y = element_text(colour = "black"),
      plot.title = element_text(colour = "black", size = 12),
      legend.position = "none",
      plot.tag = element_text(face = "bold")
    ) +
    geom_hline(yintercept = hline, colour = "black", linetype = 3, alpha = 0.8, size = 0.8)
}

Vp1 <- plot_nt_infection_pair(base_ci$g1, base_ci$g3, g_ids = c(1, 3), title = "Twice vaccinated")
Vp2 <- plot_nt_infection_pair(base_ci$g2, base_ci$g4, g_ids = c(2, 4), title = "3 - 5 times vaccinated")

VpiO <- plot_grid(Vp1, Vp2, ncol = 2)
VpiO

#make csv for infection risk estimation
est_all_nt <- bind_rows(
  lapply(seq_len(nrow(group_map)), function(i) {
    key <- paste0("g", group_map$group[i])
    base_ci[[key]] %>%
      transmute(
        vaccination = group_map$vaccination[i],
        infection   = group_map$infection[i],
        days = X,
        q025, q25, q50, q75, q975
      )
  })
) %>%
  mutate(
    vaccination = factor(vaccination, levels = vax_levels),
    infection   = factor(infection, levels = inf_levels)
  ) %>%
  relocate(vaccination, infection, days)

dplyr::glimpse(est_all_nt)
table(est_all_nt$vaccination, est_all_nt$infection)

write_excel_csv(est_all_nt, "BA5NT_est_all_by_group.csv")
