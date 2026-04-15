# app.R
# 运行：shiny::runApp("app.R")  或在 RStudio 点 Run App

library(shiny)
library(survival)

# =========================
# 1) 数据与模型（后端安全名）
# =========================
data_file <- "imputed_data.CSV"
if (!file.exists(data_file)) {
  stop("未找到数据文件：", data_file)
}

train.data <- read.csv(data_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

vars_needed <- c("V69", "V66", "V43", "V44", "V48", "V57", "V67", "V83", "V92", "V95")
missing_vars <- setdiff(vars_needed, names(train.data))
if (length(missing_vars) > 0) {
  stop("数据缺少必要变量：", paste(missing_vars, collapse = ", "))
}

clean_data <- na.omit(train.data[, vars_needed, drop = FALSE])
clean_data <- as.data.frame(clean_data)

if (nrow(clean_data) == 0) {
  stop("清洗后数据为空（全部缺失或筛除），请检查输入数据。")
}

names(clean_data) <- c(
  "time", "status",
  "Perineural_invasion_n",
  "Microvascular_invasion_n",
  "AJCC_stage_n",
  "Postoperative_hemorrhage_n",
  "Histological_subtype_n",
  "Intraoperative_blood_transfusion_ml",
  "Operative_time_min",
  "Age_y"
)

# status 尽量转成 0/1
to_binary01 <- function(x) {
  # 因子/字符 -> 两类映射为0/1
  if (is.factor(x) || is.character(x)) {
    x <- as.factor(x)
    if (nlevels(x) != 2) {
      stop("status 不是二分类（因子/字符水平不为2），无法自动转为0/1。")
    }
    return(as.numeric(x) - 1)
  }
  
  # 数值 -> 若本身是0/1直接用；若是两类值则最小映射0、最大映射1
  x_num <- as.numeric(x)
  u <- sort(unique(x_num[!is.na(x_num)]))
  if (length(u) != 2) {
    stop("status 不是二分类（唯一值个数不等于2），无法自动转为0/1。")
  }
  if (all(u %in% c(0, 1))) {
    return(x_num)
  } else {
    return(ifelse(x_num == max(u), 1, 0))
  }
}

clean_data$status <- to_binary01(clean_data$status)
clean_data$time <- as.numeric(clean_data$time)

# 分类/连续变量类型
cat_vars <- c(
  "Perineural_invasion_n",
  "Microvascular_invasion_n",
  "AJCC_stage_n",
  "Postoperative_hemorrhage_n",
  "Histological_subtype_n"
)
num_vars <- c(
  "Intraoperative_blood_transfusion_ml",
  "Operative_time_min",
  "Age_y"
)

for (v in cat_vars) clean_data[[v]] <- as.factor(clean_data[[v]])
for (v in num_vars) clean_data[[v]] <- as.numeric(clean_data[[v]])

# Cox模型
cox_model <- coxph(
  Surv(time, status) ~
    Perineural_invasion_n +
    Microvascular_invasion_n +
    AJCC_stage_n +
    Postoperative_hemorrhage_n +
    Histological_subtype_n +
    Intraoperative_blood_transfusion_ml +
    Operative_time_min +
    Age_y,
  data = clean_data,
  x = TRUE, y = TRUE
)

# 基线累计风险
bh <- basehaz(cox_model, centered = FALSE)

get_H0 <- function(t) {
  approx(
    x = bh$time, y = bh$hazard, xout = t,
    method = "linear", rule = 2, ties = "ordered"
  )$y
}

# 安全名 -> 显示名
pretty_name <- c(
  Perineural_invasion_n = "Perineural invasion, n",
  Microvascular_invasion_n = "Microvascular invasion, n",
  AJCC_stage_n = "AJCC stage, n",
  Postoperative_hemorrhage_n = "Postoperative hemorrhage, n",
  Histological_subtype_n = "Histological subtype, n",
  Intraoperative_blood_transfusion_ml = "Intraoperative blood transfusion, ml",
  Operative_time_min = "Operative time, min",
  Age_y = "Age, y"
)

# 近似 points 缩放（用于可视化，不是 regplot/rms 的原生 points）
coef_vec <- coef(cox_model)

range_effect <- function(v) {
  if (v %in% num_vars) {
    b <- coef_vec[v]
    if (is.na(b)) return(0)
    rr <- range(clean_data[[v]], na.rm = TRUE)
    abs(b) * (rr[2] - rr[1])
  } else {
    idx <- grep(paste0("^", v), names(coef_vec))
    if (length(idx) == 0) return(0)
    max(c(0, abs(coef_vec[idx])), na.rm = TRUE)
  }
}
max_eff <- max(sapply(c(cat_vars, num_vars), range_effect), na.rm = TRUE)
point_scale <- ifelse(is.finite(max_eff) && max_eff > 0, 100 / max_eff, 1)

# =========================
# 2) UI（前端显示真实名）
# =========================
ui <- fluidPage(
  titlePanel("Dynamic Nomogram for DA-OS"),
  tags$hr(),
  sidebarLayout(
    sidebarPanel(
      h4("Patient Inputs"),
      selectInput(
        "Perineural_invasion_n", "Perineural invasion, n",
        choices = levels(clean_data$Perineural_invasion_n),
        selected = levels(clean_data$Perineural_invasion_n)[1]
      ),
      selectInput(
        "Microvascular_invasion_n", "Microvascular invasion, n",
        choices = levels(clean_data$Microvascular_invasion_n),
        selected = levels(clean_data$Microvascular_invasion_n)[1]
      ),
      selectInput(
        "AJCC_stage_n", "AJCC stage, n",
        choices = levels(clean_data$AJCC_stage_n),
        selected = levels(clean_data$AJCC_stage_n)[1]
      ),
      selectInput(
        "Postoperative_hemorrhage_n", "Postoperative hemorrhage, n",
        choices = levels(clean_data$Postoperative_hemorrhage_n),
        selected = levels(clean_data$Postoperative_hemorrhage_n)[1]
      ),
      selectInput(
        "Histological_subtype_n", "Histological subtype, n",
        choices = levels(clean_data$Histological_subtype_n),
        selected = levels(clean_data$Histological_subtype_n)[1]
      ),
      sliderInput(
        "Intraoperative_blood_transfusion_ml", "Intraoperative blood transfusion, ml",
        min = floor(min(clean_data$Intraoperative_blood_transfusion_ml, na.rm = TRUE)),
        max = ceiling(max(clean_data$Intraoperative_blood_transfusion_ml, na.rm = TRUE)),
        value = round(median(clean_data$Intraoperative_blood_transfusion_ml, na.rm = TRUE)),
        step = 1
      ),
      sliderInput(
        "Operative_time_min", "Operative time, min",
        min = floor(min(clean_data$Operative_time_min, na.rm = TRUE)),
        max = ceiling(max(clean_data$Operative_time_min, na.rm = TRUE)),
        value = round(median(clean_data$Operative_time_min, na.rm = TRUE)),
        step = 1
      ),
      sliderInput(
        "Age_y", "Age, y",
        min = floor(min(clean_data$Age_y, na.rm = TRUE)),
        max = ceiling(max(clean_data$Age_y, na.rm = TRUE)),
        value = round(median(clean_data$Age_y, na.rm = TRUE)),
        step = 1
      ),
      actionButton("go", "Predict", class = "btn-primary")
    ),
    mainPanel(
      h4("Predicted Survival"),
      tableOutput("pred_table"),
      tags$br(),
      h4("Risk Stratification"),
      verbatimTextOutput("risk_text"),
      tags$br(),
      h4("Nomogram-like Variable Contribution (Approx. Points)"),
      plotOutput("contrib_plot", height = "420px"),
      tags$small("Note: Points here are custom scaled for visualization and may differ from regplot/rms native points.")
    )
  )
)

# =========================
# 3) Server
# =========================
server <- function(input, output, session) {
  
  new_patient <- eventReactive(input$go, {
    data.frame(
      Perineural_invasion_n = factor(input$Perineural_invasion_n, levels = levels(clean_data$Perineural_invasion_n)),
      Microvascular_invasion_n = factor(input$Microvascular_invasion_n, levels = levels(clean_data$Microvascular_invasion_n)),
      AJCC_stage_n = factor(input$AJCC_stage_n, levels = levels(clean_data$AJCC_stage_n)),
      Postoperative_hemorrhage_n = factor(input$Postoperative_hemorrhage_n, levels = levels(clean_data$Postoperative_hemorrhage_n)),
      Histological_subtype_n = factor(input$Histological_subtype_n, levels = levels(clean_data$Histological_subtype_n)),
      Intraoperative_blood_transfusion_ml = as.numeric(input$Intraoperative_blood_transfusion_ml),
      Operative_time_min = as.numeric(input$Operative_time_min),
      Age_y = as.numeric(input$Age_y)
    )
  }, ignoreInit = FALSE)
  
  pred_obj <- reactive({
    nd <- new_patient()
    
    lp_fit <- predict(cox_model, newdata = nd, type = "lp", se.fit = TRUE)
    lp <- as.numeric(lp_fit$fit)
    se_lp <- as.numeric(lp_fit$se.fit)
    
    tvec <- c(12, 36, 60)
    H0 <- sapply(tvec, get_H0)
    
    # S(t|x) = exp(-H0(t) * exp(lp))
    S <- exp(-H0 * exp(lp))
    
    # 近似95%CI（由 lp 的 CI 映射）
    z <- 1.96
    S_low  <- exp(-H0 * exp(lp + z * se_lp))  # 更高风险 -> 更低生存
    S_high <- exp(-H0 * exp(lp - z * se_lp))
    
    list(
      nd = nd, lp = lp, se_lp = se_lp,
      times = tvec, S = S, S_low = S_low, S_high = S_high
    )
  })
  
  output$pred_table <- renderTable({
    p <- pred_obj()
    validate(need(all(is.finite(p$S)), "预测失败：请检查输入范围或模型数据。"))
    data.frame(
      Horizon = c("1-year OS (12m)", "3-year OS (36m)", "5-year OS (60m)"),
      `Predicted OS` = sprintf("%.3f", p$S),
      `95% CI` = paste0(sprintf("%.3f", p$S_low), " - ", sprintf("%.3f", p$S_high))
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")
  
  output$risk_text <- renderText({
    p5 <- pred_obj()$S[3]
    grp <- ifelse(p5 < 0.40, "High risk", "Low risk")
    paste0(
      "5-year OS = ", sprintf("%.3f", p5),
      " | Cutoff = 0.40 | Group: ", grp
    )
  })
  
  output$contrib_plot <- renderPlot({
    p <- pred_obj()
    nd <- p$nd
    
    mm <- model.matrix(delete.response(terms(cox_model)), data = nd)
    cf <- coef(cox_model)
    
    common <- intersect(colnames(mm), names(cf))
    validate(need(length(common) > 0, "无法计算变量贡献（模型矩阵与系数不匹配）。"))
    
    contrib <- mm[1, common] * cf[common]
    
    # 汇总到变量级别
    raw_names <- names(contrib)
    var_base <- raw_names
    for (v in c(cat_vars, num_vars)) {
      var_base <- sub(paste0("^", v, ".*$"), v, var_base)
    }
    
    df <- aggregate(as.numeric(contrib), by = list(var = var_base), FUN = sum)
    names(df)[2] <- "lp_contrib"
    
    # 转 points（近似）
    df$points <- df$lp_contrib * point_scale
    df$label <- pretty_name[df$var]
    df$label[is.na(df$label)] <- df$var
    
    # 排序便于阅读
    df <- df[order(df$points), , drop = FALSE]
    
    op <- par(mar = c(5, 12, 4, 2))
    on.exit(par(op), add = TRUE)
    
    cols <- ifelse(df$points >= 0, "#D62728", "#1F77B4")
    barplot(
      df$points,
      horiz = TRUE,
      col = cols,
      border = NA,
      names.arg = df$label,
      las = 1,
      xlab = "Approx. Points Contribution",
      main = paste0("Total Points (approx): ", sprintf("%.1f", sum(df$points)))
    )
    abline(v = 0, lty = 2)
  })
}

shinyApp(ui = ui, server = server)