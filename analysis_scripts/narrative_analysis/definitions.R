# Stable data definitions and shared statistical/geometry primitives.
# Process and characteristic classification definitions.
clean_names <- function(x) trimws(sub("^\ufeff", "", x))
read_qualitative_data <- function(path) {
    if (!file.exists(path)) 
        stop("Missing input: ", path)
    dat <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
    names(dat) <- clean_names(names(dat))
    dat[] <- lapply(dat, function(x) {
        x <- trimws(as.character(x))
        x[x %in% c("", "NA", "Na", "na", "N/A", "/")] <- NA_character_
        x
    })
    dat
}
process_spec <- local({
process_metrics <- c("Blood culture completion", "Lactate testing completion", "Antibiotic process", 
    "Intravenous fluid completion", "Overall bundle completion")
metric_labels <- c(`Blood culture completion` = "Blood culture", `Lactate testing completion` = "Lactate measurement",
    `Antibiotic process` = "Antibiotic process", `Intravenous fluid completion` = "Fluid resuscitation", 
    `Overall bundle completion` = "Sepsis bundle", Mortality = "Mortality")
process_colors_deep <- NULL
process_colors_light <- NULL
count_metrics <- function(dat, metrics, side) {
    bind_rows(lapply(metrics, function(metric) {
        vals <- dat[[metric]]
        data.frame(Side = side, Metric = metric, Label = unname(metric_labels[metric]), S = sum(vals == 
            "S", na.rm = TRUE), N = sum(vals == "N", na.rm = TRUE), Reported = sum(vals %in% c("S", "N"), 
            na.rm = TRUE), stringsAsFactors = FALSE)
    }))
}
environment()
})
make_characteristic_spec <- function(feature_fields,feature_labels,feature_colors_deep,feature_colors_light,category_orders,exclude_all) {
build_long_data <- function(dat, mortality_evaluable = TRUE) {
    if (mortality_evaluable) dat <- filter(dat, .data[["Mortality improvement"]] %in% c("S", "N"))
    filter(bind_rows(lapply(feature_fields, function(feature) {
        data.frame(Study = dat$Study, Feature = feature, Feature_label = unname(feature_labels[feature]), 
            Category = dat[[feature]], Mortality = dat[["Mortality improvement"]], stringsAsFactors = FALSE)
    })), !is.na(Category), Category != "", Category != "NR", !exclude_all | tolower(Category) != "all")
}
count_categories <- function(long_dat) {
    feature_tables <- lapply(feature_fields, function(feature) {
        d <- filter(long_dat, Feature == feature)
        order_vec <- category_orders[[feature]]
        feature_counts <- d %>%
            count(Category, Mortality, name = "n") %>%
            tidyr::pivot_wider(names_from = Mortality, values_from = n, values_fill = 0) %>%
            mutate(
                Feature = feature,
                Feature_label = unname(feature_labels[feature]),
                S = ifelse(is.na(S), 0, S),
                N = ifelse(is.na(N), 0, N),
                Reported = S + N,
                Positive_rate = ifelse(Reported > 0, S / Reported, NA_real_),
                Fisher_p = mapply(function(s, n) {
                    rest_s <- sum(S, na.rm = TRUE) - s
                    rest_n <- sum(N, na.rm = TRUE) - n
                    mat <- matrix(c(s, n, rest_s, rest_n), nrow = 2, byrow = TRUE)
                    if (sum(mat) > 0 && all(rowSums(mat) > 0) && all(colSums(mat) > 0)) {
                        fisher.test(mat)$p.value
                    } else {
                        NA_real_
                    }
                }, S, N),
                Category_order = match(Category, order_vec),
                Category_order = ifelse(is.na(Category_order), 999 + row_number(), Category_order)
            )
        min_rate <- min(feature_counts$Positive_rate, na.rm = TRUE)
        max_rate <- max(feature_counts$Positive_rate, na.rm = TRUE)
        rate_range <- max_rate - min_rate
        feature_counts %>%
            mutate(
                Rate_difference = rate_range,
                Highlighted = (!is.na(Fisher_p) & Fisher_p < 0.05) |
                    (!is.na(rate_range) & rate_range >= 0.25 &
                     (abs(Positive_rate - min_rate) < 1e-12 |
                      abs(Positive_rate - max_rate) < 1e-12))
            ) %>%
            arrange(Category_order)
    })
    bind_rows(feature_tables) %>%
        mutate(Feature = factor(Feature, levels = feature_fields),
               Arc_id = paste(Feature, Category, sep = "___")) %>%
        arrange(Feature, Category_order)
}

environment()
}
alert_spec <- do.call(make_characteristic_spec, list(feature_fields = c("Deployment setting", "Alert targeting approach", 
"Intervention package type", "Algorithm class", "Delivery mode", 
"Recipient"), feature_labels = c(`Deployment setting` = "Deployment\nsetting", 
`Alert targeting approach` = "Alert targeting\napproach", `Intervention package type` = "Intervention\npackage", 
`Algorithm class` = "Algorithm\nclass", `Delivery mode` = "Delivery mode", 
Recipient = "Recipient"), feature_colors_deep = NULL, feature_colors_light = NULL, category_orders = list(`Deployment setting` = c("ICU", "Wards", "ED", "Multiple settings"
), `Alert targeting approach` = c("Suspected infection", "Sepsis with organ dysfunction", 
"Tiered sepsis alerts"), `Intervention package type` = c("Alert only", 
"Actionable alert", "Co-intervention"), `Algorithm class` = c("Rule", 
"ML"), `Delivery mode` = c("EHR display", "Page/message"), Recipient = c("Nurse", 
"Physician", "Physician and nurse")), exclude_all = FALSE))
population_spec <- do.call(make_characteristic_spec, list(feature_fields = c("Identification method", "Analytic population", 
"Baseline mortality in control group", "Total analytic sample size", 
"Center setting", "Design"), feature_labels = c(`Identification method` = "Population\nascertainment",
`Analytic population` = "Clinical population\nscope", `Baseline mortality in control group` = "Baseline mortality\nin control group",
`Total analytic sample size` = "Total analytic\nsample size", 
`Center setting` = "Centres/sites", Design = "Design"), feature_colors_deep = NULL, feature_colors_light = NULL, category_orders = list(`Identification method` = c("alert-triggered",
"alert + diagnosis", "diagnosis-based"), `Analytic population` = c("Infection", 
"Severe sepsis", "Sepsis-3"), `Baseline mortality in control group` = c("<10%", 
"10% to <20%", ">=20%"), `Total analytic sample size` = c("Small sample (<1000)", 
"Large sample (>=1000)"), `Center setting` = c("Single-center", 
"Multicenter"), Design = c("RCT", "NRS")), exclude_all = TRUE))
collapse_setting <- function(x) {
    dplyr::case_when(x %in% c("Wards", "ICU") ~ "Wards/ICU", x == "ED" ~ "ED", TRUE ~ NA_character_)
}
