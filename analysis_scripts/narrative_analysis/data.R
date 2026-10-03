# Figure 5: data module.
dat <- read_qualitative_data(INPUT_CSV)
required <- unique(c("Study", "Mortality", "Mortality improvement", process_spec$process_metrics, alert_spec$feature_fields, population_spec$feature_fields))
stopifnot(all(required %in% names(dat)), !anyDuplicated(dat$Study))
if (!identical(dat$Mortality, dat[["Mortality improvement"]])) {
  print(dat[which((is.na(dat$Mortality) != is.na(dat[["Mortality improvement"]])) |
    (!is.na(dat$Mortality) & !is.na(dat[["Mortality improvement"]]) & dat$Mortality != dat[["Mortality improvement"]])),
    c("Study", "Mortality", "Mortality improvement")]); stop("Mortality column mismatch")
}
raw_setting <- dat[["Deployment setting"]]
# Mortality-facing characteristics retain the four subgroup-analysis settings.
# The outer process layer continues to use its separate two-category definition.
stopifnot(all(is.na(raw_setting) | raw_setting %in% c("ICU","Wards","ED","Multiple settings","NR")))
dat[["Process deployment setting"]] <- collapse_setting(raw_setting)
save_csv(data.frame(Study=dat$Study,Source_setting=raw_setting,
  Mortality_category=raw_setting,Process_subgroup=dat[["Process deployment setting"]]),"deployment_setting_study_qa")
alert <- alert_spec$count_categories(alert_spec$build_long_data(dat))
population <- population_spec$count_categories(population_spec$build_long_data(dat))
process <- process_spec$count_metrics(dat, process_spec$process_metrics, "Process")
mortality <- process_spec$count_metrics(dat, "Mortality", "Outcome")
for (m in c(process_spec$process_metrics,"Mortality")) {
  ref <- if(m=="Mortality") mortality else process[process$Metric==m,]
  vals <- dat[[m]]
  stopifnot(ref$S==sum(vals=="S",na.rm=TRUE),ref$N==sum(vals=="N",na.rm=TRUE),
            ref$Reported==sum(vals %in% c("S","N")))
}
message("PASS: independent process and mortality totals")


assert_same <- function(actual, expected, keys, values, tag) {
  normalize <- function(d) {
    d <- as.data.frame(d[, c(keys, values), drop = FALSE])
    d[] <- lapply(d, function(x) if (is.factor(x)) as.character(x) else x)
    d <- d[do.call(order, d[keys]), , drop = FALSE]; rownames(d) <- NULL; d
  }
  a <- normalize(actual); b <- normalize(expected)
  if (!isTRUE(all.equal(a, b, tolerance = 1e-10, check.attributes = FALSE))) {
    print(a); print(b); stop("QA mismatch: ", tag)
  }
  message("PASS: ", tag)
}
# Independent direct CSV tabulation is compared with the shared aggregation functions.
for (side in c("alert", "population")) {
  counts <- get(side)
  direct <- bind_rows(lapply(seq_len(nrow(counts)), function(i) {
    r <- counts[i, ]; vals <- dat[["Mortality improvement"]][!is.na(dat[[as.character(r$Feature)]]) & dat[[as.character(r$Feature)]] == r$Category]
    data.frame(Feature = as.character(r$Feature), Category = r$Category,
      S = sum(vals == "S", na.rm = TRUE), N = sum(vals == "N", na.rm = TRUE))
  }))
  assert_same(direct, counts, c("Feature", "Category"), c("S", "N"), paste(side, "current script"))
}
# Study-level antibiotic administration: either evaluable endpoint S wins;
# N requires at least one evaluable endpoint and no S. Missing remains missing.
combine_delivery <- function(proportion, time) {
  has_s <- proportion %in% "S" | time %in% "S"
  has_n <- proportion %in% "N" | time %in% "N"
  ifelse(has_s, "S", ifelse(has_n, "N", NA_character_))
}
truth <- expand.grid(Proportion=c("S","N",NA_character_),Time=c("S","N",NA_character_),stringsAsFactors=FALSE)
stopifnot(identical(combine_delivery(truth$Proportion,truth$Time),
  c("S","S","S","S","N","N","S","N",NA_character_)))
individual_process_counts <- process
stopifnot("Antibiotic time" %in% names(dat))
dat[["Antibiotic administration"]] <- combine_delivery(dat[["Antibiotic process"]],dat[["Antibiotic time"]])
delivery_qa <- data.frame(Study=dat$Study,Proportion=dat[["Antibiotic process"]],
  Time=dat[["Antibiotic time"]],Delivery=dat[["Antibiotic administration"]],Mortality=dat$Mortality)
independent_delivery <- apply(delivery_qa[c("Proportion","Time")],1,function(v) {
  evaluable <- v[!is.na(v) & v %in% c("S","N")]
  if (!length(evaluable)) NA_character_ else if (any(evaluable=="S")) "S" else "N"
})
stopifnot(identical(unname(independent_delivery),dat[["Antibiotic administration"]]),!anyDuplicated(delivery_qa$Study))
save_csv(delivery_qa,"antibiotic_delivery_study_qa")
process_spec$process_metrics[process_spec$process_metrics=="Antibiotic process"] <- "Antibiotic administration"
process_spec$metric_labels["Antibiotic administration"] <- "Antibiotic administration"
process <- process_spec$count_metrics(dat,process_spec$process_metrics,"Process")
assert_same(process |> filter(Metric!="Antibiotic administration"),
  individual_process_counts |> filter(Metric!="Antibiotic process"),"Metric",c("S","N","Reported"),"four unchanged process domains")
message("PASS: delivery truth table, independent study-level union and unique studies")
links <- bind_rows(lapply(process$Metric, function(p) bind_rows(lapply(c("S", "N"), function(s) {
  vals <- dat$Mortality[!is.na(dat[[p]]) & dat[[p]] == s & dat$Mortality %in% c("S", "N")]
  data.frame(Process = p, Process_status = s, Count = length(vals), Positive = sum(vals == "S"),
    Negative = sum(vals == "N"), Positive_rate = if (length(vals)) mean(vals == "S") else NA_real_)
})))) |> filter(Count > 0)
# Compare each process status with the opposite status among joint-evaluable studies.
links <- links |> group_by(Process) |> mutate(
  Other_positive=sum(Positive)-Positive,Other_negative=sum(Negative)-Negative,
  Fisher_p=mapply(function(s,n,os,on) {
    mat <- matrix(c(s,n,os,on),nrow=2,byrow=TRUE)
    if(all(rowSums(mat)>0)&&all(colSums(mat)>0)) fisher.test(mat)$p.value else NA_real_
  },Positive,Negative,Other_positive,Other_negative),
  Rate_difference=abs(Positive_rate-Other_positive/(Other_positive+Other_negative)),
  Highlighted=(!is.na(Fisher_p)&Fisher_p<.05)|(!is.na(Rate_difference)&Rate_difference>=.25)) |>
  ungroup() |> select(-Other_positive,-Other_negative)
stopifnot(all(links$Count==links$Positive+links$Negative))

# Study subgroups are evaluated within each process.
# Each parent includes unique studies with evaluable parent-process S/N.
# NR is excluded from each subgroup; parent process S/N denominators stay unchanged.
# Process baseline comes from the exact source classifications used by S5H.
# Lower completion and longer antibiotic time indicate worse performance.
baseline_levels <- c("Low baseline","Moderate baseline","High baseline")
binary_levels <- c("Low baseline (<33.3%)","Moderate baseline (33.3% to <66.7%)","High baseline (>=66.7%)")
time_levels <- c("Long baseline time (>=6 h)","Moderate baseline time (2 to <6 h)","Short baseline time (<2 h)")
study_key <- function(x) tolower(gsub("[[:space:]]+","",trimws(x)))
stopifnot(!anyDuplicated(study_key(dat$Study)))
read_process_baseline <- function(name,time=FALSE) {
  metric <- if(time) "Antibiotic time" else switch(name,
    blood_culture="Blood culture completion",lactate_measurement="Lactate testing completion",
    antibiotic_proportion="Antibiotic process",fluid_resuscitation="Intravenous fluid completion",
    sepsis_bundle="Overall bundle completion")
  field <- paste0("Process baseline - ",metric)
  stopifnot(field %in% names(dat))
  values <- dat[[field]]
  stopifnot(all(is.na(values) | values %in% c(baseline_levels,"NR")))
  rank <- match(values,baseline_levels)
  number_field <- if(time) "Antibiotic time control baseline (h)" else
    paste0("Control baseline value (%) - ",metric)
  number <- suppressWarnings(as.numeric(dat[[number_field]]))
  numeric_rank <- if(time) ifelse(number>=6,1L,ifelse(number>=2,2L,3L)) else
    ifelse(number<33.3,1L,ifelse(number<66.7,2L,3L))
  stopifnot(all(is.na(number) | (!is.na(rank) & rank==numeric_rank)))
  rank
}
# The separate time CSV is an audit, never a restricted study universe.
time_audit <- read_qualitative_data(source_csv_path("narrative_analysis_antibiotic_time_audit.csv",root=PROJECT_ROOT))
stopifnot(!anyDuplicated(time_audit$Study),setequal(time_audit$Study,dat$Study))
time_audit <- time_audit[match(dat$Study,time_audit$Study),]
stopifnot(identical(time_audit[["Antibiotic time S/N"]],dat[["Antibiotic time"]]),
  identical(time_audit[["Baseline class"]],dat[["Process baseline - Antibiotic time"]]))
source_names <- c("blood_culture","lactate_measurement","antibiotic_proportion","fluid_resuscitation","sepsis_bundle")
baseline_ranks <- setNames(lapply(source_names,read_process_baseline),process_spec$process_metrics)
proportion_rank <- baseline_ranks[["Antibiotic administration"]]
time_rank <- read_process_baseline("antibiotic_time",TRUE)
# Rank 1 is poorest baseline / greatest headroom. One known dimension suffices.
composite_rank <- pmin(proportion_rank,time_rank,na.rm=TRUE)
composite_rank[is.infinite(composite_rank)] <- NA_integer_
reference_rank <- vapply(seq_len(nrow(dat)),function(i) {
  v <- c(proportion_rank[i],time_rank[i]); v <- v[!is.na(v)]
  if(length(v)) sort(v)[1] else NA_integer_
},integer(1))
stopifnot(identical(as.integer(composite_rank),reference_rank))
stored_composite <- match(dat[["Process baseline - Antibiotic administration"]],baseline_levels)
stopifnot(identical(as.integer(composite_rank),stored_composite))
baseline_ranks[["Antibiotic administration"]] <- composite_rank
baseline_fields <- setNames(paste0("Process baseline - ",process_spec$process_metrics),process_spec$process_metrics)
for(p in process_spec$process_metrics) dat[[baseline_fields[[p]]]] <- baseline_levels[baseline_ranks[[p]]]
save_csv(data.frame(Study=dat$Study,Proportion_baseline=baseline_levels[proportion_rank],
  Time_baseline=baseline_levels[time_rank],Composite_baseline=baseline_levels[composite_rank],
  Outcome=dat[["Antibiotic administration"]]),"antibiotic_composite_baseline_study_qa")
save_csv(data.frame(Process=process_spec$process_metrics,
  Source="source_csv/narrative_analysis_combined_conclusions.csv",
  Source_field=paste0("Process baseline - ",process_spec$process_metrics),
  Time_source=ifelse(process_spec$process_metrics=="Antibiotic administration",
    "source_csv/narrative_analysis_antibiotic_time_audit.csv; time audit, poorer available proportion/time baseline", "")),
  "process_baseline_source_mapping")
schemes <- c("Control baseline subgroup","Deployment setting","Intervention package type")
subgroup_fields <- setNames(c("Baseline mortality in control group","Process deployment setting","Intervention package type"),schemes)
subgroup_orders <- list(baseline_levels,
  c("ED","Wards/ICU"),alert_spec$category_orders[[subgroup_fields[3]]])
names(subgroup_orders) <- schemes
schemes <- c("Intervention package type","Control baseline subgroup","Deployment setting")
mapping <- expand.grid(Process=process_spec$process_metrics,Scheme=schemes,stringsAsFactors=FALSE)
mapping$Field <- unname(subgroup_fields[mapping$Scheme])
baseline_rows <- mapping$Scheme=="Control baseline subgroup"
mapping$Field[baseline_rows] <- unname(baseline_fields[mapping$Process[baseline_rows]])
mapping$Scope <- "Unique studies with evaluable parent process S/N; subgroup NR excluded"
outer <- bind_rows(lapply(seq_len(nrow(mapping)),function(i) {
  p <- mapping$Process[i]; scheme <- mapping$Scheme[i]; field <- mapping$Field[i]
  keep <- dat[[p]] %in% c("S","N")
  v <- trimws(as.character(dat[[field]][keep])); v[is.na(v) | v==""] <- "NR"
  missing_subgroup <- sum(v=="NR")
  status <- dat[[p]][keep][v!="NR"]
  v <- v[v!="NR"]
  codes <- subgroup_orders[[scheme]]
  if (any(!v %in% codes)) {print(unique(v[!v %in% codes]));stop("Unexpected subgroup category: ",field)}
  a <- data.frame(Process=p,Scheme=scheme,Field=field,Code=codes,
    Count=vapply(codes,function(code) sum(v==code),integer(1)))
  a$S <- vapply(codes,function(code) sum(v==code & status=="S"),integer(1))
  a$N <- vapply(codes,function(code) sum(v==code & status=="N"),integer(1))
  a$S_rate <- ifelse(a$Count>0,a$S/a$Count,NA_real_)
  stopifnot(all(a$Count==a$S+a$N))
  reference <- table(factor(v,levels=codes))
  assert_same(a,data.frame(Code=names(reference),Count=as.integer(reference)),"Code","Count",paste(p,scheme))
  a$Evaluable <- sum(a$Count); a$Excluded <- sum(!keep)+missing_subgroup
  a$Unclassified <- missing_subgroup
  stopifnot(a$Evaluable[1]+missing_subgroup==process$Reported[process$Metric==p],!anyDuplicated(dat$Study[keep]))
  a
}))
stopifnot(nrow(unique(outer[c("Process","Scheme")]))==nrow(process)*length(schemes))
save_csv(alert, "alert_counts"); save_csv(population, "population_counts")
save_csv(mortality, "mortality_counts"); save_csv(process, "process_counts")
save_csv(links, "process_mortality_links"); save_csv(outer, "outer_classification_counts")
save_csv(mapping, "field_mapping")
