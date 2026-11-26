library(tidyverse)

# ==========================
# 1. DATA INLADEN
# ==========================

metadata <- read.csv("sample_metadata.csv")
proteins <- read.csv("proteins_npx_raw.csv")
annotations <- read.csv("protein_annotation.csv")

# ==========================
# 2. AANTAL SAMPLES EN FEATURES
# ==========================

cat("METADATA\n")
cat("Aantal samples:", nrow(metadata), "\n")
cat("Aantal metadata variabelen:", ncol(metadata), "\n\n")

cat("PROTEIN DATA\n")
cat("Aantal samples:", nrow(proteins), "\n")
cat("Aantal eiwitten/features:", ncol(proteins) - 1, "\n\n")
# -1 omdat meestal de eerste kolom sample ID bevat

cat("ANNOTATIONS\n")
cat("Aantal geannoteerde eiwitten:", nrow(annotations), "\n")

# ==========================
# 3. VARIABLE TYPES
# ==========================

variable_types <- data.frame(
  Variable = names(metadata),
  Class = sapply(metadata, class)
)

print(variable_types)

# Nominale variabelen vinden

nominal_vars <- names(metadata)[sapply(metadata, function(x)
  is.character(x) | is.factor(x))]

cat("\nNominale variabelen:\n")
print(nominal_vars)

# Numerieke variabelen

numeric_vars <- names(metadata)[sapply(metadata, is.numeric)]

cat("\nNumerieke variabelen:\n")
print(numeric_vars)

# ==========================
# 4. MISSING VALUES IN METADATA
# ==========================

missing_metadata <- data.frame(
  Variable = names(metadata),
  N_missing = colSums(is.na(metadata)),
  Pct_missing = round(
    colSums(is.na(metadata)) / nrow(metadata) * 100,
    2
  )
)

missing_metadata <-
  missing_metadata %>%
  arrange(desc(Pct_missing))

print(head(missing_metadata, 20))

# ==========================
# 5. MISSING VALUES IN PROTEIN DATA
# ==========================

protein_missing <- data.frame(
  Protein = names(proteins),
  N_missing = colSums(is.na(proteins)),
  Pct_missing = round(
    colSums(is.na(proteins)) /
      nrow(proteins) * 100,
    2
  )
)

protein_missing <-
  protein_missing %>%
  arrange(desc(Pct_missing))

print(head(protein_missing, 20))

# ==========================
# 6. OVERZICHT HIGH-MISSING
# ==========================

high_missing <- protein_missing %>%
  filter(Pct_missing > 50)

cat("\nAantal eiwitten met >50% missing:\n")
print(nrow(high_missing))

# ==========================
# 7. ANNOTATION FILE CONTROLEREN
# ==========================

annotations %>%
  select(
    protein,
    n_missing,
    pct_missing,
    pct_below_lod,
    duplicate_olink_id
  ) %>%
  arrange(desc(pct_missing)) %>%
  head(20)

# ==========================
# 8. VISUALISATIE MISSINGNESS
# ==========================

ggplot(missing_metadata,
       aes(x = reorder(Variable,
                       Pct_missing),
           y = Pct_missing)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Missingness metadata",
    x = "",
    y = "% missing"
  )

  ggplot(protein_missing,
       aes(x = Pct_missing)) +
  geom_histogram(bins = 30) +
  labs(
    title = "Verdeling missingness eiwitten",
    x = "% missing",
    y = "Aantal eiwitten"
  )
  ggplot(protein_missing,
       aes(x = Pct_missing)) +
  geom_histogram(bins = 30) +
  labs(
    title = "Verdeling missingness eiwitten",
    x = "% missing",
    y = "Aantal eiwitten"
  )
 ggplot(protein_missing,
       aes(x = Pct_missing)) +
  geom_histogram(bins = 5) +
  labs(
    title = "Verdeling missingness eiwitten",
    x = "% missing",
    y = "Aantal eiwitten"
  )
                                      
ggplot(protein_missing,
       aes(x = Pct_missing)) +
  geom_histogram(bins = 10) +
  labs(
    title = "Verdeling missingness eiwitten",
    x = "% missing",
    y = "Aantal eiwitten"
  )
ggplot(protein_missing,
       aes(x = Pct_missing)) +
  geom_histogram(bins = 15) +
  labs(
    title = "Verdeling missingness eiwitten",
    x = "% missing",
    y = "Aantal eiwitten"
  )
#TEST 1234