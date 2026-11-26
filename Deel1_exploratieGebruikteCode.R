library(tidyverse)
meta <- read_csv("sample_metadata.csv")
npx <- read_csv("proteins_npx_raw.csv")
ann <- read_csv("protein_annotation.csv")

#Vraag 1
names(meta)
n_distinct(meta$subject_guid) #Gekeken hoeveel verschillende waarden er zijn voor id en dat was 96 dus 96 vss deelnemers
meta %>% count(subject_guid) #Heeft elke patient even veel stalen?

per_persoon <- meta %>% count(subject_guid)#minimum aantal stalen dat er was afgenomen, mediaan,gemiddelde,(maximum)
summary(per_persoon$n)

glimpse(meta)#welke kolommen zijn er en welk type data

#vraag 2 (lijn 21 kort bespreken)
sum(is.na(npx)) #10.486 ontbrekende NPX waarden. code eronder om te zien hoeveel % er mist
mean(is.na(npx)) #0,8% mist dus dit valt heel goed mee
na_per_eiwit <- colMeans(is.na(npx[,, -1])) * 100#% missend per eiwit berekenen start
sum(na_per_eiwit > 0) #Hoeveel van de eiwitten heeft meer dan 1 waarde die mist = 39
sort(na_per_eiwit, decreasing = TRUE)#Welke eiwitten missen het meest? SPECIAAL: bepaalde eiwitten die zelfde % hebben, dit wijst op technische fout wss
ann %>% filter(protein %in% c("CASC4", "CD97", "DARS"))#3 eiwitten die zelfde % hadden even beter bekijken, code eronder is om 'same assay as' te kunnen zien in spreadsheet
View(ann %>% filter(protein %in% c("CASC4", "CD97", "DARS")))#CASC4 verwijst naar GOLM2, CD97 naar ADGRE5 en DARS naar DARS1. Het zijn twee namen voor hetzelfde eiwit
na_per_eiwit[c("CD97", "ADGRE5")]#Ze zijn op gegeven moment overgeschakeld naar nieuw enaamgeving voor eiwit, dat is waarom waarde voor oude mist
sum(!is.na(npx$CD97) & !is.na(npx$ADGRE5))#Dit is de conformatie want output is 0 wat wilt zeggen dat er geen enkel staal is met namen in beide kolommen (oud en nieuw)
table(meta$olink_batch_id[!is.na(npx$CD97)])#Wanneer werd de oude naam gebruikt
#alle 156 stalen waarin de oude naam cd97 is ingevuld komen uit eenzelfde batch, batch id staat erbij,in deze meetrun gebruikte ze dus nog cd97 als naam terwijl in latere runs ADGRE5 als naam.
#HET IS DUS 1 EIWIT VERDEELD OVER 2 KOLOMMEN, DE 82% MISSENDE WAARDEN MISSEN DUS NIET ECHT ZE HEBBEN GEWOON ANDERE NAAM GEBRUIKT
# Wat gaan we doen met deze paren? We gaan dnek ik best de kolommen samenvoegen tot 1 kolom (niet verwijderen of laten staan)
ann %>% filter(!is.na(same_assay_as)) %>% select(protein, same_assay_as)#Dit zijn de paren die dubbel zijn  

#START code om paren samen te voegen
# 1. per paar de oude naam (> 50 % missend) en de nieuwe naam. HIER HOUDEN WE ENKEL DE KOLOMMEN MET NAMEN WAAR DE OUDE NAAM DE EERSTE IS, moet gedaan worden want hiervoor zagen we 18 rijen/paren maar staan er in beide richtingen in, er zijn er eigenlijk dus maar 9
paren <- ann %>%
  filter(!is.na(same_assay_as)) %>%
  select(oud = protein, nieuw = same_assay_as) %>%
  filter(na_per_eiwit[oud] > 50)

paren   # moet 9 rijen tonen

# 2. elk paar samenvoegen onder de nieuwe naam
for (i in 1:nrow(paren)) {
  npx[[paren$nieuw[i]]] <- coalesce(npx[[paren$nieuw[i]]], npx[[paren$oud[i]]])
}

# 3. de oude kolommen verwijderen
npx <- npx %>% select(-all_of(paren$oud))

ncol(npx) - 1   # aantal eiwitten na samenvoegen. RESULTAAT 1472, ZOALS MEN ZEGT IN HET ARTIKEL DAT MEN DIT AANTAL PROTEINS GAAT BEKIJKEN
#STOP code samenvoegen paren

na_per_eiwit <- colMeans(is.na(npx[, -1])) * 100
sum(na_per_eiwit > 0)#21. ER ZIJN 21 MISSENDE EIWITTEN VAN DE 1472, GROOTSTE REDEN VOOR MISSENDE WAARDE WAREN DUS DE NAAMWIJZIGINGEN
max(na_per_eiwit)#16.609. MAX 16% MIST BIJ 1 EIWIT,DUS CAVA
sort(na_per_eiwit[na_per_eiwit > 0], decreasing = TRUE)#Kijken welke eiwitten dan het meeste missen na de naamcorrectie
#Opmerkelijk, 15 vd 21 eiwitten hebben zelfde waarde, weer iets systemisch niet juist.
table(meta$olink_batch_id)#Om te zien of er een batch is met 143,922 stalen ongeveer (16,6% van 867)
#Batch Q-04064 heeft 144 stalen
table(meta$olink_batch_id[is.na(npx$CLSTN1)])#controle met 1vd 15 eiwitten
#output=144, dit wil zeggen dat alle ontbrekende waarde voor dat eiwit uit die batch komen

#Er zijn dus 15 eiwitten die in specifieke batch niet geraporteerd werden, WAT DOEN WE HIERMEE? zie foto
#Opties zijn verwijderen, invullen of llaten staan. Claude zegt beste keuze "Behouden, NA laten staan"
#Verantwoording indien akkoord behouden: Fifteen assays were not reported in batch Q-04064 (16.6% of samples); this missingness was batch-specific and therefore technical. These assays were retained; missing values were not imputed for statistical testing, and the assays were excluded from analyses requiring complete data (PCA).

#Nog even checken voor resterende 6 eiwitten die missen of er iets systematisch aan de hand is.
na_rest <- na_per_eiwit[na_per_eiwit > 0 & na_per_eiwit < 16]
sort(na_rest, decreasing = TRUE)
#Opmerkelijk, MPO/SERPINB1/CD177/GRN zijn allemaal eiwitten in d eneutrofielen (is er biologisch dus eerder iets?)
table(meta$olink_batch_id[is.na(npx$SERPINB1)])#Ligt het hier bij deze eiwitten ook weer aan de batch?
sum(is.na(npx$MPO) & is.na(npx$SERPINB1))
sum(is.na(npx$MPO))
#Alle ontbrekende waarden voor MPO komen uit dezelfde batch (20201752) en SERPINB1 ontbreekt ook in die batch, probleem weer aan de batch.
#Eiwit ontbreekt hier niet in hele batch, slechts in een deel. Misschien ligt het aan bepaalde plaat in de batch
table(meta$olink_plate_id[is.na(npx$MPO)])

#CONCLUSIEZIN VOOR IN RAPPORT EIWITTEN missende waardes: Missing NPX values (0.8% overall) were confined to 21 assays and were batch- or plate-specific (15 assays absent in batch Q-04064; 6 assays missing on Neurology/Oncology plates of batch 20201752), indicating a technical origin. These assays were retained without imputation for statistical testing and excluded from analyses requiring complete data


#Nu verder kijken naar metadata.
#Hier zijn metingen niet technisch, ontbrekende data kan gewoon doordat iets ni is gemeten tijdens bezoek

na_meta <- colMeans(is.na(meta)) * 100
sort(na_meta, decreasing = TRUE)[1:15]

table(meta$sample.visitName[!is.na(meta$cmv.igg_serology)]) #Bij welke bezoeken is Serologie ingevuld?
#Vooral gemeten bij flu year 1 dag 7, voor de rest maar af en toe dus FOUT LIGT BIJ DESIGN NIET METING
sum(is.na(meta$cmv_status)) #output=0, hieruit besluiten we dat het niet erg is dat cmv.igg_serology mist want we hebbe alle info die we nodig hebben in cmv_status
#BEST DUS OM TE BESLISSEN GEEN CMV.IGG SEROLOGY TE GEBRUIKEN

table(meta$sample.visitName[is.na(meta$am.bmi)])
table(meta$sample.visitName)
#BMI ontbreekt bij bijna elk onderzoek, niet enkel bij serologie. 
#MOGELIJKE VERKLARING:BMI werd niet bij elk bezoek opnieuw gemeten, en later in de studie gebeurde dat consequenter. Het missen hangt dus samen met het tijdstip in de studie, niet met een technische fout.
sum(meta$am.bmi == 0, na.rm = TRUE)
#Er zijn 2 stalen met BMI=0, kan niet, we moeten deze op NA zetten
meta <- meta %>% mutate(am.bmi = na_if(am.bmi, 0))#Code om 2 waarden naar NA om te zetten

#Indien we verder gaan werken met BMI in deel 2 zie foto voor opties wat te doen met missende BMI waarden

#Nu eens checken voor egfr
table(meta$sample.drawYear[is.na(meta$chem.egfr_aa)])
table(meta$sample.drawYear)
#Je ziet bij 2021 en 2022 dat aantal stalen dat waarde ontbreekt= aantal stalen.Labo is dus getsopt met dit te rapporteren
#Best dus om egfr niet te gebruiken,Wil je de nierfunctie toch meenemen, gebruik dan chem.creatinine. eGFR wordt daar oorspronkelijk zelf uit berekend.
sum(is.na(meta$chem.creatinine)) # Eens kijken hoe vaak ontbreekt creatinine.output=10 dus kunnen we wel gebruiken

#Nu chem protein checken
table(meta$sample.drawYear[is.na(meta$chem.protein)])
table(meta$sample.visitName[is.na(meta$chem.protein)])
table(meta$sample.drawYear)
#Hier zien we dat pas vanaf 2022 100% ontbreekt. Men is gewoon gestopt met dit te meten net zoals egfr
#Deze dus ook niet gebruiken in deel 2


#Vraag 3
summary(ann$pct_below_lod)
hist(ann$pct_below_lod, breaks = 20, xlab = "% stalen onder LOD", main = "") #Bekijken of er veel eiwitten onder LOD liggen, is niet zo
#Eiwitten die onder LOD liggen, wegfilteren. Waar leggen we grens?
sum(ann$pct_below_lod >= 25)
sum(ann$pct_below_lod >= 50)
sum(ann$pct_below_lod >= 75)
#Ik zou 50% als grens kiezen want dit wilt zeggen dat bij 50% van keren het eiwit onder LOD ligt
#Verantwoording in rapport: "Assays with ≥ 50% of samples below the LOD were excluded, as the majority of their values cannot be distinguished from background noise. This criterion is independent of the group labels."
#Hoeveel eiwitten houden we dan nog over zie code erronder
keep <- ann %>%
  filter(protein %in% names(npx)) %>%   # alleen eiwitten die nog in npx staan
  filter(pct_below_lod < 50) %>%        # minder dan 50 % onder LOD
  pull(protein)                         # haal de namen eruit als lijstje

length(keep)
#We houden nog 1322 eiwitten over, dus er zijn 150 eiwitten verwijderd (oorspronkelijk 1472 eiw)

summary(meta$n_missing_npx) #Checken of er stalen zijn met slechte meetkwaliteit
summary(meta$pct_below_lod)
hist(meta$pct_below_lod, breaks = 30, xlab = "% eiwitten onder LOD per staal", main = "") #Alles is eig oke

mediaan_staal <- apply(npx[, keep], 1, median, na.rm = TRUE)
summary(mediaan_staal)
#We checken nu of er voor bepaald staal algemeen veel te hoge/lage waarden worden gemeten. Aangezien het min op -0,4755 ligt moeten we wel ff verder kijken. elk staal wordt weergegeven als 1 getal namelijk de mediaan van de waarde van alle eiwitten tesamen
boxplot(mediaan_staal ~ meta$olink_batch_id, xlab = "Batch", ylab = "Mediane NPX per staal") #checken of de afwijkende stalen samen liggen of verspreid (batch effect)
#We zien dat laatste batch alles lager ligt (batch Q-04064)
#Rekening hiermee houden in verdere analyses want is batch effect


#Vraag 4 preprocessing
#log transformatie
hist(npx$GDF15, breaks = 30, xlab = "NPX GDF15", main = "") #Random eiwit een keer log transformeren en verdeling bekijken
#CONCLUSIE:NPX is al een log2-schaal, en de verdeling is vrij symmetrisch. Geen log-transformatie nodig voor de eiwitten.
#LET OP: sommige klinische waarden wel heel scheef verdeeld dus rekening mee houden indien we zo een opnemen in model

#normalisatie
#2 mogelijkheden voor corrigeren voor batch effect zie foto
table(meta$olink_batch_id, meta$subject.ageGroup)#nagaan of optie 1 foto gevaarlijk is
# in sommige batches veel meer jonge dan ouder en vice versa --> batch en leeftijd zijn verweven= confounding
#CONCLUSIE: In de statistische analyses (deel 2–3) nemen we batch op als variabele in het model. Het model scheidt dan het leeftijdseffect van het batch-effect.
#CONCLUSIE 2:Een vooraf gecorrigeerde versie van de data gebruiken we alleen voor figuren (PCA, heatmap), waar je geen model hebt. Daar komen we bij vraag 5 op terug.
#andere conclusies processing zie rapport pdf

#vraag 5: clusteren de data volgens biologie of technologische factor
pca_eiwitten <- keep[colSums(is.na(npx[, keep])) == 0]   # alleen eiwitten zonder NA, PCA aanvaardt geen ontbrekende waarden. Bij vraag 2 besloten we de 21 eiwitten met NA weg te laten voor PCA.
length(pca_eiwitten)

mat <- as.matrix(npx[, pca_eiwitten])                     # tabel -> numerieke matrix

pca <- prcomp(mat, center = TRUE, scale. = TRUE)#De pca berekenen

var_expl <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1) #hoeveel variatie heeft elke pc
var_expl[1:5]

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     batch = meta$olink_batch_id,
                     leeftijd = meta$subject.ageGroup)

ggplot(pca_df, aes(PC1, PC2, color = batch)) +
  geom_point(size = 1.5) +
  theme_classic() +
  labs(x = paste0("PC1 (", var_expl[1], " %)"),
       y = paste0("PC2 (", var_expl[2], " %)"))
#Volgende stap is figuur kleuren volgens leeftijd
ggplot(pca_df, aes(PC1, PC2, color = leeftijd)) +
  geom_point(size = 1.5) +
  theme_classic() +
  labs(x = paste0("PC1 (", var_expl[1], " %)"),
       y = paste0("PC2 (", var_expl[2], " %)"))
#Conclusie is dus dat batch de stalen beter scheidt dan leeftijd
#conclusie De grootste bron van verschillen tussen de stalen is technisch (batch), niet biologisch. Leeftijdsgroep geeft een zwakkere scheiding binnen de batches. Dat bevestigt de beslissing van vraag 4: batch moet in elk model worden meegenomen


#Vraag 6: is de opgekuiste data klaar voor volgend deel
base <- meta %>% filter(sample.daysSinceFirstVisit == 0)
nrow(base)
n_distinct(base$subject_guid)

table(base$olink_batch_id, base$subject.ageGroup)
sum(is.na(base$am.bmi))
