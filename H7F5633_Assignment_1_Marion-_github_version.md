Bioinformatics Analysis and Visualisation of Medical Genomics Data -
H7F5633 - Assignment 1
================
Marion
2024-09-19

# Task 1 - Literature

## 1.

Group 6 –\> Publication: “Single-cell analyses and host genetics
highlight the role of innate immune cells in COVID-19 severity” Edahiro
et al., Nat Genet, 2023

## 2.

### a. Medically relevant insight from the article:

Monocytes or dendritic cells are dysfunctional in patients with severe
COVID-19. Of note, the monocyte dysfunction is associated with a
particular variant of the gene IFNAR2, encoding the receptor for IFN-I,
in monocytes. In addition, there is an involvement of host genetic risk
of severe COVID-19 with innate immunity.

### b. Genomics technology that were used:

- 10 X chromium was used for the single-cell transcriptional profiling
  (scRNAseq)
- 10 X chromium was also used for the sequencing of T cell receptors
  (TCR) and B cell receptors (BCR) to analyze the T cell and B cell
  repertoires
- SNP genotyping was used to analyze the host genetics data

## 3.

### a. Questions that extend the analysis presented in the paper:

- What are the mechanisms involved in the reduced differentiation from
  classical to non-classical monocytes (ncMono) in patients with severe
  COVID-19?
- What are the mechanisms involved in the dysfunction of monocytes in
  patients with severe COVID-19?
- Could the reduced interaction between CXCL10 and CXCR3 be one of the
  factors responsible for driving COVID-19 severity?

# Task 2 - Git repositories and R Markdown

–\> creation of a GitHub account

# Task 3

``` r
# load the tidyverse and ggplot2 packages that are necessary for subsequent tasks
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggplot2)
```

Note: I cannot install (and load) the Bioconductor package because the
latest version of Bioconductor is incompatible with the R version that I
am currently using (4.4.1). To use it, I would need to downgrade to the
4.4.0 R version…

# Task 4 - Using R example data sets

## 1 and 2. Using the function help() to look at the content of the “CO2” dataset.

``` r
help("CO2")
```

Here is the description of the R internal dataset “CO2” (Carbon Dioxide
Uptake in Grass Plants): “The CO2 data frame has 84 rows and 5 columns
of data from an experiment on the cold tolerance of the grass species
Echinochloa crus-galli.”. It is an object of class c(“nfnGroupedData”,
“nfGroupedData”, “groupedData”, “data.frame”).

## 3. Code to calculate the average and median CO2 uptake of the plants from Quebec and Mississippi:

``` r
# Creates a vector which contains the characters "Quebec" and "Mississippi" by taking the levels of the factor "Type" (one of the columns of the CO2 dataset)
regions <- levels(CO2$Type)
# Creates a loop so that the actions will be done for each of the characters in the vector "regions"
for (i in regions) {
  # calculates the mean of the values in the column "uptake" from the CO2 dataset, only for the rows which have "i" as a "Type"
  average <- mean(CO2$uptake[CO2$Type==i])
  # prints "average CO2 uptake of the plants from" + i + "=" + the calculated average
  print (paste("average CO2 uptake of the plants from", i, "=", average))
  # calculates the median of the values in the column "uptake" from the CO2 dataset, only for the rows which have "i" as a "Type"
  median <- median(CO2$uptake[CO2$Type==i])
  # prints "median CO2 uptake of the plants from" + i + "=" + the calculated median
  print (paste("median CO2 uptake of the plants from", i, "=", median))
}
```

    ## [1] "average CO2 uptake of the plants from Quebec = 33.5428571428571"
    ## [1] "median CO2 uptake of the plants from Quebec = 37.15"
    ## [1] "average CO2 uptake of the plants from Mississippi = 20.8833333333333"
    ## [1] "median CO2 uptake of the plants from Mississippi = 19.3"

# Task 5 - R functions

# 1. Code to create a function that calculates the ratio of the mean and the median of a given vector

``` r
# Creates a function ratio() that calculates the ratio of the mean and the median for "vector", using the function function()
ratio <- function(vector) {
  # Calculates the mean of the vector
  x <- mean(vector)
  # Calculates the median of the vector
  y <- median(vector)
  # calculates the ratio of the mean and median and returns this as a result of the function ratio()
  return(x/y)
}
# Creates a vector "a" to test the function ratio()
a <- c(1, 3, 7, 11)
# test the function ratio
ratio(a)
```

    ## [1] 1.1

## 2. Code to create a function that ignores the lowest and highest values from a given a vector and calculate the mean

``` r
# Creates a function mean_vector() that calculates the mean of a given vector after removing the highest and lowest values of this vector
mean_vector <- function(vector) {
  # Finds the highest value from the vector
  max <- max(vector)
  # Finds the lowest value from the vector
  min <- min(vector)
  # Creates a filtered vector by removing the highest and lowest values from the initial vector
  filtered_vector <- vector[!vector %in% c(min, max)]
  # calculates the mean of the filtered vector and returns this as a result of the function mean_vector()
  return(mean(filtered_vector))
}
# Creates a vector "a" to test the function ratio()
a <- c(1, 3, 7, 11)
# test the function mean_vector()
mean_vector(a)
```

    ## [1] 5

## 3. Piping

<https://r4ds.had.co.nz/pipes.html#pipes>

Why, how and when not to use pipes: Pipes simplify complex data
manipulations by chaining operations, enhancing code readability and
reducing intermediate object creation. The %\>% operator allows for
intuitive left-to-right function composition. Pipes are ideal for linear
sequences but less suitable for complex dependency structures or when
dealing with multiple inputs/outputs. Pipes are are not recommended for
pipes longer than ten steps (it is better to use intermediate objects in
that case). And finally, they don’t work with functions that use the
current environment or rely on lazy evaluation.

## 4. Apply-family of functions

<http://uc-r.github.io/apply_family>

Why the apply-family of functions can be useful in my work: The
apply-family of functions, including apply(), lapply(), sapply(), and
tapply(), apply a given function to a data object and can replace loops
(in certain cases), which allows for a code that is more compact and
easier to read. I could use this family of functions for both
statistical analysis and visualization of my datasets, particularly the
apply() function, to apply a function to each row or column of my
matrices. It would improve my codes a lot and avoid potential errors in
the analysis.

# Task 6 - Basic visualization with R

## 1. Code to compare the distributions of the body heights of the two species from the “magic_guys.csv” dataset graphically

### a. Histograms

``` r
# read the csv file "magic_guys.csv"
magic_guys <- read_csv("~/Desktop/SciLifeLab_Riken/magic_guys.csv")
```

    ## Rows: 100 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (2): uniqId, species
    ## dbl (2): length, weight
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
body_heights <- magic_guys$length
species <-  unique(magic_guys$species)

# Creates a vector that assign the groups (species) to a color
colors <- setNames(c("sandybrown", "cadetblue"), species)

# Creates a vector containing the values of the x axis where the breaks happen
# It takes the maximum and minimum values of the "body_heights" and requires 10 breaks
my_breaks <- seq(min(body_heights), max(body_heights), length.out = 10)

# Creates an empty plot
plot(NA, xlim = range(body_heights), ylim = c(0, 20), 
     xlab = "Body height (cm)", ylab = "Frequency", 
     main = "Distributions of the body heights in different species - Base R")

# Loop to create an histogram (with 30 breaks) for each group (species) and add it to the empty plot
for (i in species) {
  current_data <- magic_guys %>% 
    filter(species == i)
  body_heights <- current_data$length
  hist(body_heights, 
       col = colors[i], 
       border = "black", 
       add = TRUE, 
       breaks = my_breaks)
}

# Adds a legend to the plot
legend("topright", legend = unique(species), fill = colors, title = "Species")
```

![](H7F5633_Assignment_1_Marion-_github_version_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# Same histogram but with 40 breaks instead of 10
body_heights <- magic_guys$length
species <-  unique(magic_guys$species)
colors <- setNames(c("sandybrown", "cadetblue"), species)
my_breaks <- seq(min(body_heights), max(body_heights), length.out = 40)
plot(NA, xlim = range(body_heights), ylim = c(0, 7), 
     xlab = "Body height (cm)", ylab = "Frequency", 
     main = "Distributions of the body heights in different species - Base R")
for (i in species) {
  current_data <- magic_guys %>% 
    filter(species == i)
  body_heights <- current_data$length
  hist(body_heights, 
       col = colors[i], 
       border = "black", 
       add = TRUE, 
       breaks = my_breaks)
}
legend("topright", legend = unique(species), fill = colors, title = "Species")
```

![](H7F5633_Assignment_1_Marion-_github_version_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
# The problem with these plots is that there is no transparency so we can't see the overlap...
```

#### Using ggplot()

``` r
body_heights <- magic_guys$length
species <-  unique(magic_guys$species)

# creates a histogram with the two groups in different colors, with 40 breaks (bins)
# The colors of the groups and the titles (plot, axis and legend) are customized
# The bars are semi-transparents to see a potential overlap of the two histograms
histo_gg <- ggplot(magic_guys, aes(x = body_heights, fill = as.factor(species))) +
# the argument "identity" og geom_histogram() allows to have two separate histograms instead of having only one and the argument apha allows to change the transparency so that we can see the overlap (instead of one histogram masking the other)
geom_histogram(color = "black", bins = 40, alpha = 0.5, position = "identity") +
scale_fill_manual(values = c("sandybrown", "cadetblue"), name = "Species") +
theme_classic() +
labs(title = "Distributions of the body heights in different species - Gggplot", x = "Body height (cm)", y = "Frequency") 

print(histo_gg)
```

![](H7F5633_Assignment_1_Marion-_github_version_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### b. Boxplots

``` r
body_heights <- magic_guys$length
species <-  unique(magic_guys$species)

box_gg <- ggplot(magic_guys, aes(x = species, y = body_heights, fill = as.factor(species))) +
  geom_boxplot(width = 0.5, color = "black") +
  scale_fill_manual(values = c("sandybrown", "cadetblue"), name = "Species") +
  theme_classic() +
  labs(title = "Distributions of the body heights in different species - Boxplot", x = NULL, y = "Body height (cm)") 

print(box_gg)
```

![](H7F5633_Assignment_1_Marion-_github_version_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### c. Saving plots in various formats

``` r
# choosing the working directory where the files will be saved
setwd("~/Desktop/SciLifeLab_Riken/")

# opens a png file and names it
png("boxplot.png")

# print the boxplot
print(box_gg)

# close the file
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# Same with pdf file
pdf("boxplot.pdf")
print(box_gg)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# Same with svg file
svg("boxplot.svg")
```

    ## Warning in grSoftVersion(): unable to load shared object '/Library/Frameworks/R.framework/Resources/modules//R_X11.so':
    ##   dlopen(/Library/Frameworks/R.framework/Resources/modules//R_X11.so, 0x0006): Library not loaded: /opt/X11/lib/libSM.6.dylib
    ##   Referenced from: <31EADEB5-0A17-3546-9944-9B3747071FE8> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/modules/R_X11.so
    ##   Reason: tried: '/opt/X11/lib/libSM.6.dylib' (no such file), '/System/Volumes/Preboot/Cryptexes/OS/opt/X11/lib/libSM.6.dylib' (no such file), '/opt/X11/lib/libSM.6.dylib' (no such file), '/Library/Frameworks/R.framework/Resources/lib/libSM.6.dylib' (no such file), '/Library/Java/JavaVirtualMachines/jdk-11.0.18+10/Contents/Home/lib/server/libSM.6.dylib' (no such file)

    ## Warning in cairoVersion(): unable to load shared object '/Library/Frameworks/R.framework/Resources/library/grDevices/libs//cairo.so':
    ##   dlopen(/Library/Frameworks/R.framework/Resources/library/grDevices/libs//cairo.so, 0x0006): Library not loaded: /opt/X11/lib/libXrender.1.dylib
    ##   Referenced from: <63619C6D-FE72-3544-BCEF-9C834A5E39D8> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/grDevices/libs/cairo.so
    ##   Reason: tried: '/opt/X11/lib/libXrender.1.dylib' (no such file), '/System/Volumes/Preboot/Cryptexes/OS/opt/X11/lib/libXrender.1.dylib' (no such file), '/opt/X11/lib/libXrender.1.dylib' (no such file), '/Library/Frameworks/R.framework/Resources/lib/libXrender.1.dylib' (no such file), '/Library/Java/JavaVirtualMachines/jdk-11.0.18+10/Contents/Home/lib/server/libXrender.1.dylib' (no such file)

    ## Warning in svg("boxplot.svg"): failed to load cairo DLL

``` r
print(box_gg)
```

![](H7F5633_Assignment_1_Marion-_github_version_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
dev.off()
```

    ## null device 
    ##           1

``` r
# Other attempts to save a svg file... still doesn't work
ggsave("boxplot.svg", plot = box_gg)
```

    ## Saving 7 x 7 in image

``` r
# Other attempts to save a svg file... Now loading the "svglite" package (after having installed it on the console) o) as it seems it was missing
library(svglite)
ggsave("boxplot.svg", plot = box_gg)
```

    ## Saving 7 x 7 in image

Notes: - I got multiple errors before I managed to save as “.svg” file.
I solved this by installing and loading the package “svglite” - I get a
problem when “knitting” if I write “install.packages(”svglite”) in the R
chunk, so I have to write it only in the console.

I would use .png for presentations (import in ppt), .svg for further
editing (import in illustrator files) and .pdf

## 2. Work with “microarray_data.tab”

### a.

``` r
# Reads the microarray file (specify the separator as "\ t")
microarray <- read.table("~/Desktop/SciLifeLab_Riken/microarray_data.tab", header=TRUE, sep="\t")
# Calculates number of columns
cols <- ncol(microarray)
# Calculates number of rows
rows <- nrow(microarray)
print(paste("The", "matrix", "has", cols, "columns", "and", rows, "rows"))
```

    ## [1] "The matrix has 1000 columns and 553 rows"

### b.

``` r
# Uses is.na() function to identify the missing values and uses colSums() to calculate this for each column
na_count <- colSums(is.na(microarray))
# print the 10 first values of "na_count"
print(head(na_count, 10))
```

    ##  g1  g2  g3  g4  g5  g6  g7  g8  g9 g10 
    ## 130 104  74  93  81  30  31  26  12  56

``` r
hist(na_count, 
     border = "black", 
     breaks = cols,
     xlim = range(na_count), 
     ylim = c(0, 25),
     xlab = "number of missing values per gene", 
     ylab = "Frequency", 
     main = "Distributions of the number of missing values per gene"
     )
```

![](H7F5633_Assignment_1_Marion-_github_version_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

### c.

``` r
# creates a new function to calculate the percentage of

na_count_df <- as.data.frame(na_count)

calcul_percent <- function(x) {
        return(x/rows*100)
    }

# Apply this new function of all elemnts of the vector na_count
na_count_df$percent <- sapply(na_count_df$na_count, function(x) calcul_percent(x))

na_count_df$gene_name <- row.names(na_count_df)

# creates a new function to detect which genes have more than 50% or 20% or 10% missing values
missing_val <- function(w, x, y, z) {
  if (x > y) {
    z <- c(z, w)
  }
  return(z)
}

# Creates a vector with the different percentages of missing values to be tested (10, 20 and 50%)
cut_off <- c(10, 20, 50)  
genes <- na_count_df$gene_name

# Creates a loop to find the genes with missing values for the cut-offs contained in the vector "cut_off"
for (y in cut_off) {
  # Creates an empty vector, which will later be filled with the names of genes
  z <- c()
  # Creates a loop so that the empty vector is filled for the genes with a % value greater than y
  for (w in genes) {
    x <- na_count_df$percent[na_count_df$gene_name == w]
    z <- missing_val(w, x, y, z)
  }
  # If the vector containing the genes with a % value greater than y are longer than 0 (meaning there is at least 1 gene that is contained in the vector), it will print a sentence along with the names of those genes
  if (length(z) > 0) {
    print(paste("Genes with greater than", y, "% missing values:", paste(z, collapse=", ")))
  } else {
    print(paste("No genes with greater than", y, "% missing values"))
  }
}
```

    ## [1] "Genes with greater than 10 % missing values: g1, g2, g3, g4, g5, g10, g11, g12, g14, g15, g16, g18, g21, g22, g23, g24, g25, g26, g28, g29, g35, g36, g37, g38, g39, g40, g41, g42, g44, g45, g46, g47, g48, g49, g50, g51, g52, g53, g54, g55, g56, g57, g58, g59, g60, g61, g62, g63, g64, g65, g66, g67, g68, g69, g70, g71, g72, g73, g74, g75, g76, g77, g78, g79, g80, g81, g82, g83, g84, g85, g86, g87, g88, g89, g90, g91, g92, g93, g94, g95, g96, g97, g98, g99, g100, g101, g102, g103, g104, g105, g106, g107, g108, g109, g110, g111, g112, g113, g114, g115, g116, g117, g118, g119, g120, g130, g131, g132, g133, g134, g135, g136, g137, g138, g139, g140, g142, g147, g148, g151, g152, g153, g154, g155, g156, g157, g158, g159, g160, g165, g171, g172, g173, g174, g175, g176, g177, g178, g179, g180, g194, g196, g200, g204, g210, g221, g223, g233, g239, g241, g242, g243, g244, g252, g258, g260, g263, g264, g268, g274, g280, g281, g284, g285, g286, g287, g288, g290, g292, g294, g295, g296, g297, g298, g299, g301, g309, g311, g312, g313, g314, g315, g316, g320, g321, g322, g323, g324, g327, g329, g331, g332, g333, g334, g335, g336, g337, g339, g344, g347, g348, g351, g352, g353, g354, g355, g356, g357, g358, g359, g360, g361, g362, g363, g364, g365, g366, g367, g368, g369, g370, g372, g374, g377, g378, g379, g380, g381, g382, g383, g384, g385, g386, g387, g388, g389, g390, g391, g392, g396, g400, g401, g402, g403, g404, g405, g406, g407, g408, g409, g410, g411, g413, g415, g416, g417, g418, g419, g421, g423, g425, g429, g430, g431, g432, g433, g434, g435, g436, g437, g438, g439, g440, g445, g450, g453, g455, g459, g460, g461, g462, g463, g464, g465, g466, g467, g468, g469, g470, g476, g477, g479, g481, g483, g491, g492, g493, g494, g495, g496, g497, g498, g499, g500, g502, g507, g510, g513, g514, g515, g518, g519, g520, g522, g524, g525, g527, g530, g531, g532, g533, g534, g535, g536, g537, g538, g539, g540, g541, g544, g547, g553, g555, g556, g558, g559, g561, g563, g564, g567, g569, g570, g571, g572, g573, g574, g575, g576, g577, g578, g579, g580, g583, g585, g586, g589, g591, g592, g595, g597, g599, g607, g610, g611, g612, g613, g614, g615, g616, g617, g618, g619, g620, g631, g638, g650, g653, g657, g660, g663, g666, g669, g672, g681, g689, g691, g694, g696, g700, g707, g709, g711, g715, g718, g719, g722, g724, g726, g743, g744, g747, g748, g749, g751, g752, g753, g754, g755, g756, g757, g758, g759, g760, g761, g762, g763, g764, g765, g766, g767, g768, g769, g770, g776, g781, g782, g783, g784, g785, g786, g787, g788, g789, g790, g795, g796, g801, g802, g803, g804, g805, g806, g807, g808, g809, g810, g812, g814, g818, g820, g822, g824, g831, g832, g833, g834, g835, g836, g837, g838, g839, g840, g843, g849, g850, g851, g854, g861, g862, g863, g864, g865, g866, g867, g868, g869, g870, g872, g874, g882, g884, g891, g892, g893, g894, g895, g896, g897, g898, g899, g900, g903, g909, g910, g911, g919, g922, g926, g931, g932, g933, g934, g935, g936, g937, g938, g939, g940, g945, g947, g948, g951, g952, g956, g963, g965, g970, g971, g972, g973, g974, g975, g976, g977, g978, g979, g980, g985, g989"
    ## [1] "Genes with greater than 20 % missing values: g1, g14, g18, g21, g22, g24, g25, g26, g29, g39, g40, g41, g48, g51, g52, g53, g54, g55, g56, g57, g58, g59, g60, g61, g62, g63, g64, g65, g66, g67, g68, g69, g70, g71, g72, g73, g74, g75, g76, g77, g78, g79, g80, g81, g82, g83, g84, g85, g86, g87, g88, g89, g90, g91, g92, g93, g94, g95, g96, g97, g98, g99, g100, g101, g102, g103, g104, g105, g106, g107, g108, g109, g110, g111, g112, g113, g114, g115, g116, g117, g118, g119, g120, g130, g131, g132, g133, g134, g135, g136, g137, g138, g139, g140, g142, g147, g148, g151, g152, g153, g154, g155, g156, g157, g158, g159, g160, g165, g171, g172, g173, g174, g175, g176, g177, g178, g179, g180, g196, g200, g204, g210, g233, g252, g260, g290, g297, g301, g321, g329, g332, g333, g334, g335, g344, g351, g352, g353, g354, g355, g356, g357, g358, g359, g360, g361, g362, g363, g364, g365, g366, g367, g368, g369, g370, g379, g381, g382, g383, g384, g385, g386, g387, g388, g389, g390, g391, g396, g400, g401, g402, g403, g404, g405, g406, g407, g408, g409, g410, g415, g417, g418, g431, g432, g433, g434, g435, g436, g437, g438, g439, g440, g445, g450, g455, g461, g462, g463, g464, g465, g466, g467, g468, g469, g470, g476, g491, g492, g493, g494, g495, g496, g497, g498, g499, g500, g510, g513, g515, g518, g519, g520, g522, g527, g531, g532, g533, g534, g535, g536, g537, g538, g539, g540, g544, g547, g558, g561, g569, g570, g571, g572, g573, g574, g575, g576, g577, g578, g579, g580, g583, g585, g586, g591, g599, g611, g612, g613, g614, g615, g616, g617, g618, g619, g620, g650, g657, g663, g669, g689, g691, g694, g744, g751, g752, g753, g754, g755, g756, g757, g758, g759, g760, g761, g762, g763, g764, g765, g766, g767, g768, g769, g770, g781, g782, g783, g784, g785, g786, g787, g788, g789, g790, g801, g802, g803, g804, g805, g806, g807, g808, g809, g810, g814, g831, g832, g833, g834, g835, g836, g837, g838, g839, g840, g849, g850, g851, g854, g861, g862, g863, g864, g865, g866, g867, g868, g869, g870, g872, g891, g892, g893, g894, g895, g896, g897, g898, g899, g900, g910, g919, g926, g931, g932, g933, g934, g935, g936, g937, g938, g939, g940, g947, g951, g970, g971, g972, g973, g974, g975, g976, g977, g978, g979, g980, g985, g989"
    ## [1] "Genes with greater than 50 % missing values: g18, g48, g55, g58, g60, g66, g73, g79, g83, g91, g93, g94, g99, g105, g109, g132, g135, g137, g138, g172, g260, g290, g301, g329, g352, g355, g362, g363, g368, g383, g388, g389, g390, g391, g406, g417, g431, g432, g440, g461, g462, g498, g519, g527, g531, g532, g538, g572, g575, g576, g577, g585, g615, g619, g663, g669, g751, g753, g766, g768, g788, g802, g804, g838, g851, g854, g864, g892, g893, g898, g919, g932, g971, g980"

### d.

``` r
# Calculate the mean of each columns of "microarray", while ignoring NAs
mean_result <- colMeans(microarray, na.rm = TRUE)

# create a vector of gene names by taking the column names of "micrpoarray"
microarray_genes <- colnames(microarray)

#print "microarray" before changing the NA values
print(head(microarray[, 1:10], 10))
```

    ##        g1     g2     g3     g4     g5     g6     g7     g8     g9    g10
    ## 1   1.802     NA -0.182  1.312  3.497  0.439  0.777  0.379  0.336  0.877
    ## 2      NA     NA  7.693     NA  0.193 -1.383 -1.309 -0.424 -0.270 -0.519
    ## 3   1.079     NA  1.556  1.652     NA  0.460  0.715  0.375 -0.138 -0.261
    ## 4   3.607     NA  1.914     NA  1.400  1.109  2.143  1.571 -0.271 -0.309
    ## 5  -1.700     NA  0.943     NA -0.170     NA -0.041     NA -0.069  1.533
    ## 6      NA     NA  0.043     NA  0.729 -0.089  0.209 -0.308  0.474 -0.655
    ## 7      NA     NA     NA     NA     NA     NA     NA  2.026 -0.567 -3.088
    ## 8      NA     NA     NA     NA     NA -1.297 -0.506 -0.228  0.595 -0.063
    ## 9      NA  1.831     NA  1.596     NA  2.656  4.831  0.677 -1.678  0.016
    ## 10     NA -2.340     NA -1.833 -2.095 -2.083 -2.977 -1.140 -1.375 -0.218

``` r
# replace NA with the average expression value for the particular gene using a loop (for each gene)
for (gene in microarray_genes) {
  microarray[is.na(microarray[,gene]), gene] <- mean_result[gene]
}

#print "microarray" after this change
print(head(microarray[, 1:10], 10))
```

    ##             g1         g2         g3          g4          g5         g6      g7
    ## 1   1.80200000  0.1656927 -0.1820000  1.31200000  3.49700000  0.4390000  0.7770
    ## 2   0.02547518  0.1656927  7.6930000 -0.06731957  0.19300000 -1.3830000 -1.3090
    ## 3   1.07900000  0.1656927  1.5560000  1.65200000 -0.01812288  0.4600000  0.7150
    ## 4   3.60700000  0.1656927  1.9140000 -0.06731957  1.40000000  1.1090000  2.1430
    ## 5  -1.70000000  0.1656927  0.9430000 -0.06731957 -0.17000000 -0.1571338 -0.0410
    ## 6   0.02547518  0.1656927  0.0430000 -0.06731957  0.72900000 -0.0890000  0.2090
    ## 7   0.02547518  0.1656927 -0.1230605 -0.06731957 -0.01812288 -0.1571338 -0.2625
    ## 8   0.02547518  0.1656927 -0.1230605 -0.06731957 -0.01812288 -1.2970000 -0.5060
    ## 9   0.02547518  1.8310000 -0.1230605  1.59600000 -0.01812288  2.6560000  4.8310
    ## 10  0.02547518 -2.3400000 -0.1230605 -1.83300000 -2.09500000 -2.0830000 -2.9770
    ##            g8     g9    g10
    ## 1   0.3790000  0.336  0.877
    ## 2  -0.4240000 -0.270 -0.519
    ## 3   0.3750000 -0.138 -0.261
    ## 4   1.5710000 -0.271 -0.309
    ## 5  -0.1771063 -0.069  1.533
    ## 6  -0.3080000  0.474 -0.655
    ## 7   2.0260000 -0.567 -3.088
    ## 8  -0.2280000  0.595 -0.063
    ## 9   0.6770000 -1.678  0.016
    ## 10 -1.1400000 -1.375 -0.218

## 3. Visualize the data in CO2 dataset in a way that gives you a deeper understanding of the data

``` r
# Plot the CO2 uptake as a function of the concentration, with different colors for "Treatment" and different line type for the "Type"
ggplot(CO2, aes(x = conc, y = uptake, color = Treatment, linetype = Type)) +
  geom_point() +
  # add a trend line
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "CO2 Concentration", 
       y = "CO2 Uptake",
       title = "CO2 Uptake as a function of the concentration") +
  # Different colors depending on Treatment
  scale_color_manual(values = c("darkgreen", "darkmagenta")) + 
  # Different line types depending on Type
  scale_linetype_manual(values = c("dotted", "solid")) + 
  theme_classic()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](H7F5633_Assignment_1_Marion-_github_version_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

What I see: - CO2 uptake is higher in non-chilled compared to chilled
treatment - CO2 uptake is higher in Quebec compared to Mississippi - CO2
uptake seems to reach a plateau after a certain concentration

# Task 7

## 1. Install Tidybiology package, using: devtools::install_github(“hirscheylab/tidybiology”)

``` r
# loading the "tidybiology" package (after installing it)
library(tidybiology)
# printing the "chromosome" dataset
print(head(chromosome, 5))
```

    ## # A tibble: 5 × 14
    ##   id    length_mm basepairs variations protein_codinggenes pseudo_genes
    ##   <fct>     <dbl>     <dbl>      <dbl>               <int>        <int>
    ## 1 1            85 248956422   12151146                2058         1220
    ## 2 2            83 242193529   12945965                1309         1023
    ## 3 3            67 198295559   10638715                1078          763
    ## 4 4            65 190214555   10165685                 752          727
    ## 5 5            62 181538259    9519995                 876          721
    ## # ℹ 8 more variables: totallongnc_rna <int>, totalsmallnc_rna <int>,
    ## #   mi_rna <int>, r_rna <int>, sn_rna <int>, sno_rna <int>, miscnc_rna <int>,
    ## #   centromereposition_mbp <dbl>

``` r
# printing the "proteins" dataset
print(head(proteins, 5))
```

    ## # A tibble: 5 × 8
    ##   uniprot_id gene_name gene_name_alt   protein_name    protein_name_alt sequence
    ##   <chr>      <chr>     <chr>           <chr>           <chr>            <chr>   
    ## 1 P04217     A1BG      <NA>            "Alpha-1B-glyc… Alpha-1-B glyco… MSMLVVF…
    ## 2 Q9NQ94     A1CF      ACF ASP         "APOBEC1 compl… APOBEC1-stimula… MESNHKS…
    ## 3 P01023     A2M       CPAMD5 FWP007   "Alpha-2-macro… Alpha-2-M) (C3 … MGKNKLL…
    ## 4 A8K2U0     A2ML1     CPAMD9          "Alpha-2-macro… C3 and PZP-like… MWAQLLL…
    ## 5 U3KPV4     A3GALT2   A3GALT2P IGBS3S "Alpha-1,3-gal… EC 2.4.1.87) (I… MALKEGL…
    ## # ℹ 2 more variables: length <dbl>, mass <dbl>

### a.

``` r
# Creates a matrix to be filled with summary data
data_summary <- matrix(NA, nrow = 3, ncol = 3)

# Set the names for the columns and rows of this matrix
colnames(data_summary) <- c("variations", "protein_codinggenes", "mi_rna")
rownames(data_summary) <- c("mean", "median", "maximum")

# Print the empty matrix
print(data_summary)
```

    ##         variations protein_codinggenes mi_rna
    ## mean            NA                  NA     NA
    ## median          NA                  NA     NA
    ## maximum         NA                  NA     NA

``` r
# Create a loop so that the statistics are calculated for the three variables of interest
for (v in colnames(data_summary)) {
  # Use the tidyverse to select the column of the variable of interest from the chromosome dataset
  current_data <- chromosome %>% 
    select(v)  %>% 
    # convert this data into a vector
    unlist()
 
# Fill the matrix with the summary data  
data_summary["mean", v] <- mean(current_data)
data_summary["median", v] <- median(current_data)
data_summary["maximum", v] <- max(current_data)
}
```

    ## Warning: Using an external vector in selections was deprecated in tidyselect 1.1.0.
    ## ℹ Please use `all_of()` or `any_of()` instead.
    ##   # Was:
    ##   data %>% select(v)
    ## 
    ##   # Now:
    ##   data %>% select(all_of(v))
    ## 
    ## See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
print(data_summary)
```

    ##         variations protein_codinggenes    mi_rna
    ## mean       6484572            849.9583  73.16667
    ## median     6172346            836.0000  75.00000
    ## maximum   12945965           2058.0000 134.00000

### b.

``` r
# Creates an histogram using ggplot
distrib <- ggplot(chromosome, aes(x = length_mm)) +  
geom_histogram(color = "deeppink4", fill = "deeppink", bins =50) +
theme_classic() +
labs(title = "Distribution of chromosome size", x = "chromosome size (mm)", y = "Frequency") 

print(distrib)
```

![](H7F5633_Assignment_1_Marion-_github_version_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

### c. 

``` r
# Plot to visualize a potential correlation between the number of protein coding genes and the length of the chromosome
blue <- ggplot(chromosome, aes(x = length_mm, y = protein_codinggenes)) +
  geom_point(size = 1, color = "black") +
  labs(x = "Chromosome length (mm)", 
       y = "Number of protein coding genes",
       title = "Number of protein coding genes as a function of chromosome length") +
  theme_classic() +
  # Add a line to visualize a linear regression
  geom_smooth(method = "lm", se = FALSE, color = "blue", size = 0.5)
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
print(blue)
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](H7F5633_Assignment_1_Marion-_github_version_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
# Plot to visualize a potential correlation between the number of miRNAs and the length of the chromosome
red <- ggplot(chromosome, aes(x = length_mm, y = mi_rna)) +
  geom_point(size = 1, color = "black") +
  labs(x = "Chromosome length (mm)", 
       y = "Number of miRNAs",
       title = "Number of miRNAs as a function of chromosome length") +
  theme_classic() +
  # Add a line to visualize a linear regression
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 0.5)

print(red)
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](H7F5633_Assignment_1_Marion-_github_version_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

Outcome: Both the number of protein coding genes and the number of
miRNAs correlate with the chromosome length

### d. 

``` r
summarized_data <- matrix(NA, nrow = 3, ncol = 2)

# Set the names for the columns and rows of this matrix
colnames(summarized_data) <- c("length", "mass")
rownames(summarized_data) <- c("mean", "median", "maximum")

# Create a loop so that the statistics are calculated for the three variables of interest
for (v in colnames(summarized_data)) {
  # Use the tidyverse to select the column of the variable of interest from the chromosome dataset
  current_data <- proteins %>% 
    select(v)  %>% 
    # convert this data into a vector
    unlist()
 
# Fill the matrix with the summarized data  
summarized_data["mean", v] <- mean(current_data)
summarized_data["median", v] <- median(current_data)
summarized_data["maximum", v] <- max(current_data)
}

print(summarized_data)
```

    ##             length       mass
    ## mean      557.1603   62061.38
    ## median    414.0000   46140.50
    ## maximum 34350.0000 3816030.00

``` r
#
length_vs_mass <- ggplot(proteins, aes(x = length, y = mass)) +
  geom_point(size = 3, color = "navyblue", fill = "mediumturquoise", shape = 21) +
  labs(x = "Protein length", 
       y = "Protein mass",
       title = "Protein mass as a function of protein length") +
  theme_classic() +
  geom_smooth(method = "lm", se = FALSE, color = "navyblue", size = 1)

print(length_vs_mass)
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](H7F5633_Assignment_1_Marion-_github_version_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

### 
