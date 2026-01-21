#' Example datasets for saemvs
#'
#' Simulated datasets used to demonstrate the usage of the `saemvs` package
#' for variable selection in nonlinear mixed-effects models.
#'
#' The `.rda` file contains two data frames:
#' - `df_long`: longitudinal measurements in long format
#' - `df_cov`: covariates in wide format
#'
#' @docType data
#' @name example_data_df
#' @usage data(example_data_df)
NULL

#' Longitudinal measurements dataset
#'
#' Longitudinal measurements in long format for 200 individuals.
#'
#' @docType data
#' @name df_long
#' @format A data frame with columns:
#' \describe{
#'   \item{id}{Individual identifier}
#'   \item{time}{Measurement time}
#'   \item{y}{Observed response}
#' }
NULL

#' Covariates dataset
#'
#' Covariates in wide format, one row per individual.
#'
#' @docType data
#' @name df_cov
#' @format A data frame with 201 columns:
#' \describe{
#'   \item{id}{Individual identifier}
#'   \item{x1}{Numeric covariate 1}
#'   \item{x2}{Numeric covariate 2}
#'   \item{x3}{Numeric covariate 3}
#'   \item{x4}{Numeric covariate 4}
#'   \item{x5}{Numeric covariate 5}
#'   \item{x6}{Numeric covariate 6}
#'   \item{x7}{Numeric covariate 7}
#'   \item{x8}{Numeric covariate 8}
#'   \item{x9}{Numeric covariate 9}
#'   \item{x10}{Numeric covariate 10}
#'   \item{x11}{Numeric covariate 11}
#'   \item{x12}{Numeric covariate 12}
#'   \item{x13}{Numeric covariate 13}
#'   \item{x14}{Numeric covariate 14}
#'   \item{x15}{Numeric covariate 15}
#'   \item{x16}{Numeric covariate 16}
#'   \item{x17}{Numeric covariate 17}
#'   \item{x18}{Numeric covariate 18}
#'   \item{x19}{Numeric covariate 19}
#'   \item{x20}{Numeric covariate 20}
#'   \item{x21}{Numeric covariate 21}
#'   \item{x22}{Numeric covariate 22}
#'   \item{x23}{Numeric covariate 23}
#'   \item{x24}{Numeric covariate 24}
#'   \item{x25}{Numeric covariate 25}
#'   \item{x26}{Numeric covariate 26}
#'   \item{x27}{Numeric covariate 27}
#'   \item{x28}{Numeric covariate 28}
#'   \item{x29}{Numeric covariate 29}
#'   \item{x30}{Numeric covariate 30}
#'   \item{x31}{Numeric covariate 31}
#'   \item{x32}{Numeric covariate 32}
#'   \item{x33}{Numeric covariate 33}
#'   \item{x34}{Numeric covariate 34}
#'   \item{x35}{Numeric covariate 35}
#'   \item{x36}{Numeric covariate 36}
#'   \item{x37}{Numeric covariate 37}
#'   \item{x38}{Numeric covariate 38}
#'   \item{x39}{Numeric covariate 39}
#'   \item{x40}{Numeric covariate 40}
#'   \item{x41}{Numeric covariate 41}
#'   \item{x42}{Numeric covariate 42}
#'   \item{x43}{Numeric covariate 43}
#'   \item{x44}{Numeric covariate 44}
#'   \item{x45}{Numeric covariate 45}
#'   \item{x46}{Numeric covariate 46}
#'   \item{x47}{Numeric covariate 47}
#'   \item{x48}{Numeric covariate 48}
#'   \item{x49}{Numeric covariate 49}
#'   \item{x50}{Numeric covariate 50}
#'   \item{x51}{Numeric covariate 51}
#'   \item{x52}{Numeric covariate 52}
#'   \item{x53}{Numeric covariate 53}
#'   \item{x54}{Numeric covariate 54}
#'   \item{x55}{Numeric covariate 55}
#'   \item{x56}{Numeric covariate 56}
#'   \item{x57}{Numeric covariate 57}
#'   \item{x58}{Numeric covariate 58}
#'   \item{x59}{Numeric covariate 59}
#'   \item{x60}{Numeric covariate 60}
#'   \item{x61}{Numeric covariate 61}
#'   \item{x62}{Numeric covariate 62}
#'   \item{x63}{Numeric covariate 63}
#'   \item{x64}{Numeric covariate 64}
#'   \item{x65}{Numeric covariate 65}
#'   \item{x66}{Numeric covariate 66}
#'   \item{x67}{Numeric covariate 67}
#'   \item{x68}{Numeric covariate 68}
#'   \item{x69}{Numeric covariate 69}
#'   \item{x70}{Numeric covariate 70}
#'   \item{x71}{Numeric covariate 71}
#'   \item{x72}{Numeric covariate 72}
#'   \item{x73}{Numeric covariate 73}
#'   \item{x74}{Numeric covariate 74}
#'   \item{x75}{Numeric covariate 75}
#'   \item{x76}{Numeric covariate 76}
#'   \item{x77}{Numeric covariate 77}
#'   \item{x78}{Numeric covariate 78}
#'   \item{x79}{Numeric covariate 79}
#'   \item{x80}{Numeric covariate 80}
#'   \item{x81}{Numeric covariate 81}
#'   \item{x82}{Numeric covariate 82}
#'   \item{x83}{Numeric covariate 83}
#'   \item{x84}{Numeric covariate 84}
#'   \item{x85}{Numeric covariate 85}
#'   \item{x86}{Numeric covariate 86}
#'   \item{x87}{Numeric covariate 87}
#'   \item{x88}{Numeric covariate 88}
#'   \item{x89}{Numeric covariate 89}
#'   \item{x90}{Numeric covariate 90}
#'   \item{x91}{Numeric covariate 91}
#'   \item{x92}{Numeric covariate 92}
#'   \item{x93}{Numeric covariate 93}
#'   \item{x94}{Numeric covariate 94}
#'   \item{x95}{Numeric covariate 95}
#'   \item{x96}{Numeric covariate 96}
#'   \item{x97}{Numeric covariate 97}
#'   \item{x98}{Numeric covariate 98}
#'   \item{x99}{Numeric covariate 99}
#'   \item{x100}{Numeric covariate 100}
#'   \item{x101}{Numeric covariate 101}
#'   \item{x102}{Numeric covariate 102}
#'   \item{x103}{Numeric covariate 103}
#'   \item{x104}{Numeric covariate 104}
#'   \item{x105}{Numeric covariate 105}
#'   \item{x106}{Numeric covariate 106}
#'   \item{x107}{Numeric covariate 107}
#'   \item{x108}{Numeric covariate 108}
#'   \item{x109}{Numeric covariate 109}
#'   \item{x110}{Numeric covariate 110}
#'   \item{x111}{Numeric covariate 111}
#'   \item{x112}{Numeric covariate 112}
#'   \item{x113}{Numeric covariate 113}
#'   \item{x114}{Numeric covariate 114}
#'   \item{x115}{Numeric covariate 115}
#'   \item{x116}{Numeric covariate 116}
#'   \item{x117}{Numeric covariate 117}
#'   \item{x118}{Numeric covariate 118}
#'   \item{x119}{Numeric covariate 119}
#'   \item{x120}{Numeric covariate 120}
#'   \item{x121}{Numeric covariate 121}
#'   \item{x122}{Numeric covariate 122}
#'   \item{x123}{Numeric covariate 123}
#'   \item{x124}{Numeric covariate 124}
#'   \item{x125}{Numeric covariate 125}
#'   \item{x126}{Numeric covariate 126}
#'   \item{x127}{Numeric covariate 127}
#'   \item{x128}{Numeric covariate 128}
#'   \item{x129}{Numeric covariate 129}
#'   \item{x130}{Numeric covariate 130}
#'   \item{x131}{Numeric covariate 131}
#'   \item{x132}{Numeric covariate 132}
#'   \item{x133}{Numeric covariate 133}
#'   \item{x134}{Numeric covariate 134}
#'   \item{x135}{Numeric covariate 135}
#'   \item{x136}{Numeric covariate 136}
#'   \item{x137}{Numeric covariate 137}
#'   \item{x138}{Numeric covariate 138}
#'   \item{x139}{Numeric covariate 139}
#'   \item{x140}{Numeric covariate 140}
#'   \item{x141}{Numeric covariate 141}
#'   \item{x142}{Numeric covariate 142}
#'   \item{x143}{Numeric covariate 143}
#'   \item{x144}{Numeric covariate 144}
#'   \item{x145}{Numeric covariate 145}
#'   \item{x146}{Numeric covariate 146}
#'   \item{x147}{Numeric covariate 147}
#'   \item{x148}{Numeric covariate 148}
#'   \item{x149}{Numeric covariate 149}
#'   \item{x150}{Numeric covariate 150}
#'   \item{x151}{Numeric covariate 151}
#'   \item{x152}{Numeric covariate 152}
#'   \item{x153}{Numeric covariate 153}
#'   \item{x154}{Numeric covariate 154}
#'   \item{x155}{Numeric covariate 155}
#'   \item{x156}{Numeric covariate 156}
#'   \item{x157}{Numeric covariate 157}
#'   \item{x158}{Numeric covariate 158}
#'   \item{x159}{Numeric covariate 159}
#'   \item{x160}{Numeric covariate 160}
#'   \item{x161}{Numeric covariate 161}
#'   \item{x162}{Numeric covariate 162}
#'   \item{x163}{Numeric covariate 163}
#'   \item{x164}{Numeric covariate 164}
#'   \item{x165}{Numeric covariate 165}
#'   \item{x166}{Numeric covariate 166}
#'   \item{x167}{Numeric covariate 167}
#'   \item{x168}{Numeric covariate 168}
#'   \item{x169}{Numeric covariate 169}
#'   \item{x170}{Numeric covariate 170}
#'   \item{x171}{Numeric covariate 171}
#'   \item{x172}{Numeric covariate 172}
#'   \item{x173}{Numeric covariate 173}
#'   \item{x174}{Numeric covariate 174}
#'   \item{x175}{Numeric covariate 175}
#'   \item{x176}{Numeric covariate 176}
#'   \item{x177}{Numeric covariate 177}
#'   \item{x178}{Numeric covariate 178}
#'   \item{x179}{Numeric covariate 179}
#'   \item{x180}{Numeric covariate 180}
#'   \item{x181}{Numeric covariate 181}
#'   \item{x182}{Numeric covariate 182}
#'   \item{x183}{Numeric covariate 183}
#'   \item{x184}{Numeric covariate 184}
#'   \item{x185}{Numeric covariate 185}
#'   \item{x186}{Numeric covariate 186}
#'   \item{x187}{Numeric covariate 187}
#'   \item{x188}{Numeric covariate 188}
#'   \item{x189}{Numeric covariate 189}
#'   \item{x190}{Numeric covariate 190}
#'   \item{x191}{Numeric covariate 191}
#'   \item{x192}{Numeric covariate 192}
#'   \item{x193}{Numeric covariate 193}
#'   \item{x194}{Numeric covariate 194}
#'   \item{x195}{Numeric covariate 195}
#'   \item{x196}{Numeric covariate 196}
#'   \item{x197}{Numeric covariate 197}
#'   \item{x198}{Numeric covariate 198}
#'   \item{x199}{Numeric covariate 199}
#'   \item{x200}{Numeric covariate 200}
#' }
NULL
