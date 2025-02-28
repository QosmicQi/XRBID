###################################################################################
##########		Keeping track of all headers needed 		########### 
##########		across different Python scripts, etc. 		########### 
###################################################################################

###-----------------------------------------------------------------------------------------------------

### Data types from CSCView ###
# Master - best estimates based on combination of data from observations
# Stack - properties derived from a stacked observation
# Per-Observation - detections from a signle observation

mas = "master_source."		# no prefix in CSCView, but I use 'm'
stack = "stack_source."		# denoted with a prefix of 's' in CSCView
perob = "observation_source."	# denoted with a prefix of 'o' in CSCView

### Declaring headers retrieved from CSCView ###
sep = ".separation"
name = "master_source.name"
significance = "master_source.significance"
ra = "master_source.ra"
dec = "master_source.dec"
errel_r0 = "master_source.err_ellipse_r0"
errel_r1 = "master_source.err_ellipse_r1"
errel_ang = "master_source.err_ellipse_ang" 
signif = "master_source.significance" 
likeli = "master_source.likelihood_class" 
conf = "master_source.conf_flag" 
sat = "master_source.sat_src_flag"
streak = "master_source.streak_src_flag" 
flux_b = "master_source.flux_aper_b" 
flux_blo = "master_source.flux_aper_lolim_b" 
flux_bhi = "master_source.flux_aper_hilim_b" 
hm = "master_source.hard_hm" 
hm_lo = "master_source.hard_hm_lolim" 
hm_hi = "master_source.hard_hm_hilim" 
hs = "master_source.hard_hs" 
hs_lo = "master_source.hard_hs_lolim" 
hs_hi = "master_source.hard_hs_hilim" 
ms = "master_source.hard_ms" 
ms_lo = "master_source.hard_ms_lolim" 
ms_hi = "master_source.hard_ms_hilim" 
flux_bb = "master_source.flux_bb" 
flux_bblo = "master_source.flux_bb_lolim" 
flux_bbhi = "master_source.flux_bb_hilim" 
flux_pow = "master_source.flux_powlaw" 
flux_powlo = "master_source.flux_powlaw_lolim" 
flux_powhi = "master_source.flux_powlaw_hilim" 
var_flag = "master_source.var_flag"
ext_flag = "master_source.extent_flag"
t_counts = "stack_source.src_cnts_aper_b" 
t_counts_lo = "stack_source.src_cnts_aper_lolim_b" 
t_counts_hi = "stack_source.src_cnts_aper_hilim_b" 
h_counts = "stack_source.src_cnts_aper_h" 
h_counts_lo = "stack_source.src_cnts_aper_lolim_h" 
h_counts_hi = "stack_source.src_cnts_aper_hilim_h" 
m_counts = "stack_source.src_cnts_aper_m" 
m_counts_lo = "stack_source.src_cnts_aper_lolim_m" 
m_counts_hi = "stack_source.src_cnts_aper_hilim_m" 
s_counts = "stack_source.src_cnts_aper_s"
s_counts_lo = "stack_source.src_cnts_aper_lolim_s" 
s_counts_hi = "stack_source.src_cnts_aper_hilim_s"
u_counts = "stack_source.src_cnts_aper_u" 
u_counts_lo = "stack_source.src_cnts_aper_lolim_u" 
u_counts_hi = "stack_source.src_cnts_aper_hilim_u"
acis_time = "acis_time"

### Declaring possible headers used in source DataFrame ###
ID = "ID"
CSCID = "CSC ID"
LongID = "Long ID"
LongCSCID = "Long CSC ID"
LehmerID = "Lehmer ID"
Sig = "Significance"
RA = "RA"
Dec = "Dec"
X = "x"
Y = "y" 
PosErr = "Pos Err"
Bounds = "Bounds"
T_counts = "Total Counts"
T_countslo = "Total Counts Lo"
T_countshi = "Total Counts Hi"
U_counts = "Ultrasoft Counts"
U_countslo = "Ultrasoft Counts Lo"
U_countshi = "Ultrasoft Counts Hi"
S_counts = "Soft Counts"
S_countslo = "Soft Counts Lo"
S_countshi = "Soft Counts Hi"
M_counts = "Medium Counts"
M_countslo = "Medium Counts Lo"
M_countshi = "Medium Counts Hi"
H_counts = "Hard Counts"
H_countslo = "Hard Counts Lo"
H_countshi = "Hard Counts Hi"
HSRat = "HS Ratio"
HSlo = "HS Lo"
HShi = "HS Hi" 
HMRat = "HM Ratio"
HMerr = "HM Err"
HMlo = "HM Lo" 
HMhi = "HM Hi"
MSRat = "MS Ratio"
MSerr = "MS Err"
MSlo = "MS Lo" 
MShi = "MS Hi"
Flux_bb = "Flux (Blackbody)"
Flux_bblo = "Flux Lo (Blackbody)" 
Flux_bbhi = "Flux Hi (Blackbody)"
Flux_pow = "Flux (Power Law)"  # retrieved from CSCView
Flux_powlo = "Flux Lo (Power Law)"
Flux_powhi = "Flux Hi (Power Law)"
Spectra = "Spectra"            # spectral shape reported by Long (T, P, D) 
Variable = "Variable"          # reported by Long (Yes/No)
Var_flag = "Variability Flag"  # retrieved from CSCView (True/False)
Ext_flag = "Extended"          # Does CSC observe object as extended?
Class = "Class"
Cand = "Candidate Class"
LongClass = "Long Class"
Overlap = "Overlap"            # Overlap with Long sources
Crossrefs = "Crossrefs"        # Alternative IDs from other catalogs
Sep = "Dist from Center [arcsec]"
F1 = "0.35-8keV"
F2 = "0.35-1.1keV"
F3 = "1.1-2.6keV"
F4 = "2.6-8keV"
F5 = "0.2-2keV"
F6 = "2-8keV"
# More headers from CSCView
Fluxaper_b = "Aperture Corrected Flux (B)"
Fluxaper_blo = "Aperture Corrected Flux Lo (B)"
Fluxaper_bhi = "Aperture Corrected Flux Hi (B)"
Fluxaper_h = "Aperture Corrected Flux (H)"
Fluxaper_hlo = "Aperture Corrected Flux Lo (H)"
Fluxaper_hhi = "Aperture Corrected Flux Hi (H)"
Fluxaper_m = "Aperture Corrected Flux (M)"
Fluxaper_mlo = "Aperture Corrected Flux Lo (M)"
Fluxaper_mhi = "Aperture Corrected Flux Hi (M)"
Fluxaper_s = "Aperture Corrected Flux (S)"
Fluxaper_slo = "Aperture Corrected Flux Lo (S)"
Fluxaper_shi = "Aperture Corrected Flux Hi (S)"
Fluxaper_u = "Aperture Corrected Flux (U)"
Fluxaper_ulo = "Aperture Corrected Flux Lo (U)"
Fluxaper_uhi = "Aperture Corrected Flux Hi (U)"
Fluxaper_w = "Aperture Corrected Flux (W)"
Fluxaper_wlo = "Aperture Corrected Flux Lo (W)"
Fluxaper_whi = "Aperture Corrected Flux Hi (W)"
CSCTime = "Chandra Time (ks)"

# Headers from the evolutionary tracks (Padova?)
Zini = "Zini"
LogAge = "logAge"
Mini = "Mini"
IntIMF = "int_IMF"
Mass = "Mass"
LogL = "logL"
LogTe = "logTe"
Logg = "logg"
Label = "label"
McoreTP = "McoreTP"
C_O = "C_O"
Period0 = "period0"
Period1 = "period1"
Pmode = "pmode"
Mloss = "Mloss"
Tau = "tau1m"
Xc = "Xc"
Xn = "Xn"
Xo = "Xo"
Cexc = "Cexcess"
Z = "Z"
Mbol = "mbolmag"
# The following headers are from DAO Find
Mag = "Magnitude" 
Sharp = "Sharpness"
# The following headers are for by-eye classifications of sources
Class2 = "Secondary Class"
Conf = "Confidence Flag" 
Notes = "Notes"
# The following headers are for CMD files created with default settings
Umag = "Umag"
Bmag = "Bmag"
Vmag = "Vmag"
Rmag = "Rmag"
Imag = "Imag"
Jmag = "Jmag"
Hmag = "Hmag"
Kmag = "Kmag"
# The following headers are for CMD files created with WFC3 medium/wide filters
F390M = "F390M"
F410M = "F410M"
FQ422M = "FQ422M"
F467M = "F467M" 
F547M = "F547M" 
F621M = "F621M"
F689M = "F689M"
F763M = "F763M"
F845M = "F845M"
F098M = "F098M"
F127M = "F127M"
F149M = "F139M"
F153M = "F153M"
F218W = "F218W"
F225W = "F225W"
F275W = "F275W"
F336W = "F336W"
F290W = "F390W"
F438W = "F438W"
F475W = "F475W"
F555W = "F555W"
F606W = "F606W"
F625W = "F625W"
F775W = "F775W"
F814W = "F814W"
F105W = "F105W"
F110W = "F110W"
F125W = "F125W"
F140W = "F140W"
F160W = "F160W"
# The following are headers from the evolutionary tracks (https://people.sissa.it/~sbressan/parsec.html)
Model = "Model"
Age = "Age" 
LogL = "Luminosity (log)" 
LogTe = "Effective Temp (log)"
LogR = "Radius (log)" 
LogRat = "Mass loss rate (log)"
McoreHe = "Mass of H-exhausted core (M_sun)"
McoreC = "Mass of  He-exhauseted core (M_sun)"
HCen = "Central H composition"
HeCen = "Central He composition" 
CCen = "Central C composition"
OCen = "Central O composition" 
LX = "Fractional Lx" 
LY = "Fractional Ly" 
LC = "Fractional Lc" 
LNeut = "Fractional Lneutr"
HSur = "Surface H composition" 
HeSur = "Surface He composition"
CSur = "Surface C composition" 
NSur = "Surface N composition" 
OSur = "Surface S composition"
Phase = "Phase"

# Headers used in Hardness Ratio Files
TH = "Total-Hard"
TM = "Total-Medium"
TS = "Total-Soft"
TU = "Total-Ultrasoft"
HM = "Hard-Medium"
HS = "Hard-Soft"
HU = "Hard-Ultrasoft"
MS = "Medium-Soft"
MU = "Medium-Ultrasoft"
L = "Luminosity"
AbsHS = "Abs. " + HS
AbsHM = "Abs. " + HM
AbsMS = "Abs. " + MS

# Headers for mass tracks
M1 = "1 Msun" 
M2 = "2 Msun"
M3 = "3 Msun" 
M5 = "5 Msun" 
M8 = "8 Msun" 
M20 = "20 Msun" 
M40 = "40 Msun"

M1x = "1 Msun (x)" 
M2x = "2 Msun (x)"
M3x = "3 Msun (x)" 
M5x = "5 Msun (x)" 
M8x = "8 Msun (x)" 
M20x = "20 Msun (x)" 
M40x = "40 Msun (x)"

M1y = "1 Msun (y)" 
M2y = "2 Msun (y)"
M3y = "3 Msun (y)" 
M5y = "5 Msun (y)" 
M8y = "8 Msun (y)" 
M20y = "20 Msun (y)" 
M40y = "40 Msun (y)"

# Headers for apertures/photometry
Aphalf = "Aper. 0.5 mag"
Ap1 = "Aper. 1 mag"
Ap2 = "Aper. 2 mag" 
Ap3 = "Aper. 3 mag" 
Ap4 = "Aper. 4 mag" 
Ap5 = "Aper. 5 mag" 
Ap6 = "Aper. 6 mag" 
Ap7 = "Aper. 7 mag" 
Ap8 = "Aper. 8 mag" 
Ap9 = "Aper. 9 mag"
Ap10 = "Aper. 10 mag"
Ap15 = "Aper. 15 mag"
Ap20 = "Aper. 20 mag"
Filter = "Filter"

# CMD analysis headers
V = "V" 
B = "B" 
I = "I"
U = "U"
VI = "V-I"
BV = "B-V"
BI = "B-I"
UB = "U-B"
Verr = "V Err" 
Berr = "B Err" 
Ierr = "I Err"
Uerr = "U Err"

# Headers from Lehmer et. al (2019)
LogF = "Flux (log)"
Gal = "Galaxy" 
Offset = "Offset" 
NetCount = "Net Counts" # from Lehmer, 0.5-8 keV
NetCount_err = "Net Count Uncertainty" # 1sig
Gamma = "Gamma" # Best-fit photon index
Gamma_err = "Gamma Uncertainty" 
Loc = "Location Flag" 


# Misc. headers
Blair = "Blair"
HLA = "HLA"
Cat = "Catalog"
Catalog = Cat
Radius = "Radius"
DaoNo = "Dao No"
Greenx = "Green x"
Greeny = "Green y" 
Redx = "Red x"
Redy = "Red y" 
Bluex = "Blue x"
Bluey = "Blue y"

Field = "Field"
F = "Field"
Sig1 = "1sig Radius" 
Sig2 = "2sig Radius" 
Location = "Location" 

# list of all known headers; add more as needed
heads = [ID, CSCID, CSCID, LongID, LongCSCID, LehmerID, Sig, RA, Dec, X, Y, PosErr, Bounds, T_counts, T_countslo, T_countshi, \
	U_counts, U_countslo, U_countshi, S_counts, S_countslo, S_countshi, M_counts, M_countslo, M_countshi, \
	H_counts, H_countslo, H_countshi, HS, HSlo, HShi, HM, HMerr, HMlo, HMhi, MS, MSerr, MSlo, MShi, \
	Flux_bb, Flux_bblo, Flux_bbhi, Flux_pow, Flux_powlo, Flux_powhi, Spectra, \
	Variable, Var_flag, Ext_flag, Class, LongClass, Overlap, Crossrefs, Sep, F1, F2, F3, F4, F5, F6, Zini, LogAge, \
	Mini, IntIMF, Mass, LogL, LogTe, Logg, Label, McoreTP, C_O, Period0, Period1, Pmode, Mloss, Tau, Xc, \
	Xn, Xo, Cexc, Z, Mbol, Mag, Sharp, Class2, Conf, Notes, Umag, Bmag, Vmag, Rmag, Imag, Jmag, Hmag, Kmag, \
	F390M, F410M, FQ422M, F467M, F547M, F621M, F689M, F763M, F845M, F098M, F127M, F149M, F153M, F218W, F225W, \
	F275W, F336W, F290W, F438W, F475W, F555W, F606W, F625W, F775W, F814W, F105W, F110W, F125W, F140W, F160W, \
	Model, Age, LogL, LogTe, LogR, LogRat, McoreHe, McoreC, HCen, HeCen, CCen, OCen, LX, LY, LC, LNeut, HSur, HeSur, CSur, NSur, OSur, Phase, \
	TH, TM, TS, TU, HM, HS, HU, MS, MU, L, Fluxaper_b, Fluxaper_blo, Fluxaper_bhi, Fluxaper_h, \
	Fluxaper_hlo, Fluxaper_hhi, Fluxaper_m, Fluxaper_mlo, Fluxaper_mhi, Fluxaper_s, Fluxaper_slo, \
	Fluxaper_shi, Fluxaper_u, Fluxaper_ulo, Fluxaper_uhi, Fluxaper_w, Fluxaper_wlo, Fluxaper_whi, AbsHM, AbsHS, AbsMS, Aphalf, \
	Ap1, Ap2, Ap3, Ap4, Ap5, Ap6, Ap7, Ap8, Ap9, Ap10, Ap15, Ap20, V, B, I, U, VI, BV, BI, UB, Verr, Berr, Ierr, Uerr, Filter, M1, M2, M3, M5, M8, M20, M40, Cand, LogF, Gal, Cat, Offset, NetCount, NetCount_err, Gamma, Gamma_err, Loc, Blair, HLA, Catalog, Radius, DaoNo, Greenx, Greeny, Redx, Redy, Bluex, Bluey, CSCTime, Sig1, Sig2, Location]

lowerheads = [x.lower() for x in heads]

# tabheads are for headers from CSCView. The indices should match those of the corresponding headers in 'heads' above.
tabheads = [None, name, "name", None, None, None, significance, ra, dec, None, None, None, "Out", t_counts, t_counts_lo, t_counts_hi, \
	u_counts, u_counts_lo, u_counts_hi, s_counts, s_counts_lo, s_counts_hi, m_counts, m_counts_lo, m_counts_hi, \
	h_counts, h_counts_lo, h_counts_hi, hs, hs_lo, hs_hi, hm, None, hm_lo, hm_hi, ms, None, ms_lo, ms_hi, flux_bb, \
	flux_bblo, flux_bbhi, flux_pow, flux_powlo, flux_powhi, "N/A", "No", var_flag, ext_flag, None, None, "No", None, sep, \
	None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, \
	None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, \
	None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, \
	None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, \
	None, None, None, None, None, None, None, None, None, None, None, None, None, \
	None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, \
	None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, \
	None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, \
	None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, \
	None, None, None, None, None, None, None, acis_time, None, None, None]

# For a simpler way of dealing with these, define a dictionary
# Depending on how the table is formatted, may need to remove the split. 
# For now, the headers read in from CSC do not have the master_source prefix.
tabheads_split = [tabheads[i].split(".")[-1] for i in range(2, len(tabheads)) if tabheads[i] != None]

headers_dict = dict([(tabheads[i].split(".")[-1], heads[i]) for i in range(len(tabheads)) if tabheads[i] != None])


