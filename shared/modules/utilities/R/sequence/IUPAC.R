
# declare tables of lookup values for all base combinations that have either:
#       an output IUPAC code other than the default of N
#       and output SW score other than the default of mismatchPenalty

# three-base degenerate codes are not used, only bases and 2-base degeneracy
# thus, a base set with three possible bases defaults to output N
# and code BDHV will never match anything

# e.g.  ACGTTCAARCC ----> ACRTTYAANCC
#       ACATTYAAYCC

# paired base values for incoming raw bases (A C G T)
baseMatches <-  list( # score with full matchScore in SW
    AA = 'A',
    CC = 'C',
    GG = 'G',
    TT = 'T'    
)
baseMismatches <- list( # score with full mismatchPenalty in SW
    AG = 'R',
    GA = 'R',
    
    CT = 'Y',
    TC = 'Y',
    
    GC = 'S',
    CG = 'S',
    
    AT = 'W',
    TA = 'W',
    
    GT = 'K',
    TG = 'K',
    
    AC = 'M',
    CA = 'M'
)

# paired base values that include RYSWKM, i.e. partial degeneracy
ryswkmMatches <- list( # score with half matchScore in SW
    RR = 'R',
    RA = 'R',
    RG = 'R',
    AR = 'R',
    GR = 'R',
    
    YY = 'Y',
    YC = 'Y',
    YT = 'Y',
    CY = 'Y',
    TY = 'Y',
    
    SS = 'S',
    SG = 'S',
    SC = 'S',
    GS = 'S',
    CS = 'S',
    
    WW = 'W',
    WA = 'W',
    WT = 'W',
    AW = 'W',
    TW = 'W',
    
    KK = 'K',
    KG = 'K',
    KT = 'K',
    GK = 'K',
    TK = 'K',
    
    MM = 'M',
    MA = 'M',
    MC = 'M',
    AM = 'M',
    CM = 'M'
)

# paired base values for any combination that includes N, i.e. full degeneracy
nMatches <- list( # score with neutral=0 in SW (neither penalize nor promote)
    NN = 'N',    # thus, a base previously declared uninformative is an M operation but has no alignment value
    
    'NA' = 'N',
    NC = 'N',
    NG = 'N',
    NT = 'N',
    AN = 'N',
    CN = 'N',
    GN = 'N',
    TN = 'N',
    
    NR = 'N',
    NY = 'N',
    NS = 'N',
    NW = 'N',
    NK = 'N',
    NM = 'N',
    RN = 'N',
    YN = 'N',  
    SN = 'N',
    WN = 'N',
    KN = 'N',
    MN = 'N'
)

# all combinations that do not default to output base N
# used to generate consensus sequences
mOperations <- c(baseMatches, baseMismatches, ryswkmMatches)

# combinations for doing ungapped comparisons via sub noIndelMatch
ACGTNMatches   <- c(baseMatches, nMatches)
noIndelMatches <- c(baseMatches, baseMismatches, nMatches)

# combinations for assembling consesnuses with missing data
spaceMatches <- list(
    ' A' = 'A',
    ' C' = 'C',
    ' G' = 'G',
    ' T' = 'T',
    'A ' = 'A',
    'C ' = 'C',
    'G ' = 'G',
    'T ' = 'T',
    '  ' = ' '
)

#A	Adenine
#C	Cytosine
#G	Guanine
#T (or U)	Thymine (or Uracil)
#R	A or G
#Y	C or T
#S	G or C
#W	A or T
#K	G or T
#M	A or C
#B	C or G or T
#D	A or G or T
#H	A or C or T
#V	A or C or G
#N	any base
#. or -	gap
