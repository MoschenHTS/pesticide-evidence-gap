"""PubMed query blocks and the text rules used to classify records: model groups, outcome domains, methods, design features, exclusions."""

import re

from .names import AMBIG, CLASS_WORDS_SAFE, NAMES
from .pubmed import OR, tiab

Y0, Y1 = 2005, 2025

MODEL_QUERY = {
    "Zebrafish": tiab("zebrafish", "danio rerio", "zebra fish", "Brachydanio", "D. rerio"),
    "Other fish": tiab(
        "fish",
        "fishes",
        "carp",
        "trout",
        "medaka",
        "Oryzias",
        "tilapia",
        "croaker",
        "salmon",
        "fathead minnow",
        "catfish",
        "goldfish",
        "killifish",
        "stickleback",
        "sturgeon",
        "eel",
        "bass",
        "guppy",
        "Poecilia",
        "loach",
        "minnow",
        "flounder",
        "seabream",
        "sea bream",
        "grouper",
        "pangasius",
    ),
    "Amphibians": tiab("Xenopus", "tadpole*", "amphibian*", "frog", "frogs", "toad", "toads", "salamander*", "newt", "newts", "Rana", "Lithobates"),
    "Rodents": tiab("mouse", "mice", "rat", "rats", "rodent*", "murine", "Wistar", "Sprague-Dawley", "C57BL/6", "hamster*", "guinea pig*"),
    "Other mammals": tiab(
        "rabbit*",
        "swine",
        "pig",
        "pigs",
        "porcine",
        "piglet*",
        "sheep",
        "ovine",
        "cattle",
        "bovine",
        "calf",
        "calves",
        "dog",
        "dogs",
        "canine",
        "monkey*",
        "macaque*",
        "primate*",
        "goat",
        "goats",
        "horse",
        "horses",
        "bat",
        "bats",
    ),
    "Birds": tiab("bird", "birds", "avian", "chicken*", "quail", "duck", "ducks", "poultry", "chick embryo*", "broiler*", "raptor*", "passerine*", "pigeon*"),
    "Aquatic invertebrates": tiab(
        "Daphnia",
        "Gammarus",
        "shrimp",
        "prawn*",
        "crustacean*",
        "mussel*",
        "oyster*",
        "Chironomus",
        "copepod*",
        "Hyalella",
        "amphipod*",
        "clam",
        "clams",
        "bivalve*",
        "snail",
        "snails",
        "gastropod*",
        "crab",
        "crabs",
        "sea urchin*",
        "Ceriodaphnia",
        "mysid*",
        "rotifer*",
        "Artemia",
        "planarian*",
        "Hydra",
        "mosquito larvae",
    ),
    "Bees and pollinators": tiab(
        "honey bee*",
        "honeybee*",
        "Apis mellifera",
        "bumble bee*",
        "bumblebee*",
        "Bombus",
        "stingless bee*",
        "Melipona",
        "bee",
        "bees",
        "pollinator*",
        "Osmia",
        "solitary bee*",
    ),
    "Soil and other invertebrates": tiab(
        "earthworm*",
        "Eisenia",
        "Caenorhabditis",
        "C. elegans",
        "nematode*",
        "Drosophila",
        "springtail*",
        "Folsomia",
        "collembola*",
        "isopod*",
        "soil invertebrate*",
        "Enchytraeid*",
        "predatory mite*",
        "parasitoid*",
        "ladybird*",
        "lacewing*",
        "ground beetle*",
        "carabid*",
        "natural enem*",
    ),
    "Cells": "("
    + OR(
        tiab(
            "cell line*",
            "cell culture*",
            "cultured cells",
            "primary culture*",
            "in vitro",
            "HepG2",
            "SH-SY5Y",
            "hepatocyte*",
            "fibroblast*",
            "Caco-2",
            "HEK293",
            "HeLa",
            "neuroblastoma",
            "PC12",
            "Neuro-2a",
            "BV2",
            "organoid*",
            "iPSC",
            "stem cell*",
            "lymphocyte*",
            "macrophage*",
            "keratinocyte*",
            "cardiomyocyte*",
            "myoblast*",
            "C2C12",
            "3T3-L1",
            "A549",
            "HaCaT",
            "Sertoli",
            "Leydig",
            "granulosa",
            "spermatozoa",
            "sperm",
        ),
        '"Cell Line"[mh]',
        '"Cells, Cultured"[mh]',
    )
    + ")",
}
ANY_MODEL = "(" + OR(*MODEL_QUERY.values(), '"Animals"[mh]', '"Disease Models, Animal"[mh]') + ")"

TOX = tiab(
    "toxic*",
    "ecotox*",
    "non-target",
    "nontarget",
    "adverse effect*",
    "adverse outcome",
    "neurotox*",
    "hepatotox*",
    "cytotox*",
    "genotox*",
    "oxidative stress",
    "bioaccumulat*",
    "bioconcentrat*",
    "developmental toxicity",
    "mitochondrial dysfunction",
    "exposed to",
    "exposure to",
    "effects of",
    "impact",
    "lethal",
    "mortality",
    "LC50",
    "EC50",
    "sublethal",
    "sub-lethal",
    "endocrine disrupt*",
    "teratogen*",
    "behavio*",
    "apoptosis",
    "DNA damage",
    "histopatholog*",
    "immunotox*",
    "reproductive",
)
FR = ' NOT ("flame retardant*"[tiab] OR TPhP[tiab] OR TDCPP[tiab] OR "organophosphate ester*"[tiab] OR "organophosphorus flame"[tiab])'
GEN_Q = (
    "("
    + OR(
        tiab(
            "pesticide*",
            "herbicide*",
            "insecticide*",
            "fungicide*",
            "agrochemical*",
            "agrotoxic*",
            "acaricide*",
            "nematicide*",
            "rodenticide*",
            "molluscicide*",
            "phytosanitary",
            "crop protection",
            "plant protection product*",
            "agricultural chemical*",
        ),
        '"Pesticides"[mh]',
        '"Herbicides"[mh]',
        '"Insecticides"[mh]',
        '"Fungicides, Industrial"[mh]',
        '"Acaricides"[mh]',
        '"Rodenticides"[mh]',
        '"Molluscacides"[mh]',
        '"Nematocides"[mh]',
        '"Insect Repellents"[mh]',
        '"Pesticide Residues"[mh]',
        '"Pesticides"[pa]',
        '"Herbicides"[pa]',
        '"Insecticides"[pa]',
        '"Fungicides, Industrial"[pa]',
        '"Acaricides"[pa]',
        '"Rodenticides"[pa]',
        '"Molluscacides"[pa]',
    )
    + ")"
)


def names_for_class(g, c, include_ambiguous=False):
    """Ingredient names and class words for a group and class, leaving out ambiguous names unless requested."""
    L = [n for n in NAMES[(g, c)] if (include_ambiguous or n not in AMBIG) and len(n) >= 3]
    L += CLASS_WORDS_SAFE.get(c, [])
    return L


def batches(L, k):
    """List split into consecutive batches of k items."""
    return [L[i : i + k] for i in range(0, len(L), k)]


OUTCOME_RE = {
    "Acute toxicity, survival and development": r"mortality|survival|lethal|\blc50\b|\bec50\b|\bld50\b|malformation|teratogen|hatching|growth (inhibition|retardation)|body (weight|length)|developmental (toxicity|abnormal|delay)|morphological",
    "Mitochondria and bioenergetics": r"mitochondri|bioenerget|oxidative phosphorylation|respirometr|oxygen consumption|electron transport|complex (i|ii|iii|iv)\b|membrane potential|atp (content|production|level|synthesis)",
    "Metabolism (lipid, glucose, obesity)": r"obes|obesogen|adipos|adipogen|lipid metabolism|lipid accumulation|lipid homeostasis|steatosis|fatty liver|triglycerid|cholesterol|glucose metabolism|glycolipid|hyperglyc|diabet|lipidom|insulin|energy metabolism",
    "Behavior and cognition": r"cogniti|learning|memory|neurobehavio|anxiety|behavio|locomot|swimming|shoaling|novel tank|startle|open field|maze|forced swim|depress",
    "Neurotoxicity": r"neurodevelop|neurotox|neurogenesis|axon|motor neuron|neural development|acetylcholinesterase|dopaminergic|neurodegenerat|\bache\b|cholinesterase|neuronal",
    "Endocrine disruption and reproduction": r"endocrine|thyroid|estrogen|androgen|estradiol|testosterone|vitellogenin|hormone|reproduct|spermatogen|oogenesis|fertil|gonad|sperm quality|ovar|testi",
    "Offspring and multigenerational effects": r"transgeneration|intergeneration|multigeneration|multi-generation|parental exposure|maternal exposure|paternal|f1 generation|f2 generation|offspring|progeny|prenatal exposure|gestation|lactation|perinatal",
    "Gut microbiota": r"microbiota|microbiome|gut microb|intestinal flora|16s|dysbiosis",
    "Immunity and inflammation": r"inflammat|immunotox|immune|cytokine|macrophage|neutrophil|interleukin|\btnf",
    "Oxidative stress": r"oxidative stress|reactive oxygen|\bros\b|malondialdehyde|\bmda\b|superoxide dismutase|catalase|glutathione|antioxidant|lipid peroxidation|\bnrf2",
    "Genotoxicity and cell death": r"genotox|comet assay|micronucle|dna damage|apoptos|tunel|caspase|ferroptosis|pyroptosis|cell death|cytotox|chromosom|autophag",
    "Organ toxicity and histopathology": r"hepatotox|nephrotox|cardiotox|liver (damage|injury)|kidney (damage|injury)|histopatholog|histolog|heart rate|pericardial|cardiac|renal|hepatic|pathological (change|alteration)|lung injury|pulmonary|gill",
}
METHOD_RE = {
    "Behavioral tracking": r"behavio|locomot|swimming|anxiety|novel tank|light-dark|shoaling|open field|photomotor|startle|elevated plus|forced swim",
    "Learning and memory tests": r"learning|memory|cogniti|avoidance|y-maze|t-maze|novel object|morris|water maze",
    "Mitochondrial function assays": r"oxygen consumption|respirometr|seahorse|mitochondrial respiration|membrane potential|mitochondrial function|oxidative phosphorylation|atp (content|production|level)|complex (i|ii|iii|iv)\b|electron transport",
    "Oxidative stress markers": r"oxidative stress|reactive oxygen|\bros\b|malondialdehyde|\bmda\b|superoxide dismutase|catalase|glutathione",
    "RT-qPCR": r"qpcr|rt-pcr|real-time (pcr|polymerase)|quantitative (real-time )?pcr|qrt-pcr|rt-qpcr",
    "RNA-seq or transcriptomics": r"rna-seq|rna sequencing|transcriptom|microarray",
    "Metabolomics or lipidomics": r"metabolom|lipidom",
    "Proteomics": r"proteom",
    "Integrative multi-omics": r"multi-omics|multiomics|integrative omics|integrated omics|integrative analysis",
    "Gut microbiota (16S)": r"16s|microbiota|microbiome|gut microb",
    "Epigenetic assays": r"dna methylation|epigenom|epigenetic|bisulfite|histone",
    "Histopathology": r"histolog|histopatholog|hematoxylin|h&e",
    "Chemical analysis (LC-MS, GC-MS)": r"lc-ms|gc-ms|uplc|hplc|mass spectrometry|bioaccumulat|bioconcentrat|tissue residue|toxicokinetic",
    "Molecular docking or dynamics": r"docking|molecular dynamics",
    "Mutant, transgenic or CRISPR lines": r"crispr|knockout|knock-out|mutant|transgenic|tg\(|morpholino",
    "Imaging (confocal, in situ, fluorescence)": r"confocal|fluorescen|whole-mount|in situ hybridi|immunofluorescen|light-sheet|immunohistochem",
    "Cardiac assays": r"heart rate|cardiac|cardiotox|electrocardio",
    "Hormone and endocrine assays": r"thyroid|estradiol|testosterone|vitellogenin|hormone|endocrine|elisa",
    "Apoptosis and cell death": r"apoptos|tunel|caspase|ferroptosis|pyroptosis|annexin",
    "Genotoxicity assays": r"genotox|comet assay|micronucle|dna damage",
    "Machine learning or in silico prediction": r"machine learning|qsar|deep learning|in silico|adverse outcome pathway|network toxicology",
    "Mixture designs": r"mixture|co-exposure|coexposure|combined exposure|joint toxicity|binary",
}
DESIGN_RE = {
    "Early-life exposure (embryo, larva, neonate)": r"embryo|larva|\bhpf\b|\bdpf\b|hatching|neonat|prenatal|in utero|gestation|pup\b|pups\b|tadpole|nestling",
    "Adult or juvenile exposure": r"\badult|juvenile|young (rats|mice)|weaned",
    "Chronic exposure (28 days or longer)": r"chronic|long-term|28[- ]day|30[- ]day|60[- ]day|90[- ]day|\b(4|8|12|16) weeks|months|subchronic|sub-chronic",
    "Environmentally relevant concentrations": r"environmentally relevant|environmental concentration|ecologically relevant|low-dose|low dose|sublethal|sub-lethal|µg/l|μg/l|ng/l|field-realistic|field realistic",
    "Both sexes analysed": r"male and female|both sexes|sex-specific|sex differen|sexually dimorph|sex-dependent|sex-based",
    "Multigenerational design": r"transgeneration|intergeneration|multigeneration|parental exposure|maternal exposure|paternal exposure|f1 generation|f2 generation|offspring",
    "Commercial formulation or product": r"commercial formulation|formulated product|commercial product|technical grade|roundup|formulation|commercial-grade",
    "Water or tissue concentration measured": r"measured concentration|actual concentration|tissue concentration|bioaccumulat|bioconcentrat|residue",
}
EXTRA_RE = {"hfd_diet": r"high-fat|high fat|obesogenic diet|western diet|diet-induced|overfeed"}

MODEL_RE = {
    "Zebrafish": r"zebrafish|danio rerio|zebra fish|brachydanio|\bd\. rerio",
    "Other fish": r"\bfish(es)?\b|\bcarp\b|\btrout\b|medaka|oryzias|tilapia|croaker|salmon|fathead minnow|catfish|goldfish|killifish|stickleback|sturgeon|\beels?\b|\bbass\b|guppy|poecilia|loach|minnow|flounder|sea ?bream|grouper|pangasius|carassius|cyprinus|oncorhynchus|\bgambusia|\bpimephales",
    "Amphibians": r"xenopus|tadpole|amphibian|\bfrogs?\b|\btoads?\b|salamander|\bnewts?\b|\brana\b|lithobates|bufo|hyla\b",
    "Rodents": r"\bmouse\b|\bmice\b|\brats?\b|rodent|murine|wistar|sprague|c57bl|hamster|guinea pig|\bmus musculus|rattus",
    "Other mammals": r"rabbit|\bswine\b|\bpigs?\b|porcine|piglet|\bsheep\b|ovine|cattle|bovine|\bcalf\b|\bcalves\b|\bdogs?\b|canine|monkey|macaque|primate|\bgoats?\b|\bhorses?\b|\bbats?\b|\bcows?\b|dairy",
    "Birds": r"\bbirds?\b|avian|chicken|\bquail|\bducks?\b|poultry|chick embryo|broiler|raptor|passerine|pigeon|\bgallus|songbird|\bgeese\b|goose|\bhens?\b",
    "Aquatic invertebrates": r"daphni|gammarus|shrimp|prawn|crustacean|mussel|oyster|chironomus|copepod|hyalella|amphipod|\bclams?\b|bivalve|\bsnails?\b|gastropod|\bcrabs?\b|sea urchin|ceriodaphnia|mysid|rotifer|artemia|planarian|\bhydra\b|crayfish|lobster|mollus|\bscallop|cladocera|aquatic invertebrate|macroinvertebrate|ephemeropter|trichopter|\bmayfl",
    "Bees and pollinators": r"honey ?bee|apis mellifera|bumble ?bee|\bbombus|stingless bee|melipona|\bbees?\b|pollinator|\bosmia|solitary bee|\bmegachile|butterfl|monarch",
    "Soil and other invertebrates": r"earthworm|eisenia|caenorhabditis|c\. elegans|nematode|drosophila|springtail|folsomia|collembol|isopod|soil invertebrate|enchytraeid|predatory mite|parasitoid|ladybird|ladybeetle|lacewing|ground beetle|carabid|natural enem|spider|soil fauna|mesofauna|\bants?\b|silkworm|bombyx",
    "Cells": r"cell lines?\b|cell culture|cultured cell|primary culture|in vitro|hepg2|sh-sy5y|hepatocyte|fibroblast|caco-2|hek293|\bhela\b|neuroblastoma|\bpc12|neuro-2a|\bbv2\b|organoid|ipsc|stem cell|lymphocyte|keratinocyte|cardiomyocyte|myoblast|c2c12|3t3-l1|\ba549\b|hacat|sertoli|leydig|granulosa|\bcells\b",
}
HUMAN_RE = r"epidemiolog|cohort|case-control|poisoning|poisoned|patients?\b|self-harm|suicid|occupational exposure|farmworker|pesticide applicators?|biomonitoring|urinary metabolite|children|pregnan[ct]\w* wom[ae]n|cross-sectional|odds ratio|hazard ratio"
TARGET_PEST_RE = r"aedes|anopheles|culex|mosquito|aphid|spodoptera|plutella|helicoverpa|tribolium|whitefl|bemisia|planthopper|nilaparvata|tetranychus|spider mite|thrip|leafhopper|cotton bollworm|armyworm|diamondback|sitophilus|callosobruchus|musca domestica|house ?fl|cockroach|blattella|\bticks?\b|rhipicephalus|\blice\b|louse|bed ?bug|termite|\bweeds?\b|\bfusarium|botrytis|phytophthora|magnaporthe|\brust\b|powdery mildew|plasmopara"

DELIV = re.compile(
    r"fluorescen\w*[ -](probe|detection|sensor)|colorimetric|biosensor|electrochemical sensor|pesticide delivery|nanocarrier|controlled[- ]release|metal-organic framework|\bMOF\b|molecularly imprinted|stimuli-responsive|smart delivery",
    re.I,
)
SYNTH = re.compile(
    r"design, synthesis|synthesi[sz]ed|synthesis and (biological|antifungal|herbicidal|insecticidal|bioactivity)|structure-activity relationship|novel (derivatives|compounds) (as|with)|bioactivity of (novel|new)",
    re.I,
)
EFFIC = re.compile(
    r"insecticide resistance|herbicide resistance|fungicide resistance|resistan(ce|t) (mechanism|to|management|strain|population|evolution)|bio-?efficacy|efficacy (of|against)|field efficacy|larvicidal|ovicidal|adulticidal|knockdown resistance|\bkdr\b|susceptibility (to|of)|control (of|efficacy)|mode of action in (insects|pests)|toxicity to (aphid|mosquito|larvae of|pests?)|lethal concentration for (pest|mosquito)",
    re.I,
)
RESID = re.compile(
    r"residue|multiresidue|quechers|analytical method|determination of|detection of|extraction method|monitoring of|dissipation|photodegradation|biodegradation|degradation (of|by)|sorption|adsorption|removal of|remediation|occurrence of|contamination (of|in)|risk assessment of dietary|dietary exposure",
    re.I,
)
STRONG_TOX = re.compile(
    r"toxic|ecotox|lethal|\\blc50|\\bec50|sublethal|sub-lethal|oxidative stress|genotox|cytotox|neurotox|hepatotox|apoptos|biomarker|malformation|teratogen|behavio|mortality|histopath|bioassay|embryo|larva|developmental",
    re.I,
)
BENEFICIAL = re.compile(
    r"selectivity|non-?target|beneficial (insect|arthropod|organism)|natural enem|predatory (mite|insect|arthropod|bug|beetle)|\bpredators?\b|parasitoid|pollinator|biological control agent|side[- ]effects? (on|of)",
    re.I,
)
PD_TOOL_T = re.compile(r"parkinson", re.I)
PD_TOOL_C = re.compile(r"rotenone|paraquat|mptp", re.I)
PEST_CTX = re.compile(
    r"pesticid|insecticid|herbicid|fungicid|agrochemical|agrotoxic|acaricid|nematicid|rodenticid|molluscicid|oomycetic|miticid|larvicid|phytosanitary|crop protection|plant protection product|agricultural chemical",
    re.I,
)
