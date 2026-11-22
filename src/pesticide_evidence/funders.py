"""Funder-name classification used for the funder analysis.

N Northern public funder, S Southern public funder, I agrochemical industry, O unclassified or ambiguous.
Names are matched case-insensitively; a name matching both a Northern and a Southern pattern is ambiguous (O).
The 150 most frequent names carry the class assigned by hand in config/funder_classes_checked.csv.
"""

import csv
import re

from pesticide_evidence.config import ROOT

_NORTH = re.compile(
    r"National Institute of (Environmental|Allergy|General|Neurological|Child|Diabetes|Aging|Mental|Drug|Heart|Cancer|Arthritis|Occupational|Dental|Nursing|Deafness|Minority|Biomedical)|Gates|"
    r"National Institutes? of Health|NIH|NIEHS|Environmental Protection|\bEPA\b|Department of Agriculture|USDA|National Institute of Food and Agriculture|National Science Foundation|"
    r"Agricultural Research Service|Centers for Disease|Wellcome|Medical Research Council|Biotechnology and Biological|UK Research|Natural Environment Research|Engineering and Physical|"
    r"European|Horizon|Natural Sciences and Engineering|Fisheries and Oceans|Canadian Institutes|Japan Society|Japan Science|Agence Nationale|Deutsche Forsch|Fundação para a Ci|Fundacao para a Ci|"
    r"Ministerio de Econom|Ministerio de Ciencia|Ministerio de Educaci|Narodowe|National Research Foundation of Korea|Rural Development Administration|Australian Research|Grains Research|"
    r"Humboldt|Swiss National|Research Council of Norway|Danish|Swedish|Vetenskapsr|Formas|Academy of Finland|National Health and Medical|Israel Science|CSIRO|"
    r"Ministry of Education, Culture|Ministry of Agriculture, Forestry and Fisheries|Ministry of Health, Labour|Bundesministerium|Fonds|Netherlands Organisation|Research Foundation - Flanders|"
    r"Austrian Science|Agencia Estatal|Instituto de Salud Carlos|Junta de|Generalitat|Ministero|Consiglio Nazionale delle Ricerche|Department of Defense|U\.S\.|United States|Veterans|"
    r"National Aeronautics|Department of Energy|Medical Research Future|Marsden|Ministry of Business, Innovation|Korea Research|Korea Institute|Ministry of Science and ICT|"
    r"Ministry of Environment \(Korea|Leibniz|Max Planck|Helmholtz|Fraunhofer|CNRS|INRAE|ANSES|Region|Innovate UK|Department for Environment",
    re.I,
)
_SOUTH = re.compile(
    r"National Natural Science Foundation of China|Natural Science Foundation of|National Key R|People's Republic of China|Fundamental Research Funds|China Postdoctoral|Chinese Academy|China Agriculture|"
    r"Ministry of Science and Technology of the People|Ministry of Education of (the )?(People|China)|of China|Guangdong|Jiangsu|Zhejiang|Shandong|Conselho Nacional|CNPq|Coordenação de Aperfeiçoamento|CAPES|"
    r"Fundação de Amparo|FAPESP|FAPEMIG|FINEP|Consejo Nacional de Investigaciones|CONICET|Agencia Nacional de Promoci|Consejo Nacional de Ciencia|CONACYT|University Grants Commission|"
    r"Department of Science and Technology|Indian Council|Science and Engineering Research Board|Council of Scientific and Industrial|Department of Biotechnology|King Saud|King Abdul|Deanship|Egyptian|"
    r"Science, Technology|T[ÜU]B[İI]TAK|Higher Education Commission|Thailand|Thai|Vietnam|Serbia|Iran|Pakistan|Indonesia|Malaysia|Mexico|Chile|ANID|Colombia|Colciencias|Peru|Nigeria|South Africa|"
    r"National Research Foundation \(South|Kenya|Philippine|Russian|Ukrain|Turkey|Taiwan|Saudi|Brazil|India|Bangladesh|Sri Lanka|Hacettepe|Cairo|Tehran",
    re.I,
)
_INDUSTRY = re.compile(r"Bayer|Syngenta|BASF|Corteva|Dow AgroSciences|FMC|Sumitomo Chemical|Monsanto|DuPont|Nufarm|Adama|UPL|CropLife|Novartis|Ciba", re.I)


def _checked():
    """Hand-checked funder classes keyed by funder name."""
    with open(ROOT / "config" / "funder_classes_checked.csv", newline="") as f:
        return {row["name"]: row["class"] for row in csv.DictReader(f)}


CHECKED = _checked()


def classify_funder_by_rule(name):
    """Funder class from name patterns: I industry, S Southern, N Northern, O unclassified."""
    if _INDUSTRY.search(name):
        return "I"
    south, north = bool(_SOUTH.search(name)), bool(_NORTH.search(name))
    if south and not north:
        return "S"
    if north and not south:
        return "N"
    return "O"


def classify_funder(name):
    """Hand-checked class for the most frequent funder names, name-pattern class for the rest."""
    return CHECKED.get(name) or classify_funder_by_rule(name)
