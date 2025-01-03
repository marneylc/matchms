import os
from matchms.importing import load_from_mgf

# parameters to tune
# intensity_from=0.1
# tolerance=0.005

#############################
## import a library
#############################
path_data = "."
file_mgf = os.path.join(path_data, "GNPS-NIH-NATURALPRODUCTSLIBRARY.mgf")
spectrums = list(load_from_mgf(file_mgf))
import matchms.filtering as ms_filters
def metadata_processing(spectrum):
    spectrum = ms_filters.default_filters(spectrum)
    spectrum = ms_filters.repair_inchi_inchikey_smiles(spectrum)
    spectrum = ms_filters.derive_inchi_from_smiles(spectrum)
    spectrum = ms_filters.derive_smiles_from_inchi(spectrum)
    spectrum = ms_filters.derive_inchikey_from_inchi(spectrum)
    spectrum = ms_filters.harmonize_undefined_smiles(spectrum)
    spectrum = ms_filters.harmonize_undefined_inchi(spectrum)
    spectrum = ms_filters.harmonize_undefined_inchikey(spectrum)
    spectrum = ms_filters.add_precursor_mz(spectrum)
    return spectrum

def peak_processing(spectrum):
    spectrum = ms_filters.default_filters(spectrum)
    spectrum = ms_filters.normalize_intensities(spectrum)
    spectrum = ms_filters.select_by_intensity(spectrum, intensity_from=0.01)
    spectrum = ms_filters.select_by_mz(spectrum, mz_from=40, mz_to=1300)
    return spectrum

spectrums = [metadata_processing(s) for s in spectrums] 
spectrums = [peak_processing(s) for s in spectrums]

########################################
## import experimental spectra
########################################
path_data = "."
file_mgf_experimental = os.path.join(path_data, "Emily_WS_Tea_PeakTableFragments.mgf")
experimental_spectrums = list(load_from_mgf(file_mgf))
experimental_spectrums = [peak_processing(s) for s in experimental_spectrums]

# plot individual spectra
from matplotlib import pyplot as plt
i = 35
experimental_spectrums[i].plot()
experimental_spectrums[i].metadata
plt.savefig("experimental_spectrum", dpi=300)

################################################################################
# need a check for peaks that have no MSMS spectra (an artifact from progenesis processing), Num_Peaks: 0, is a way to catch them
# can do with a vim macro, save it as "@q" and run it for the whole file with ":%:normal @q"
################################################################################

# room for implementing any of the cleaning tools here

###########################
## run cosine
###########################
from matchms import calculate_scores
from matchms.similarity import CosineGreedy
from matchms.similarity import ModifiedCosine 

# select the algorithm
similarity_measure = ModifiedCosine(tolerance=0.005)
#similarity_measure = CosineGreedy(tolerance=0.05)

# calculate scores
scores = calculate_scores(experimental_spectrums, spectrums, similarity_measure, is_symmetric=False)

# inspect results by score
scores_array = scores.scores.to_array()
scores_array[:5, :5]["ModifiedCosine_score"]
scores_array[:5, :5]["ModifiedCosine_matches"]
#scores_array[:5, :5]["CosineGreedy_score"]
#scores_array[:5, :5]["CosineGreedy_matches"]

# plot scores as a matrix
from matplotlib import pyplot as plt
plt.figure(figsize=(26,26), dpi=150)
#plt.imshow(scores_array[:50, :50][type(similarity_measure).__name__ + "_score"], cmap="viridis")
plt.imshow(scores_array[type(similarity_measure).__name__ + "_score"], cmap="viridis")
plt.colorbar(shrink=0.7)
plt.title("Spectra similarities by " + type(similarity_measure).__name__)
plt.xlabel("Spectrum #ID")
plt.ylabel("Spectrum #ID")
plt.savefig("similarity_matrix.png", dpi=300)


# return scores for the top-10 candidates (Cosine score + number of matching peaks):
# first one is match to itself
best_matches = scores.scores_by_query(spectrums[5], name="CosineGreedy_score", sort=True)[:10]
print([x[1] for x in best_matches])
[x[0].get("smiles") for x in best_matches]

## to draw the stuctures
#from rdkit import Chem
#from rdkit.Chem import Draw
#for i, smiles in enumerate([x[0].get("smiles") for x in    best_matches]):
#    m = Chem.MolFromSmiles(smiles)
#    Draw.MolToFile(m, f"compound_{i}.png")

## apply a minimum number of matching peaks filter
min_match = plt.figure(figsize=(6,6), dpi=150)
plt.imshow(scores_array[:50, :50]["ModifiedCosine_score"] \
           * (scores_array[:50, :50]["ModifiedCosine_matches"] >= min_match), cmap="viridis")
plt.colorbar(shrink=0.7)
plt.title("Modified Cosine spectra similarities (min_match=5)")
plt.xlabel("Spectrum #ID")
plt.ylabel("Spectrum #ID")