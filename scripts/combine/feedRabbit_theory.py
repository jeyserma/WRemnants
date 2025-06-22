import argparse
import hist
import numpy as np
from combinetf2 import tensorwriter
import pickle
import lz4.frame
import combinetf2
import combinetf2.io_tools
from wremnants import theory_corrections
from wums import boostHistHelpers as hh
from utilities import common

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", default="./", help="output directory")
parser.add_argument("--outname", default="test_tensor", help="output file name")
parser.add_argument(
    "--postfix",
    default=None,
    type=str,
    help="Postfix to append on output file name",
)
parser.add_argument(
    "--sparse",
    default=False,
    action="store_true",
    help="Make sparse tensor",
)
parser.add_argument(
    "--symmetrizeAll",
    default=False,
    action="store_true",
    help="Make fully symmetric tensor",
)
parser.add_argument(
    "--skipMaskedChannels",
    default=False,
    action="store_true",
    help="Skip adding masked channels",
)
parser.add_argument(
    "--systematicType",
    choices=["log_normal", "normal"],
    default="log_normal",
    help="probability density for systematic variations",
)

args = parser.parse_args()


theory_input_file = f"{common.data_dir}/TheoryCorrections/scetlib_dyturboCT18ZVarsCorrZ.pkl.lz4"
with lz4.frame.open(theory_input_file, "rb") as f:
    data = pickle.load(f)

theory_corrs = [
    "scetlib_dyturbo",
    "scetlib_dyturboCT18ZVars",
    "scetlib_dyturboCT18Z_pdfas",
    #"powhegFOEW",
    #"pythiaew_ISR",
]
corr_helpers = theory_corrections.load_corr_helpers(
    "Z", theory_corrs, make_tensor=False, minnlo_ratio=False
)

# load in data histogram and covariance matrix
# TODO: fix hardcoding
data_input_file = "/ceph/submit/data/group/cms/store/user/lavezzo/alphaS//250617_histmaker//ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile/fitresults_asimov.hdf5"
fitresult, meta = combinetf2.io_tools.get_fitresult(
    data_input_file, result='asimov', meta=True
)
h_data = fitresult['physics_models']['Select helicitySig:slice(0,1)']['channels']['ch0_masked']['hist_postfit_inclusive'].get()[:,:,0] # helicity UL
h_data_cov = fitresult['physics_models']['Select helicitySig:slice(0,1)']['hist_postfit_inclusive_cov'].get()

# Build tensor
writer = tensorwriter.TensorWriter(
    sparse=args.sparse,
    systematic_type=args.systematicType,
)

# add data
writer.add_channel(h_data.axes, "ch0")
writer.add_data(h_data, "ch0")
writer.add_data_covariance(h_data_cov) # run with --externalCovariance

proc = "Z"
generator = "scetlib_dyturbo"
fname = f"{common.data_dir}/TheoryCorrections/{generator}Corr{proc[0]}.pkl.lz4"
h_sig = theory_corrections.load_corr_hist(fname, proc[0], f"{generator}_hist")[{'vars': 'pdf0'}].project("qT", "absY")
h_sig = hh.rebinHist(h_sig, 'qT', h_data.axes[0].edges)
h_sig = hh.rebinHist(h_sig, 'absY', h_data.axes[1].edges)

lumi = 16800
scale = lumi

# add signal
h_sig = h_sig * scale
writer.add_process(h_sig, "Zmumu", "ch0", signal=False)

# add systematic variations

corr_coeffs = theory_corrections.make_qcd_uncertainty_helper_by_helicity(
    is_z=True,
    filename=f"/ceph/submit/data/group/cms/store/user/lavezzo/alphaS//250617_gen/w_z_helicity_xsecs_scetlib_dyturboCorr_maxFiles_m1.hdf5",
    rebi_ptVgen=False,
    return_tensor=False,
)
corr_coeffs = corr_coeffs.project("ptVgen", "absYVgen", "vars", "helicity", "corr")
corr_coeffs = corr_coeffs[{'corr': 1}]
corr_coeffs = corr_coeffs[{'helicity': -1}]
corr_coeffs = corr_coeffs[{'vars': 'pythia_shower_kt'}]
corr_coeffs = hh.rebinHist(corr_coeffs, 'ptVgen', h_data.axes[0].edges)
corr_coeffs = hh.rebinHist(corr_coeffs, 'absYVgen', h_data.axes[1].edges)
writer.add_systematic(
    corr_coeffs,
    "pythia_shower_kt",
    "Zmumu",
    "ch0",
    mirror=True, 
    groups=["helicity_shower_kt", "angularCoeffs", "theory"],
)

print("Now at scetlib_dyturbo")
h = corr_helpers['Z']['scetlib_dyturbo']*scale
h = h.project("qT", "absY", "vars")
h = hh.rebinHist(h, 'qT', h_data.axes[0].edges)
h = hh.rebinHist(h, 'absY', h_data.axes[1].edges)
vars = [

    # gamma NP uncertainties
    ["omega_nu0.5","c_nu-0.1-omega_nu0.5", "scetlibNPgamma"],

    # correlated NP uncertainties
    ["Lambda20.25","Lambda2-0.25", "chargeVgenNP0scetlibNPZLambda2"],
    ["Lambda4.01","Lambda4.16", "chargeVgenNP0scetlibNPZLambda4"],
    ["Delta_Lambda2-0.02","Delta_Lambda20.02", "chargeVgenNP0scetlibNPZDelta_Lambda2"],
]
for var in vars:
    print(var)

    writer.add_systematic(
        [h[{'vars': var[1]}],  h[{'vars': var[0]}]],
        var[2],
        "Zmumu",
        "ch0",
        symmetrize='average', 
        groups=["resumTNP", "resum", "pTModeling", "theory"]
    )

vars = [
    # TNP
    ["gamma_cusp-1.","gamma_cusp1."],
    ["gamma_mu_q-1.","gamma_mu_q1."],
    ["gamma_nu-1.", "gamma_nu1."],
    ["h_qqV-1.", "h_qqV1."],
    ["s-1.", "s1."],
    ["b_qqV-0.5", "b_qqbarV-0.5"],
    ["b_qqS-0.5", "b_qqDS-0.5"],
    ["b_qqV0.5", "b_qqbarV0.5"],
    ["b_qqS0.5", "b_qqDS0.5"],
    ["b_qg-0.5", "b_qg0.5"],

]
for var in vars:
    print(var)

    writer.add_systematic(
        [h[{'vars': var[1]}],  h[{'vars': var[0]}]],
        "resumTNP_" + var[0].split('-')[0],
        "Zmumu",
        "ch0",
        symmetrize='average', 
        groups=["resumTNP", "resum", "pTModeling", "theory"]
    )

vars = [
    # transition FO scale uncertainties
    ["transition_points0.2_0.35_1.0", "transition_points0.2_0.75_1.0", "resumTransitionZ"],
    ["renorm_scale_pt20_envelope_Down", "renorm_scale_pt20_envelope_Up", "resumFOScaleZ"],
]
for var in vars:
    print(var)

    writer.add_systematic(
        [h[{'vars': var[1]}],  h[{'vars': var[0]}]],
        var[2],
        "Zmumu",
        "ch0",
        symmetrize='quadratic',
        groups=["resumTransitionFOScale", "resum", "pTModeling", "theory"],
    )

print("Now at scetlib_dyturboCT18ZVars")
h = corr_helpers['Z']['scetlib_dyturboCT18ZVars']*scale* 1 / 1.645
h = h.project("qT", "absY", "vars")
h = hh.rebinHist(h, 'qT', h_data.axes[0].edges)
h = hh.rebinHist(h, 'absY', h_data.axes[1].edges)
for ivar in range(1, len(h.axes[-1]), 2):
    print(ivar)
    h_syst_up = h[{'vars': ivar}]
    h_syst_down = h[{'vars': ivar + 1}]
    writer.add_systematic(
        [h_syst_up, h_syst_down],
        f"pdf{int((ivar+1)/2)}CT18Z",
        "Zmumu",
        "ch0",
        symmetrize='quadratic',
        groups=["pdfCT18Z", f"pdfCT18ZNoAlphaS", "theory"],
    )

print("Now at scetlib_dyturboCT18Z_pdfas")
h = corr_helpers['Z']['scetlib_dyturboCT18Z_pdfas']*scale
h = h.project("qT", "absY", "vars")
h = hh.rebinHist(h, 'qT', h_data.axes[0].edges)
h = hh.rebinHist(h, 'absY', h_data.axes[1].edges)
h_syst_up = h[{'vars': 2}]
h_syst_down = h[{'vars': 1}]
writer.add_systematic(
    [h_syst_up, h_syst_down],
    'pdfAlphaS',
    "Zmumu",
    "ch0",
    noi=True,
    constrained=False,
    symmetrize='average',
    #groups=["slopes", "slopes_background"],
)


print("Now at powhegFOEW")
theory_corrs = [
    "powhegFOEW",
]
corr_helpers = theory_corrections.load_corr_helpers(
    "Z", theory_corrs, make_tensor=False
)
h = corr_helpers['Z']['powhegFOEW']
h = h.project("absYVlhe", "weak")
h = hh.rebinHist(h, 'absYVlhe', [h_data.axes[1].edges[0], *h_data.axes[1].edges[2:]])
vals = h.values()
vals = np.vstack([vals[0, :], vals])
vals = np.broadcast_to(vals, (len(h_data.axes[0]), vals.shape[0], vals.shape[1]))
h = hist.Hist(
    h_data.axes[0],
    h_data.axes[1],
    h.axes[1],
)
h[...] = vals

writer.add_systematic(
    h[{'weak': 'weak_ps'}],
    'weak_ps',
    "Zmumu",
    "ch0",
    mirror=True,
    groups=[f"theory_ew_virtZ_scheme", "theory_ew", "theory"],
)
writer.add_systematic(
    h[{'weak': 'weak_aem'}],
    'weak_aem',
    "Zmumu",
    "ch0",
    mirror=True,
    groups=[f"theory_ew_virtZ_scheme", "theory_ew", "theory"],
)
writer.add_systematic(
    h[{'weak': 'weak_aem'}],
    'weak_aem',
    "Zmumu",
    "ch0",
    mirror=True,
    groups=[f"theory_ew_virtZ_corr", "theory_ew", "theory"],
)

# TODO need to re-run with consistent pT binning
# print("Now at pythiaew_ISR")
# theory_corrs = [
#     "pythiaew_ISR",
# ]
# corr_helpers = theory_corrections.load_corr_helpers(
#     "Z", theory_corrs, make_tensor=False
# )
# h = corr_helpers['Z']['pythiaew_ISR']
# print(h)
# h = h.project("ptVgen", "absYVgen", "systIdx")
# h = h[{"systIdx": 1}]
# h = hh.rebinHist(h, 'ptVgen', h_data.axes[0].edges)
# h = hh.rebinHist(h, 'absYVlhe', [h_data.axes[1].edges[0], *h_data.axes[1].edges[2:]])
# print("pythiaew_ISRCorr")
# h_syst = h[{'vars': var}]
# writer.add_systematic(
#     h_syst,
#     "pythiaew_ISRCorr",
#     "Zmumu",
#     "ch0",
#     mirror=True,
#     kfactor=2.0,
#     groups=[f"theory_ew_pythiaew_ISR", "theory_ew", "theory"],
# )


directory = args.output
if directory == "":
    directory = "./"
filename = args.outname
if args.postfix:
    filename += f"_{args.postfix}"
writer.write(outfolder=directory, outfilename=filename)
print(f"Written to {directory}{filename}.hdf5")