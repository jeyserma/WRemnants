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
parser.add_argument(
    "infile",
    type=str,
    help="Input unfolded fit result",
)
parser.add_argument(
    "--predFile",
    type=str,
    default=f"{common.data_dir}/TheoryCorrections/scetlib_dyturboCorrZ.pkl.lz4",
    help="File containing theory predictions for the Z ptVgen and absYVgen",
)
parser.add_argument(
    "--predGenerator",
    type=str,
    default=f"scetlib_dyturbo",
    help="Generator used for the predictions. Expect the prediction histogram to be named as <generator>_hist",
)
parser.add_argument("-o", "--output", default="./", help="output directory")
parser.add_argument("--outname", default="carrot", help="output file name")
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


# load in data histogram and covariance matrix
# TODO: test difference between files
fitresult, meta = combinetf2.io_tools.get_fitresult(
    args.infile, result='asimov', meta=True
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

# scaling for background and variations
lumi = 16800
scale = lumi

# add background
h_sig = theory_corrections.load_corr_hist(args.predFile, "Z", f"{args.predGenerator}_hist")[{'vars': 'pdf0'}].project("qT", "absY")
h_sig = hh.rebinHist(h_sig, 'qT', h_data.axes[0].edges)
h_sig = hh.rebinHist(h_sig, 'absY', h_data.axes[1].edges)
h_sig = h_sig * scale
writer.add_process(h_sig, "Zmumu", "ch0", signal=False)

# add systematic variations

# TODO the corrections (in theory_corections.py) are implemented as a ratio to the UL. Does that mean that there is no uncert for UL?
# corr_coeffs = theory_corrections.make_qcd_uncertainty_helper_by_helicity(
#     is_z=True,
#     filename=f"/ceph/submit/data/group/cms/store/user/lavezzo/alphaS//250617_gen/w_z_helicity_xsecs_scetlib_dyturboCorr_maxFiles_m1.hdf5",
#     rebi_ptVgen=False,
#     return_tensor=False,
# )
# corr_coeffs = corr_coeffs*scale
# corr_coeffs = corr_coeffs.project("ptVgen", "absYVgen", "vars", "helicity", "corr")
# corr_coeffs = corr_coeffs[{'corr': 1}]
# corr_coeffs = corr_coeffs[{'helicity': -1.0j}]
# corr_coeffs = corr_coeffs[{'vars': 'pythia_shower_kt'}]
# corr_coeffs = hh.rebinHist(corr_coeffs, 'ptVgen', h_data.axes[0].edges)
# corr_coeffs = hh.rebinHist(corr_coeffs, 'absYVgen', h_data.axes[1].edges)
# writer.add_systematic(
#     corr_coeffs,
#     "pythia_shower_kt",
#     "Zmumu",
#     "ch0",
#     mirror=True, 
#     groups=["helicity_shower_kt", "angularCoeffs", "theory"],
# )

print("Now at scetlib_dyturbo")
theory_corrs = [
    "scetlib_dyturbo",
]
corr_helpers = theory_corrections.load_corr_helpers(
    "Z", theory_corrs, make_tensor=False, minnlo_ratio=False
)
h = corr_helpers['Z']['scetlib_dyturbo']*scale
h = h.project("qT", "absY", "vars")
h = hh.rebinHist(h, 'qT', h_data.axes[0].edges)
h = hh.rebinHist(h, 'absY', h_data.axes[1].edges)

# correlated NP uncertainties
corr_NP_uncs = [
    ["Lambda20.25","Lambda2-0.25", "chargeVgenNP0scetlibNPZLambda2"],
    ["Lambda4.01","Lambda4.16", "chargeVgenNP0scetlibNPZLambda4"],
    ["Delta_Lambda2-0.02","Delta_Lambda20.02", "chargeVgenNP0scetlibNPZDelta_Lambda2"],
]
for var in corr_NP_uncs:
    print(var)

    writer.add_systematic(
        [h[{'vars': var[1]}],  h[{'vars': var[0]}]],
        var[2],
        "Zmumu",
        "ch0",
        symmetrize='average', 
        groups=["resumTNP", "resum", "pTModeling", "theory"]
    )

# gamma NP uncertainties
gamma_NP_uncs = [
    ["omega_nu0.5","c_nu-0.1-omega_nu0.5", "scetlibNPgamma"],
]
for var in gamma_NP_uncs:
    print(var)

    writer.add_systematic(
        [h[{'vars': var[1]}],  h[{'vars': var[0]}]],
        var[2],
        "Zmumu",
        "ch0",
        symmetrize='average', 
        groups=["resumTNP", "resum", "pTModeling", "theory"]
    )

# TNP
TNP_uncs = [
    ["gamma_cusp-1.","gamma_cusp1."],
    ["gamma_mu_q-1.","gamma_mu_q1."],
    ["gamma_nu-1.", "gamma_nu1."],
    ["h_qqV-1.", "h_qqV1."],
    ["s-1.", "s1."],
    ["b_qqV-0.5", "b_qqV0.5"],
    ["b_qqbarV-0.5", "b_qqbarV0.5"],
    ["b_qqS-0.5", "b_qqS0.5"],
    ["b_qqDS-0.5", "b_qqDS0.5"],
    ["b_qg-0.5", "b_qg0.5"],
]
for var in TNP_uncs:
    print(var)

    writer.add_systematic(
        [h[{'vars': var[1]}],  h[{'vars': var[0]}]],
        "resumTNP_" + var[0].split('-')[0],
        "Zmumu",
        "ch0",
        symmetrize='average', 
        groups=["resumTNP", "resum", "pTModeling", "theory"]
    )

# transition FO scale uncertainties
transition_FO_uncs = [
    ["transition_points0.2_0.35_1.0", "transition_points0.2_0.75_1.0", "resumTransitionZ"],
    ["renorm_scale_pt20_envelope_Down", "renorm_scale_pt20_envelope_Up", "resumFOScaleZ"],
]
for var in transition_FO_uncs:
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
theory_corrs = [
    "scetlib_dyturboCT18ZVars",
]
corr_helpers = theory_corrections.load_corr_helpers(
    "Z", theory_corrs, make_tensor=False, minnlo_ratio=False
)
h = corr_helpers['Z']['scetlib_dyturboCT18ZVars']*scale
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
        kfactor=1/1.645,
        groups=["pdfCT18Z", f"pdfCT18ZNoAlphaS", "theory"],
    )

print("Now at scetlib_dyturboCT18Z_pdfas")
theory_corrs = [
    "scetlib_dyturboCT18Z_pdfas",
]
corr_helpers = theory_corrections.load_corr_helpers(
    "Z", theory_corrs, make_tensor=False, minnlo_ratio=False
)
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
    kfactor=1.5/2.0
)

print("Now at MSHT20mcrangeCorrZ")
theory_corrs = [
    "scetlib_dyturboMSHT20mcrange",
    "scetlib_dyturboMSHT20mbrange"
]
corr_helpers = theory_corrections.load_corr_helpers(
    "Z", theory_corrs, make_tensor=False, minnlo_ratio=False
)
h = corr_helpers['Z']['scetlib_dyturboMSHT20mbrange']*scale  
h = h.project("qT", "absY", "vars")
h = hh.rebinHist(h, 'qT', h_data.axes[0].edges)
h = hh.rebinHist(h, 'absY', h_data.axes[1].edges)  
writer.add_systematic(
    [h[{'vars': 1}], h[{'vars': -1}]],
    "pdfMSHT20mbrange",
    "Zmumu",
    "ch0",
    symmetrize='quadratic',
    groups=["bcQuarkMass", "pTModeling", "theory"],
)
h = corr_helpers['Z']['scetlib_dyturboMSHT20mcrange']*scale  
h = h.project("qT", "absY", "vars")
h = hh.rebinHist(h, 'qT', h_data.axes[0].edges)
h = hh.rebinHist(h, 'absY', h_data.axes[1].edges) 
writer.add_systematic(
    [h[{'vars': 1}], h[{'vars': -1}]],
    "pdfMSHT20mcrange",
    "Zmumu",
    "ch0",
    symmetrize='quadratic',
    groups=["bcQuarkMass", "pTModeling", "theory"],
)

# TODO fix this
print("Now at powhegFOEW")
theory_corrs = [
    "powhegFOEW",
]
corr_helpers = theory_corrections.load_corr_helpers(
    "Z", theory_corrs, make_tensor=False
)
h = corr_helpers['Z']['powhegFOEW']
print(h[{'weak': 'weak_default'}])
print(h_data)
h = h.project("absYVlhe", "weak")[{'weak': 'weak_ps'}]
print(h)
h = hh.rebinHist(h, 'absYVlhe', [0, 0.5, 1, 1.5, 2.5])
print(h_sig.project("absY"))
denom = h_sig.project("absY")
denom = hh.rebinHist(denom, 'absY', [0, 0.5, 1, 1.5, 2.5])
print(denom.values()/h.values())
exit()

vals = h.values()
vals = np.vstack([vals[0, :], vals])
vals = np.broadcast_to(vals, (len(h_data.axes[0]), vals.shape[0], vals.shape[1]))
denom = np.broadcast_to(denom, vals.shape[0:2])
vals = np.array([vals[:,:,i]/denom for i in range(vals.shape[2])])
print(vals)
# vals = (1+vals) * h_data.values().reshape(*h_data.values().shape, 1)
h = hist.Hist(
    h.axes[1],      # vars
    h_data.axes[0], # pT
    h_data.axes[1], # absY
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
    h[{'weak': 'weak_default'}],
    'weak_default',
    "Zmumu",
    "ch0",
    mirror=True,
    groups=[f"theory_ew_virtZ_corr", "theory_ew", "theory"],
)

print("Now at pythiaew_ISR")
theory_corrs = [
    "pythiaew_ISR",
]
fname = f"{common.data_dir}/TheoryCorrections/pythiaew_ISRCorr{proc[0]}.pkl.lz4"
h = theory_corrections.load_corr_hist(fname, proc[0], f"{theory_corrs[0]}_minnlo_ratio") # not actually minnlo ratio
h = h.project("ptVgen", "absYVgen", "systIdx")
h = h[{'systIdx': 1}] 
h = hh.rebinHist(h, 'ptVgen', h_data.axes[0].edges)
h = hh.rebinHist(h, 'absYVgen', h_data.axes[1].edges)
h[...] = h.values() * h_sig.values() # TODO does this make sense?
writer.add_systematic(
    h,
    "pythiaew_ISRCorr0",
    "Zmumu",
    "ch0",
    mirror=True,
    kfactor=2.0,
    groups=[f"theory_ew_pythiaew_ISR", "theory_ew", "theory"],
)




directory = args.output
if directory == "":
    directory = "./"
filename = args.outname
if args.postfix:
    filename += f"_{args.postfix}"
writer.write(outfolder=directory, outfilename=filename)
print(f"Written to {directory}{filename}.hdf5")