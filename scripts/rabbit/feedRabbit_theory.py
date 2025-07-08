import argparse
import copy

import h5py
import hist
import numpy as np

import rabbit
import rabbit.io_tools
from rabbit import tensorwriter
from utilities import common
from utilities.io_tools import input_tools
from wremnants import theory_corrections
from wums import boostHistHelpers as hh


class HistFormatter:

    def __init__(self, boson, channel, hist_ref, lumi=16800, scale=1.0):

        supported_channels = ["chSigmaUL", "chAis"]
        supported_bosons = ["Z"]
        if channel not in supported_channels:
            raise Exception(f"Channel {channel} is not supported.")
        if boson not in supported_bosons:
            raise Exception(f"Boson {boson} is not supported.")

        self.boson = boson
        self.channel = channel
        self.hist_ref = hist_ref
        self.norm = lumi * scale

        self.pt_bins = self.hist_ref.axes["ptVGen"].edges
        self.y_bins = self.hist_ref.axes["absYVGen"].edges

    def __call__(self, h, rebin_pt=True, rebin_y=True, normalize=True):
        """
        Apply mass, charge selections to input histogram.
        Re-order axes to match the predictions (pt, absY, ...).
        Optionally, re-bin pt and y.
        Returns histogram with changes applied.
        """

        h = copy.deepcopy(h)  # pass by value

        axes_names = list(h.axes.name)

        h = self.apply_selections(h)

        if "ptVgen" in axes_names:
            pt_axis_name = "ptVgen"
        elif "ptVGen" in axes_names:
            pt_axis_name = "ptVGen"
        elif "qT" in axes_names:
            pt_axis_name = "qT"
        else:
            raise ValueError("Neither 'ptVgen' nor 'ptVGen' found in histogram axes.")

        if "absYVgen" in axes_names:
            absY_axis_name = "absYVgen"
        elif "absYVGen" in axes_names:
            absY_axis_name = "absYVGen"
        elif "absY" in axes_names:
            absY_axis_name = "absY"
        else:
            raise ValueError(
                "Neither 'absYVgen' nor 'absYVGen' found in histogram axes."
            )

        # optionally, rebin pt and y, and scale
        if rebin_pt:
            h = hh.rebinHist(h, pt_axis_name, self.pt_bins)
        if rebin_y:
            h = hh.rebinHist(h, absY_axis_name, self.y_bins)
        if normalize:
            h *= self.norm

        # re-order remaining axes to match expected output
        remaining_axes = list(h.axes.name)
        remaining_axes.remove(pt_axis_name)
        remaining_axes.remove(absY_axis_name)
        h = h.project(pt_axis_name, absY_axis_name, *remaining_axes)

        # calculate Ai's from helicities, or select the sigmaUL
        if self.channel == "chAis":
            h = calculate_ais_from_helicities_hist(h)
        elif self.channel == "chSigmaUL":
            if "helicity" in axes_names:
                h = h[{"helicity": -1.0j}]

        return h

    def apply_selections(self, h):
        axes_names = list(h.axes.name)
        if self.boson == "Z":

            # mass selection
            if "massVgen" in axes_names:
                h = h[{"massVgen": 90.0j}]
            elif "Q" in axes_names:
                h = h[{"Q": 90.0j}]

            # charge selection
            if "chargeVgen" in axes_names:
                h = h[{"chargeVgen": 0.0j}]
            elif "charge" in axes_names:
                h = h[{"charge": 0.0j}]

        return h


def calculate_ais_from_helicities_hist(h_hels):
    """
    Calculate A_i histogram from helicities histogram.
    Assuming as input a histogram with the helicities,
    with the first element being the sigma_UL.
    Returns a the A_i's histogram.
    """

    h_hels = copy.deepcopy(h_hels)  # pass by value

    # re-order histogram axes such that the last one is the helicity
    if h_hels.axes[-1].name != "helicity":
        axes = list(h_hels.axes.name)
        axes.remove("helicity")
        axes.append("helicity")
        h_hels = h_hels.project(*axes)

    vals = copy.copy(h_hels.values())
    vars = copy.copy(
        h_hels.variances()
    )  # need to copy in case the storage is hist.storage.Double

    # A_i = sigma_i/sigma_UL
    num = vals[..., 1:]
    den = vals[..., 0][..., np.newaxis]
    vals[..., 1:] = np.where(den == 0, 0.0, num / den)

    # treat these as poisson uncorrelated rv's
    if np.any(vars):

        vars[..., 1:] = np.where(
            vals[..., 0][..., np.newaxis] == 0,
            0.0,
            den**-2 * vars[..., 1:] + num**2 * den**-4 * vars[..., 0][..., np.newaxis],
        )

    h_out = hist.Hist(
        *h_hels.axes[:-1],
        hist.axis.Integer(0, 8, name="ais"),
        storage=hist.storage.Weight() if np.any(vars) else hist.storage.Double(),
    )
    h_out.values()[...] = vals[..., 1:]
    if np.any(vars):
        h_out.variances()[...] = vars[..., 1:]

    return h_out


parser = argparse.ArgumentParser()
parser.add_argument(
    "infile",
    type=str,
    help="Input unfolded fit result",
)
parser.add_argument(
    "--fitresultModel", type=str, default="Select helicitySig:slice(0,1)"
)
parser.add_argument(
    "--predGenerator",
    type=str,
    default=f"scetlib_dyturbo",
    help="Generator used for the sigmaUL predictions. "
    "Expect the prediction file to be named as <generator>CorrZ.pkl.lz4. "
    "Expect the prediction histogram to be named as <generator>_hist",
)
parser.add_argument(
    "--predAiFile",
    type=str,
    default=f"{common.data_dir}/w_z_helicity_xsecs_scetlib_dyturboCorr_maxFiles_m1_alphaSunfoldingBinning_helicity.hdf5",
    help="Gen file used for the Ai predictions."
    "Will be stitched with the --predGenerator file.",
)
parser.add_argument(
    "--angularCoeffs",
    action="store_true",
    default=False,
    help="Fit the angular coefficients."
    "Predictions of the Ai's from --predAiFile."
    "Predictions of the sigma_UL from infile.",
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

# Build tensor
writer = tensorwriter.TensorWriter(
    sparse=args.sparse,
    systematic_type="normal" if args.angularCoeffs else args.systematicType,
    allow_negative_expectation=args.angularCoeffs,
)

# load in data histogram and covariance matrix
fitresult, meta = rabbit.io_tools.get_fitresult(args.infile, result="asimov", meta=True)
if args.angularCoeffs:
    try:
        h_data_ai = fitresult["physics_models"][args.fitresultModel]["channels"][
            "AngularCoefficients ch0_masked ptVGen:rebin(0,3,6,9,12,16,20,24,28,33,44)_ch0_masked"
        ]["hist_postfit_inclusive"].get()
        h_data = fitresult["physics_models"][args.fitresultModel]["channels"][
            "Select helicitySig:slice(0,1)_ch0_masked"
        ]["hist_postfit_inclusive"].get()[:, :, 0]
        h_data_cov = fitresult["physics_models"][args.fitresultModel][
            "hist_postfit_inclusive_cov"
        ].get()
    except KeyError:
        raise KeyError(
            f"Couldn't find {args.fitresultModel}. Available physics models are: {",".join(fitresult['physics_models'].keys())}"
        )
    writer.add_channel(h_data_ai.axes, "chAis")
    writer.add_data(h_data_ai, "chAis")
    writer.add_channel(h_data.axes, "chSigmaUL")
    writer.add_data(h_data, "chSigmaUL")
else:
    h_data = fitresult["physics_models"][args.fitresultModel]["channels"]["ch0_masked"][
        "hist_postfit_inclusive"
    ].get()[
        :, :, 0
    ]  # grabbing the unpolarized term
    h_data_cov = fitresult["physics_models"][args.fitresultModel][
        "hist_postfit_inclusive_cov"
    ].get()
    writer.add_channel(h_data.axes, "chSigmaUL")
    writer.add_data(h_data, "chSigmaUL")
writer.add_data_covariance(h_data_cov)  # N.B: run fit with --externalCovariance

# now that we have loaded in data, we can define our histgram formatters
chSigmaULFormatter = HistFormatter("Z", "chSigmaUL", h_data)
if args.angularCoeffs:
    chAisFormatter = HistFormatter("Z", "chAis", h_data_ai)

# add background, Z->mumu
h_sig_sigmaUL = theory_corrections.load_corr_hist(
    f"{common.data_dir}/TheoryCorrections/{args.predGenerator}CorrZ.pkl.lz4",
    "Z",
    f"{args.predGenerator}_hist",
)
h_sig_sigmaUL = h_sig_sigmaUL[{"vars": "pdf0"}]  # select baseline variation
h_sig_sigmaUL = chSigmaULFormatter(h_sig_sigmaUL)
writer.add_process(h_sig_sigmaUL, "Zmumu", "chSigmaUL", signal=False)

# if set, load in the helicity cross sections predictions from MINNLO
if args.angularCoeffs:
    with h5py.File(args.predAiFile, "r") as ff:
        inputs = input_tools.load_results_h5py(ff)
        h_sig_hels = inputs["Z"]
    h_sig_hels = h_sig_hels[{"muRfact": 1.0j}][
        {"muFfact": 1.0j}
    ]  # select baseline variation
    h_sig_ais = chAisFormatter(h_sig_hels)
    writer.add_process(h_sig_ais, "Zmumu", "chAis", signal=False)

# add systematic variations

# alphaS variation
if args.predGenerator == "scetlib_dyturbo":
    theory_corrs = [
        "scetlib_dyturboCT18Z_pdfas",
    ]
    corr_helpers = theory_corrections.load_corr_helpers(
        "Z", theory_corrs, make_tensor=False, minnlo_ratio=False
    )
    h = corr_helpers["Z"]["scetlib_dyturboCT18Z_pdfas"]
    h = chSigmaULFormatter(h)
    h_syst_down = h[{"vars": 1}]
    h_syst_up = h[{"vars": 2}]

    if args.angularCoeffs:
        with h5py.File(
            args.predAiFile.replace("w_z_helicity_xsecs", "w_z_gen_dists"), "r"
        ) as ff:
            inputs = input_tools.load_results_h5py(ff)
            h_hels = inputs["ZmumuPostVFP"]["output"][
                "nominal_gen_helicity_nominal_gen_pdfCT18ZalphaS002"
            ].get()
        h_ais = chAisFormatter(h_hels)
        h_syst_up_ais = h_ais[{"alphasVar": "as0120"}]
        h_syst_down_ais = h_ais[{"alphasVar": "as0116"}]
        writer.add_systematic(
            [h_syst_up_ais, h_syst_down_ais],
            "pdfAlphaS",
            "Zmumu",
            "chAis",
            noi=True,
            constrained=False,
            symmetrize="average",
            kfactor=1.5 / 2.0,
        )

    writer.add_systematic(
        [h_syst_up, h_syst_down],
        "pdfAlphaS",
        "Zmumu",
        "chSigmaUL",
        noi=True,
        constrained=False,
        symmetrize="average",
        kfactor=1.5 / 2.0,
    )
elif args.predGenerator == "scetlib_nnlojetN4p0LLN3LO":
    # TODO need to check this more carefully
    # alphaS variation N4p0LLN3LO
    fname = (
        f"{common.data_dir}/TheoryCorrections//{args.predGenerator}_pdfasCorrZ.pkl.lz4"
    )
    corrh = theory_corrections.load_corr_hist(
        fname, "Z", f"{args.predGenerator}_pdfas_minnlo_ratio"
    ).project("qT", "absY", "vars")
    minnlo_ref = theory_corrections.load_corr_hist(
        fname, "Z", "minnlo_ref_hist"
    ).project("qT", "absY")
    vals = np.einsum("ijk,ij->ijk", corrh.values(), minnlo_ref.values())
    corrh[...] = vals * 16800
    corrh = hh.rebinHist(corrh, "qT", h_data.axes[0].edges)
    corrh = hh.rebinHist(corrh, "absY", h_data.axes[1].edges)
    h_syst_up = corrh[{"vars": 2}]
    h_syst_down = corrh[{"vars": 1}]
    writer.add_systematic(
        [h_syst_up, h_syst_down],
        "pdfAlphaS",
        "Zmumu",
        "chSigmaUL",
        noi=True,
        constrained=False,
        symmetrize="average",
        # kfactor=1.5/2.0
    )
else:
    raise Exception("No valid configuration found for alphaS variation.")

if args.angularCoeffs:

    qcd_helper = theory_corrections.make_qcd_uncertainty_helper_by_helicity(
        is_z=True,
        filename=args.predAiFile,
        rebi_ptVgen=h_data_ai.axes["ptVGen"].edges,
        rebi_absYVgen=h_data_ai.axes["absYVGen"].edges,
        return_tensor=False,
    )

    qcd_helper_ais = chAisFormatter(
        qcd_helper, rebin_pt=False, rebin_y=False, normalize=False
    )
    qcd_helper_ais = hh.divideHists(
        qcd_helper_ais[{"corr": 1}], qcd_helper_ais[{"corr": 0}]
    )
    qcd_helper_ais = qcd_helper_ais.project("ptVgen", "absYVgen", "ais", "vars")
    h_sig_bd = hh.broadcastSystHist(
        h_sig_ais, qcd_helper_ais, by_ax_name=False, flow=False
    )
    qcd_helper_ais = hh.multiplyHists(qcd_helper_ais, h_sig_bd)

    # pythia showering uncertainties
    print("Now at pythia_shower_kt")
    writer.add_systematic(
        qcd_helper_ais[{"vars": "pythia_shower_kt"}],
        "pythia_shower_kt",
        "Zmumu",
        "chAis",
        mirror=True,
        groups=["helicity_shower_kt", "angularCoeffs", "theory"],
    )

    # QCD scales
    print("Now at QCD scales")

    # prepare fine binning hists
    fine_pt_binning = qcd_helper.axes["ptVgen"].edges
    nptfine = len(fine_pt_binning) - 1
    scale_inclusive = np.sqrt((nptfine - 1) / nptfine)

    for hel in range(0, 7 + 1):  # no correction on sigma_UL

        # fine binning
        qcd_scales_hel_up = qcd_helper_ais[{"vars": f"helicity_{hel}_Up"}].project(
            "ptVgen", "absYVgen", "ais"
        )
        qcd_scales_hel_down = qcd_helper_ais[{"vars": f"helicity_{hel}_Down"}].project(
            "ptVgen", "absYVgen", "ais"
        )

        for bin in range(len(fine_pt_binning) - 1):

            ptl = fine_pt_binning[bin]
            pth = fine_pt_binning[bin + 1]

            qcd_scales_hel_pt_up = copy.deepcopy(qcd_scales_hel_up)
            qcd_scales_hel_pt_up.values()[...] = h_sig_ais.values()
            qcd_scales_hel_pt_up.values()[bin, ...] = qcd_scales_hel_up[
                {"ptVgen": bin}
            ].values()

            qcd_scales_hel_pt_down = copy.deepcopy(qcd_scales_hel_down)
            qcd_scales_hel_pt_down.values()[...] = h_sig_ais.values()
            qcd_scales_hel_pt_down.values()[bin, ...] = qcd_scales_hel_up[
                {"ptVgen": bin}
            ].values()

            writer.add_systematic(
                [qcd_scales_hel_pt_up, qcd_scales_hel_pt_down],
                f"QCDscaleZfine_Pt{ptl}_{pth}helicity_{hel}",
                "Zmumu",
                "chAis",
                symmetrize="quadratic",
                groups=["QCDScaleZMiNNLO", "QCDscale", "angularCoeffs", "theory"],
            )

        # inclusive
        inclusive_pt_binning = [
            qcd_scales_hel_up.axes["ptVgen"].edges[0],
            qcd_scales_hel_up.axes["ptVgen"].edges[-1],
        ]
        qcd_scales_hel_int_up = hh.rebinHist(
            qcd_scales_hel_up, "ptVgen", inclusive_pt_binning
        )
        qcd_scales_hel_int_down = hh.rebinHist(
            qcd_scales_hel_down, "ptVgen", inclusive_pt_binning
        )

        for bin in range(len(fine_pt_binning) - 1):
            qcd_scales_hel_up.values()[bin, ...] = qcd_scales_hel_int_up.values()
            qcd_scales_hel_down.values()[bin, ...] = qcd_scales_hel_int_down.values()

        writer.add_systematic(
            [qcd_scales_hel_up, qcd_scales_hel_down],
            f"QCDscaleZinclusive_Pt{inclusive_pt_binning[0]}_{inclusive_pt_binning[1]}helicity_{hel}",
            "Zmumu",
            "chAis",
            kfactor=scale_inclusive,
            symmetrize="quadratic",
            groups=["QCDScaleZMiNNLO", "QCDscale", "angularCoeffs", "theory"],
        )


print(f"Now at variations from {args.predGenerator}")
theory_corrs = [
    args.predGenerator,
]
corr_helpers = theory_corrections.load_corr_helpers(
    "Z", theory_corrs, make_tensor=False, minnlo_ratio=False
)
h = corr_helpers["Z"][args.predGenerator]
h = chSigmaULFormatter(h)

# correlated NP uncertainties
corr_NP_uncs = [
    ["Lambda20.25", "Lambda2-0.25", "chargeVgenNP0scetlibNPZLambda2"],
    ["Lambda4.16", "Lambda4.01", "chargeVgenNP0scetlibNPZLambda4"],
    ["Delta_Lambda20.02", "Delta_Lambda2-0.02", "chargeVgenNP0scetlibNPZDelta_Lambda2"],
]
for var in corr_NP_uncs:
    writer.add_systematic(
        [h[{"vars": var[0]}], h[{"vars": var[1]}]],
        var[2],
        "Zmumu",
        "chSigmaUL",
        symmetrize="average",
        groups=["resumNonpert", "resum", "pTModeling", "theory"],
    )

# gamma NP uncertainties
gamma_NP_uncs = [
    ["omega_nu0.5", "c_nu-0.1-omega_nu0.5", "scetlibNPgamma"],
]
for var in gamma_NP_uncs:
    writer.add_systematic(
        [h[{"vars": var[0]}], h[{"vars": var[1]}]],
        var[2],
        "Zmumu",
        "chSigmaUL",
        symmetrize="average",
        groups=["resumTNP", "resum", "pTModeling", "theory"],
    )

# TNP
TNP_uncs = [
    ["gamma_cusp1.", "gamma_cusp-1."],
    ["gamma_mu_q1.", "gamma_mu_q-1."],
    ["gamma_nu1.", "gamma_nu-1."],
    ["h_qqV1.", "h_qqV-1."],
    ["s1.", "s-1."],
    ["b_qqV0.5", "b_qqV-0.5"],
    ["b_qqbarV0.5", "b_qqbarV-0.5"],
    ["b_qqS0.5", "b_qqS-0.5"],
    ["b_qqDS0.5", "b_qqDS-0.5"],
    ["b_qg0.5", "b_qg-0.5"],
]
for var in TNP_uncs:
    var_name = "resumTNP_" + var[1].split("-")[0]
    writer.add_systematic(
        [h[{"vars": var[0]}], h[{"vars": var[1]}]],
        var_name,
        "Zmumu",
        "chSigmaUL",
        symmetrize="average",
        groups=["resumTNP", "resum", "pTModeling", "theory"],
    )

# transition FO scale uncertainties
transition_FO_uncs = [
    [
        "transition_points0.2_0.75_1.0",
        "transition_points0.2_0.35_1.0",
        "resumTransitionZ",
    ],
    [
        "renorm_scale_pt20_envelope_Up",
        "renorm_scale_pt20_envelope_Down",
        "resumFOScaleZ",
    ],
]
for var in transition_FO_uncs:
    writer.add_systematic(
        [h[{"vars": var[0]}], h[{"vars": var[1]}]],
        var[2],
        "Zmumu",
        "chSigmaUL",
        symmetrize="quadratic",
        groups=["resumTransitionFOScale", "resum", "pTModeling", "theory"],
    )

# PDF uncertainties
print("Now at PDF variations")
if args.angularCoeffs:
    # for Ai's, we have MINNLO, so use it for sigmaUL + Ai's to be consistent
    with h5py.File(
        args.predAiFile.replace("w_z_helicity_xsecs", "w_z_gen_dists"), "r"
    ) as ff:
        inputs = input_tools.load_results_h5py(ff)
        h = inputs["ZmumuPostVFP"]["output"]["nominal_gen_helicity_pdfCT18Z"].get()
    h_sigmaUL = chSigmaULFormatter(h)
    h_ais = chAisFormatter(h)
    for ivar in range(1, len(h.axes[-1]), 2):
        # sigmaUL
        h_syst_up = h_sigmaUL[{"pdfVar": ivar + 1}]
        h_syst_up = hh.divideHists(h_syst_up, h_sigmaUL[{"pdfVar": "pdf0CT18Z"}])
        h_syst_up = hh.multiplyHists(h_syst_up, h_sig_sigmaUL, flow=False)
        h_syst_down = h_sigmaUL[{"pdfVar": ivar}]
        h_syst_down = hh.divideHists(h_syst_down, h_sigmaUL[{"pdfVar": "pdf0CT18Z"}])
        h_syst_down = hh.multiplyHists(h_syst_down, h_sig_sigmaUL, flow=False)
        writer.add_systematic(
            [h_syst_up, h_syst_down],
            f"pdf{int((ivar+1)/2)}CT18Z",
            "Zmumu",
            "chSigmaUL",
            symmetrize="quadratic",
            kfactor=1 / 1.645,
            groups=["pdfCT18Z", f"pdfCT18ZNoAlphaS", "theory"],
        )
        # Ai's
        h_syst_up = h_ais[{"pdfVar": ivar + 1}]
        h_syst_up = hh.divideHists(h_syst_up, h_ais[{"pdfVar": "pdf0CT18Z"}])
        h_syst_up = hh.multiplyHists(h_syst_up, h_sig_ais, flow=False)
        h_syst_down = h_ais[{"pdfVar": ivar}]
        h_syst_down = hh.divideHists(h_syst_down, h_ais[{"pdfVar": "pdf0CT18Z"}])
        h_syst_down = hh.multiplyHists(h_syst_down, h_sig_ais, flow=False)
        writer.add_systematic(
            [h_syst_up, h_syst_down],
            f"pdf{int((ivar+1)/2)}CT18Z",
            "Zmumu",
            "chAis",
            symmetrize="quadratic",
            kfactor=1 / 1.645,
            groups=["pdfCT18Z", f"pdfCT18ZNoAlphaS", "theory"],
        )
else:
    # for sigmaUL, scetlib+dyturbo has the latest & greatest PDFs
    theory_corrs = [
        "scetlib_dyturboCT18ZVars",
    ]
    corr_helpers = theory_corrections.load_corr_helpers(
        "Z", theory_corrs, make_tensor=False, minnlo_ratio=False
    )
    h = corr_helpers["Z"]["scetlib_dyturboCT18ZVars"]
    h = chSigmaULFormatter(h)
    for ivar in range(1, len(h.axes[-1]), 2):
        h_syst_up = h[{"vars": ivar + 1}]
        h_syst_down = h[{"vars": ivar}]
        writer.add_systematic(
            [h_syst_up, h_syst_down],
            f"pdf{int((ivar+1)/2)}CT18Z",
            "Zmumu",
            "chSigmaUL",
            symmetrize="quadratic",
            kfactor=1 / 1.645,
            groups=["pdfCT18Z", f"pdfCT18ZNoAlphaS", "theory"],
        )

# quark mass effects
print("Now at MSHT20mcrangeCorrZ")
theory_corrs = ["scetlib_dyturboMSHT20mcrange", "scetlib_dyturboMSHT20mbrange"]
corr_helpers = theory_corrections.load_corr_helpers(
    "Z", theory_corrs, make_tensor=False, minnlo_ratio=False
)
h = corr_helpers["Z"]["scetlib_dyturboMSHT20mbrange"]
h = chSigmaULFormatter(h)
varup = hh.divideHists(h[{"vars": -1}], h[{"vars": 0}])
varup = hh.multiplyHists(h_sig_sigmaUL, varup)
vardown = hh.divideHists(h[{"vars": 1}], h[{"vars": 0}])
vardown = hh.multiplyHists(h_sig_sigmaUL, vardown)
writer.add_systematic(
    [varup, vardown],
    "pdfMSHT20mbrange",
    "Zmumu",
    "chSigmaUL",
    symmetrize="quadratic",
    groups=["bcQuarkMass", "pTModeling", "theory"],
)
h = corr_helpers["Z"]["scetlib_dyturboMSHT20mcrange"]
h = chSigmaULFormatter(h)
varup = hh.divideHists(h[{"vars": -1}], h[{"vars": 0}])
varup = hh.multiplyHists(h_sig_sigmaUL, varup)
vardown = hh.divideHists(h[{"vars": 1}], h[{"vars": 0}])
vardown = hh.multiplyHists(h_sig_sigmaUL, vardown)
writer.add_systematic(
    [varup, vardown],
    "pdfMSHT20mcrange",
    "Zmumu",
    "chSigmaUL",
    symmetrize="quadratic",
    groups=["bcQuarkMass", "pTModeling", "theory"],
)

# ISR corrections
print("Now at pythiaew_ISR")
theory_corrs = [
    "pythiaew_ISR",
]
fname = f"{common.data_dir}/TheoryCorrections/pythiaew_ISRCorrZ.pkl.lz4"
h = theory_corrections.load_corr_hist(
    fname, "Z", f"{theory_corrs[0]}_minnlo_ratio"
)  # not actually minnlo ratio
h = h[{"systIdx": 1}]
h = chSigmaULFormatter(h, rebin_pt=False, rebin_y=False, normalize=False)
h = hh.multiplyHists(h, h_sig_sigmaUL, flow=False)
writer.add_systematic(
    h,
    "pythiaew_ISRCorr0",
    "Zmumu",
    "chSigmaUL",
    mirror=True,
    kfactor=2.0,
    groups=[f"theory_ew_pythiaew_ISR", "theory_ew", "theory"],
)

# write output
directory = args.output
if directory == "":
    directory = "./"
filename = args.outname
if args.angularCoeffs:
    filename += f"_sigmaUL_Ais"
else:
    filename += f"_sigmaUL"
if args.postfix:
    filename += f"_{args.postfix}"
writer.write(outfolder=directory, outfilename=filename)
print(f"Written to {directory}/{filename}.hdf5")
