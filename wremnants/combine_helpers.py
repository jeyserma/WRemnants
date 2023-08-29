from utilities import boostHistHelpers as hh, common, logging, input_tools
from wremnants import syst_tools,theory_tools,recoil_tools, unfolding_tools

from wremnants import histselections as sel
import hist
import numpy as np
import uproot
import h5py

logger = logging.child_logger(__name__)

def add_recoil_uncertainty(card_tool, samples, passSystToFakes=False, pu_type="highPU", flavor="", group_compact=True):
    met = input_tools.args_from_metadata(card_tool, "met")
    if flavor == "":
        flavor = input_tools.args_from_metadata(card_tool, "flavor")
    rtag = f"{pu_type}_{flavor}_{met}"
    if not rtag in recoil_tools.recoil_cfg:
        logger.warning(f"Recoil corrections for {pu_type}, {flavor}, {met} not available.")
        return
    recoil_cfg = recoil_tools.recoil_cfg[rtag]
    recoil_vars = list(recoil_cfg['corr_z'].keys()) + list(recoil_cfg['unc_z'].keys())
    recoil_grps = recoil_vars
    if group_compact:
        recoil_grps = ["CMS_recoil"]*len(recoil_cfg)
    for i, tag in enumerate(recoil_vars):
        card_tool.addSystematic("recoilUnc_%s" % tag,
            processes=samples,
            mirror = False,
            group = recoil_grps[i],
            systAxes = ["recoilVar"],
            passToFakes=passSystToFakes,
        )

def setSimultaneousABCD(cardTool, fakerate_integration_axes=[], fakerate_axes=["pt","eta","charge"],
    thresholdMT=40, integrateMT=False, axis_name_mt="mt", 
    variation_fakerate=0.5, variation_normalization_fake=0.1
):
    # Having 1 process for fakes, for each bin 3 free floating parameters, 2 normalization for lowMT and highMT and one fakerate between iso and anti iso
    logger.info(f"Set processes for simultaneous ABCD fit")

    # expected fake contribution
    hist_fake = sum([group.hists[cardTool.nominalName] if name == cardTool.getDataName() else -1*group.hists[cardTool.nominalName] for name, group in cardTool.datagroups.groups.items()])

    # setting errors to 0
    hist_fake.view(flow=True)[...] = np.stack((hist_fake.values(flow=True), np.zeros_like(hist_fake.values(flow=True))), axis=-1)

    failMT, passMT = sel.get_mt_selection(hist_fake, thresholdMT, axis_name_mt, integrateMT)
    passMTName = [k for k in passMT.keys()][0]

    if common.passIsoName not in hist_fake.axes.name or passMTName not in hist_fake.axes.name:
        raise RuntimeError(f'{common.passIsoName} and {passMTName} expected to be found in histogram, but only have axes {hist_fake.axes.name}')

    # axes in the correct ordering
    axes = cardTool.project[:]
    axes += [ax for ax in [common.passIsoName, passMTName] if ax not in axes]

    if set(hist_fake.axes.name) != set(axes) or hist_fake.axes.name[-2] != common.passIsoName or hist_fake.axes.name[-1] != passMTName:
        logger.info(f"Axes in histogram '{hist_fake.axes.name}' are not the same as required '{axes}' or in a different order than expected, try to project")
        hist_fake = hist_fake.project(*axes)

    # set the expected values in the signal region
    hist_fake.values(flow=True)[...,1,1] = sel.fakeHistABCD(hist_fake, fakerate_integration_axes=fakerate_integration_axes).values(flow=True)
    
    fakename = cardTool.getFakeName()

    cardTool.datagroups.addGroup(fakename, label = "Nonprompt", color = "grey", members=[],) #TODO check if existing group can be used
    cardTool.datagroups.groups[fakename].hists[f"{cardTool.nominalName}"] = hist_fake
    
    # axes in low MT
    fakerate_axes = [ax.name for ax in hist_fake.axes if ax.name not in [*fakerate_integration_axes, common.passIsoName, passMTName]]
    fakerate_bin_sizes = [ax.size for ax in hist_fake.axes if ax.name not in [*fakerate_integration_axes, common.passIsoName, passMTName]]
    # axes in high MT
    highMT_axes = [ax.name for ax in hist_fake.axes if ax.name in axes and ax.name not in [common.passIsoName, passMTName]]
    highMT_bin_sizes = [ax.size for ax in hist_fake.axes if ax.name in axes and ax.name not in [common.passIsoName, passMTName]]
    # axes in high MT that are not in lowMT
    other_axes = [n for n in highMT_axes if n not in fakerate_axes]
    other_bin_sizes = [s for n, s in zip(highMT_axes, highMT_bin_sizes) if n not in fakerate_axes]
    # all axes
    all_axes = fakerate_axes + other_axes

    if any(a in hist_fake.axes.name for a in fakerate_integration_axes):
        hist_failMT_failIso = hist_fake[{**common.failIso, **failMT}].project(fakerate_axes)
        hist_failMT_passIso = hist_fake[{**common.passIso, **failMT}].project(fakerate_axes)
    else:
        hist_failMT_failIso = hist_fake[{**common.failIso, **failMT}]
        hist_failMT_passIso = hist_fake[{**common.passIso, **failMT}]        
    
    hist_passMT = hist_fake[{**common.failIso, **passMT}] + hist_fake[{**common.passIso, **passMT}]

    # helper function to get indices from a 'flat' bin count
    def get_ax_idx(bin_sizes, i):
        current_i = i
        ax_idx = []
        for num in reversed(bin_sizes):
            ax_idx.insert(0, current_i % num)
            current_i //= num
        return ax_idx

    # loop over all fakerate axes bins
    for i in range(np.product(fakerate_bin_sizes)):

        fakerate_ax_idx = get_ax_idx(fakerate_bin_sizes, i)
        fakerate_bin_name = "_".join([f"{ax}{fakerate_ax_idx[j]}" for j, ax in enumerate(fakerate_axes)])
        fakerate_indices = {ax: fakerate_ax_idx[j] for j, ax in enumerate(fakerate_axes)}
        logger.debug(f"Now at fakerate bin {i}/{np.product(fakerate_bin_sizes)}: {fakerate_bin_name}")

        n_failMT_failIso = hist_failMT_failIso[fakerate_indices].value
        n_failMT = n_failMT_failIso + hist_failMT_passIso[fakerate_indices].value
        fr = n_failMT_failIso / n_failMT

        # systematic variation for fakerate, should be smaller 1 and bigger 0
        diff = min(fr, 1-fr)
        frUp = fr + variation_fakerate * diff
        hist_var_fakerate = hist.Hist(*hist_fake.axes, storage=hist.storage.Double())
        hist_var_fakerate.view(flow=True)[...] = hist_fake.values(flow=True)

        cardTool.datagroups.groups[fakename].hists[f"{cardTool.nominalName}_r{fakename}_{fakerate_bin_name}"] = hist_var_fakerate

        cardTool.addSystematic(f"r{fakename}_{fakerate_bin_name}",
            processes=[fakename],
            group=f"r{fakename}",
            noConstraint=True,
            outNames=[f"r{fakename}_{fakerate_bin_name}"],
            mirror=True
        )

        # systematic variation for fake normalization in low MT
        hist_var_lowMT = hist.Hist(*hist_fake.axes, storage=hist.storage.Double())
        hist_var_lowMT.view(flow=True)[...] = hist_fake.values(flow=True)

        cardTool.datagroups.groups[fakename].hists[f"{cardTool.nominalName}_N{fakename}LowMT_{fakerate_bin_name}"] = hist_var_lowMT

        cardTool.addSystematic(f"N{fakename}LowMT_{fakerate_bin_name}",
            processes=[fakename],
            group=f"{fakename}LowMT",
            noConstraint=True,
            outNames=[f"N{fakename}LowMT_{fakerate_bin_name}"],
            mirror=True
        )

        # loop over all other axes bins
        jMax = np.product(other_bin_sizes) if len(other_bin_sizes) else 1
        for j in range(jMax):
            other_ax_idx = get_ax_idx(other_bin_sizes, j)
            other_indices = {ax: other_ax_idx[k] for k, ax in enumerate(other_axes)}

            all_ax_idx = fakerate_ax_idx + other_ax_idx
            highMT_bin_name = "_".join([f"{ax}{all_ax_idx[k]}" for k, ax in enumerate(all_axes) if ax in highMT_axes])
            highMT_indices = {ax: all_ax_idx[k] for k, ax in enumerate(all_axes) if ax in highMT_axes}
            logger.debug(f"Now at other bin {j}/{np.product(highMT_indices)}: {highMT_bin_name}")

            # systematic variation for fake normalization in high MT, only define one time for each high mT bin
            if len(other_indices) == 0 or all(other_indices.values() == 0):
                hist_var_highMT = hist.Hist(*hist_fake.axes, storage=hist.storage.Double())
                hist_var_highMT.view(flow=True)[...] = hist_fake.values(flow=True)

                cardTool.datagroups.groups[fakename].hists[f"{cardTool.nominalName}_N{fakename}HighMT_{highMT_bin_name}"] = hist_var_highMT

                cardTool.addSystematic(f"N{fakename}HighMT_{highMT_bin_name}",
                    processes=[fakename],
                    group=f"{fakename}HighMT",
                    noConstraint=True,
                    outNames=[f"N{fakename}HighMT_{highMT_bin_name}"],
                    mirror=True
                )
            else:
                hist_var_highMT = cardTool.datagroups.groups[fakename].hists[f"{cardTool.nominalName}_N{fakename}HighMT_{highMT_bin_name}"]

            n_passMT = hist_passMT[{**fakerate_indices, **other_indices}].value

            # set fakerate variation histogram
            hist_var_fakerate[{**common.failIso, **failMT, **fakerate_indices, **other_indices}] = n_failMT * frUp
            hist_var_fakerate[{**common.passIso, **failMT, **fakerate_indices, **other_indices}] = n_failMT * (1-frUp)
            hist_var_fakerate[{**common.failIso, **passMT, **fakerate_indices, **other_indices}] = n_passMT * frUp
            hist_var_fakerate[{**common.passIso, **passMT, **fakerate_indices, **other_indices}] = n_passMT * (1-frUp)

            # set lowMT variations histogram
            hist_var_lowMT[{**common.failIso, **failMT, **fakerate_indices, **other_indices}] = (1+variation_normalization_fake) * n_failMT * fr
            hist_var_lowMT[{**common.passIso, **failMT, **fakerate_indices, **other_indices}] = (1+variation_normalization_fake) * n_failMT * (1-fr)

            # set highMT variations histogram
            hist_var_highMT[{**common.failIso, **passMT, **fakerate_indices, **other_indices}] = (1+variation_normalization_fake) * n_passMT * fr
            hist_var_highMT[{**common.passIso, **passMT, **fakerate_indices, **other_indices}] = (1+variation_normalization_fake) * n_passMT * (1-fr)


def getTheoryFitData(fitresult, axes=None, base_processes = "W", poi_type="pmaskedexp"):
    logger.info(f"Prepare theory fit: load measured differential cross secction distribution and covariance matrix")

    if fitresult.endswith(".root"):
        if project is None:
            raise RuntimeError("When fitresult is provided as root file the axes need to be specified")
        
        rfile = uproot.open(fitresult)
        df = unfolding_tools.get_results(rfile, poi_type)

        # write out unfolded data as 1D hist
        hist_xsec = hist.Hist(
            hist.axis.Regular(bins=len(df), start=0.5, stop=len(df)+0.5, underflow=False, overflow=False), storage=hist.storage.Weight())
        hist_xsec.view(flow=False)[...] = np.stack([df["value"].values, (df["err_total"].values)**2], axis=-1)

        data = hist_xsec.values(flow=False).flatten()

        # write out covariance as 2D hist
        cov = unfolding_tools.matrix_poi(rfile, poi_type, base_process=base_processes[0], axes=axes).values(flow=False)
    elif fitresult.endswith(".hdf5"):
        hfile = h5py.File(fitresult, mode='r')

        outvals = hfile[f"{poi_type}_outvals"][...]
        npoi = len(outvals)
        # make matrix between POIs only; assume POIs come first (which should be the case in combinetf)
        outcov = hfile[f"{poi_type}_outcov"][:npoi,:npoi]

        # select POIs for each base process, assume correct order
        all_indices = np.zeros(npoi, dtype=bool)
        data = []
        for p in base_processes:
            indices = np.array([s.startswith(p) for s in hfile[f"{poi_type}_names"][...].astype(str)], dtype=bool)

            data.append(outvals[indices])

            all_indices = all_indices | indices
        #select rows and columns
        cov = outcov[all_indices][:, all_indices]

    else:
        raise NotImplementedError(f"Unkown data type for fitresult {fitresult}")
    
    return data, cov
