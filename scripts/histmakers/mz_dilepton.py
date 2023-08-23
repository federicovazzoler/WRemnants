from utilities import boostHistHelpers as hh, common, output_tools, logging, differential

parser,initargs = common.common_parser(True)

import ROOT
import narf
import wremnants
from wremnants import theory_tools,syst_tools,theory_corrections, muon_validation, muon_calibration, muon_selections, unfolding_tools, helicity_utils
from wremnants.histmaker_tools import scale_to_data, aggregate_groups
from wremnants.datasets.dataset_tools import getDatasets
import hist
import lz4.frame
import math
import time
import os


parser.add_argument("--csVarsHist", action='store_true', help="Add CS variables to dilepton hist")
parser.add_argument("--axes", type=str, nargs="*", default=["mll", "ptll"], help="")
parser.add_argument("--finePtBinning", action='store_true', help="Use fine binning for ptll")
parser.add_argument("--noAuxiliaryHistograms", action="store_true", help="Remove auxiliary histograms to save memory (removed by default with --unfolding or --theoryAgnostic)")
parser.add_argument("--addHelicityHistos", action='store_true', help="Add V qT,Y axes and helicity axes for W samples")
parser.add_argument("--onlyTheorySyst", action="store_true", help="Keep only theory systematic variations, mainly for tests")

parser = common.set_parser_default(parser, "genVars", ["ptVGen", "absYVGen"])
parser = common.set_parser_default(parser, "pt", [34,26.,60.])
parser = common.set_parser_default(parser, "eta", [48,-2.4,2.4])
parser = common.set_parser_default(parser, "aggregateGroups", ["Diboson", "Top", "Wtaunu", "Wmunu"])

args = parser.parse_args()
logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

thisAnalysis = ROOT.wrem.AnalysisType.Dilepton

datasets = getDatasets(maxFiles=args.maxFiles,
                        filt=args.filterProcs,
                        excl=args.excludeProcs, 
                        nanoVersion="v8" if args.v8 else "v9", base_path=args.dataPath)

era = args.era

# dilepton invariant mass cuts
mass_min = 60
mass_max = 120

ewMassBins = theory_tools.make_ew_binning(mass = 91.1535, width = 2.4932, initialStep=0.010)

# available axes for dilepton validation plots
all_axes = {
    "mll": hist.axis.Regular(60, 60., 120., name = "mll", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "yll": hist.axis.Regular(20, -2.5, 2.5, name = "yll", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "absYll": hist.axis.Regular(10, 0., 2.5, name = "absYll", underflow=False, overflow=not args.excludeFlow),
    "ptll": hist.axis.Variable(common.ptV_binning if not args.finePtBinning else range(60), name = "ptll", underflow=False, overflow=not args.excludeFlow),
    "etaPlus": hist.axis.Regular(int(args.eta[0]), args.eta[1], args.eta[2], name = "etaPlus"),
    "etaMinus": hist.axis.Regular(int(args.eta[0]), args.eta[1], args.eta[2], name = "etaMinus"),
    "etaSum": hist.axis.Regular(12, -4.8, 4.8, name = "etaSum"),
    "etaDiff": hist.axis.Variable([-4.8, -1.0, -0.6, -0.2, 0.2, 0.6, 1.0, 4.8], name = "etaDiff"),
    "ptPlus": hist.axis.Regular(int(args.pt[0]), args.pt[1], args.pt[2], name = "ptPlus"),
    "ptMinus": hist.axis.Regular(int(args.pt[0]), args.pt[1], args.pt[2], name = "ptMinus"),
    "cosThetaStarll": hist.axis.Regular(20, -1., 1., name = "cosThetaStarll", underflow=False, overflow=False),
    "phiStarll": hist.axis.Regular(20, -math.pi, math.pi, circular = True, name = "phiStarll"),
    #"charge": hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge") # categorical axes in python bindings always have an overflow bin, so use a regular
    "massVgen": hist.axis.Variable(ewMassBins, name = "massVgen", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "ewMll": hist.axis.Variable(ewMassBins, name = "ewMll", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "ewMlly": hist.axis.Variable(ewMassBins, name = "ewMlly", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
    "ewLogDeltaM": hist.axis.Regular(100, -10, 4, name = "ewLogDeltaM", overflow=not args.excludeFlow, underflow=not args.excludeFlow),
}

auxiliary_gen_axes = ["massVgen", # preFSR variables
    "ewMll", "ewMlly", "ewLogDeltaM" # ew variables
    ]

for a in args.axes:
    if a not in all_axes.keys():
        logger.error(f" {a} is not a known axes! Supported axes choices are {list(axes.keys())}")

nominal_cols = args.axes

if args.csVarsHist:
    nominal_cols += ["cosThetaStarll", "phiStarll"]

nominal_axes = [all_axes[a] for a in nominal_cols] 

groups_to_aggregate = args.aggregateGroups

gen_axes = {
    "ptVGen": hist.axis.Variable(common.ptV_binning if not args.finePtBinning else range(60), name = "ptVGen", underflow=False, overflow=False),
    "absYVGen": hist.axis.Regular(10, 0, 2.5, name = "absYVGen", underflow=False, overflow=False),  
}

if args.unfolding:
    unfolding_axes, unfolding_cols, unfolding_selections = differential.get_dilepton_axes(args.genVars, gen_axes)
    datasets = unfolding_tools.add_out_of_acceptance(datasets, group = "Zmumu")

elif args.addHelicityHistos:

    #list(range(0,50,5)).append(np.inf) ,
    axis_ptVgen = hist.axis.Variable(
        [0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50.],
        #[0., 2.5, 5., 7.5, 10., 12.5, 15., 17.5, 20., 22.5, 25., 30., 35., 40., 45., 50.],
        name = "ptVgenSig", underflow=False, overflow=False
    )
    #axis_ptVgen.append(np.inf)
    #Taken from w-z gen histomaker     
    axis_absYVgen = hist.axis.Variable(
        [0, 0.5, 1., 1.5, 2.0, 2.5],
        #[0.25*i for i in range(11)],
        #[0, 1.25, 2.5],
        name = "absYVgenSig", underflow=False, overflow=False
    )
    theoryAgnostic_axes = [axis_absYVgen, axis_ptVgen]
    theoryAgnostic_cols = ["absYVgen", "ptVgen"] # name of the branch, not of the axis
    axis_helicity = helicity_utils.axis_helicity_multidim
    # axis_helicity = hist.axis.Integer(-1, 5, name="helicity", overflow=False, underflow=False)
    # the following just prepares the existence of the group for out-of-acceptance signal, but doesn't create or define the histogram yet
    datasets = unfolding_tools.add_out_of_acceptance(datasets, group = "Zmumu")
    groups_to_aggregate.append("BkgZmumu")

# define helpers
muon_prefiring_helper, muon_prefiring_helper_stat, muon_prefiring_helper_syst = wremnants.make_muon_prefiring_helpers(era = era)

qcdScaleByHelicity_helper = wremnants.makeQCDScaleByHelicityHelper(is_w_like = True)

# extra axes which can be used to label tensor_axes
if args.binnedScaleFactors:
    logger.info("Using binned scale factors and uncertainties")
    # might never use it really anymore, but let's warn the user that this is obsolete
    logger.warning("Only SF with no uT dependence are implemented, and the treatment for trigger is like Wlike")
    # add usePseudoSmoothing=True for tests with Asimov
    muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_binned(filename = args.sfFile,
                                                                                                                                     era = era,
                                                                                                                                     max_pt = args.pt[2],
                                                                                                                                     is_w_like = True) 
else:
    logger.info("Using smoothed scale factors and uncertainties")
    muon_efficiency_helper, muon_efficiency_helper_syst, muon_efficiency_helper_stat = wremnants.make_muon_efficiency_helpers_smooth(filename = args.sfFile,
                                                                                                                                     era = era,
                                                                                                                                     max_pt = args.pt[2],
                                                                                                                                     what_analysis = thisAnalysis, isoEfficiencySmoothing=args.isoEfficiencySmoothing, smooth3D=args.smooth3dsf)
logger.info(f"SF file: {args.sfFile}")

pileup_helper = wremnants.make_pileup_helper(era = era)

calib_filepaths = common.calib_filepaths
closure_filepaths = common.closure_filepaths
diff_weights_helper = ROOT.wrem.SplinesDifferentialWeightsHelper(calib_filepaths['tflite_file']) if (args.muonScaleVariation == 'smearingWeightsSplines' or args.validationHists) else None
mc_jpsi_crctn_helper, data_jpsi_crctn_helper, mc_jpsi_crctn_unc_helper, data_jpsi_crctn_unc_helper = muon_calibration.make_jpsi_crctn_helpers(args, calib_filepaths, make_uncertainty_helper=True)
z_non_closure_parametrized_helper, z_non_closure_binned_helper = muon_calibration.make_Z_non_closure_helpers(args, calib_filepaths, closure_filepaths)

mc_calibration_helper, data_calibration_helper, calibration_uncertainty_helper = muon_calibration.make_muon_calibration_helpers(args)

smearing_helper, smearing_uncertainty_helper = (None, None) if args.noSmearing else muon_calibration.make_muon_smearing_helpers()

bias_helper = muon_calibration.make_muon_bias_helpers(args) 

corr_helpers = theory_corrections.load_corr_helpers([d.name for d in datasets if d.name in common.vprocs], args.theoryCorr)

######################################################
######################################################
######################################################
## FIXME:
## move next functions below to some other place to import from
def select_fiducial_space_theoryAgnostic(df, ptVgenMax, absYVgenMax, accept=True):
    # might need to account for some Ai which we neglect, we may have an out-oc-acceptance per helicity state
    selection = f"ptVgen < {ptVgenMax} && absYVgen < {absYVgenMax}"
    df = df.Define("fiducial", selection)
    if accept:
        df = df.Filter("fiducial")
        logger.debug(f"Theory agnostic fiducial cut: {selection}")
    else:
        df = df.Filter("fiducial == 0")
        logger.debug(f"Theory agnostic fiducial cut (out-of-acceptance): not ({selection})")
    return df

def define_helicity_weights(df):
    # define the helicity tensor, here nominal_weight will only have theory weights, no experimental pieces, it is defined in theory_tools.define_theory_weights_and_corrs
    weightsByHelicity_helper = wremnants.makehelicityWeightHelper()
    df = df.Define("helWeight_tensor", weightsByHelicity_helper, ["massVgen", "yVgen", "ptVgen", "chargeVgen", "csSineCosThetaPhi", "nominal_weight"])
    df = df.Define("nominal_weight_helicity", "wrem::scalarmultiplyHelWeightTensor(nominal_weight,helWeight_tensor)")
    return df

def add_xnorm_histograms_theoryAgnostic(results, df, args, dataset_name, corr_helpers, qcdScaleByHelicity_helper, theoryAgnostic_axes, theoryAgnostic_cols):
    # add histograms before any selection
    df_xnorm = df
    df_xnorm = df_xnorm.DefinePerSample("exp_weight", "1.0")
    df_xnorm = theory_tools.define_theory_weights_and_corrs(df_xnorm, dataset_name, corr_helpers, args)
    # define the helicity tensor, here nominal_weight will only have theory weights, no experimental pieces, it is defined in theory_tools.define_theory_weights_and_corrs
    df_xnorm = define_helicity_weights(df_xnorm)
    df_xnorm = df_xnorm.DefinePerSample("xnorm", "0.5")
    axis_xnorm = hist.axis.Regular(1, 0., 1., name = "count", underflow=False, overflow=False)
    xnorm_axes = [axis_xnorm, *theoryAgnostic_axes]
    xnorm_cols = ["xnorm", *theoryAgnostic_cols]
    xnormByHelicity = df_xnorm.HistoBoost("xnorm", xnorm_axes, [*xnorm_cols, "nominal_weight_helicity"], tensor_axes=[axis_helicity])
    results.append(xnormByHelicity)
    if not args.onlyMainHistograms:
        syst_tools.add_theory_hists(results, df_xnorm, args, dataset_name, corr_helpers, qcdScaleByHelicity_helper, xnorm_axes, xnorm_cols, base_name="xnorm", for_wmass=True, addhelicity=True)
    else:
        #FIXME: hardcoded to keep mass weights, this would be done in add_theory_hists
        df_xnorm = syst_tools.define_mass_weights(df_xnorm, dataset_name)
        syst_tools.add_massweights_hist(results, df_xnorm, xnorm_axes, xnorm_cols, base_name="xnorm", proc=dataset_name, addhelicity=True)
        
# graph building for W sample with helicity weights
def whistosbyHelicity(df, results, dataset, nominal_axes_thAgn, nominal_cols_thAgn, args):
    logger.info("Entered function whistobyHelicity")
    df = define_helicity_weights(df)
    nominalByHelicity = df.HistoBoost("nominal", nominal_axes_thAgn, [*nominal_cols_thAgn, "nominal_weight_helicity"], tensor_axes=[axis_helicity])
    results.append(nominalByHelicity)

    if not args.onlyMainHistograms:
        if not args.onlyTheorySyst:
            df = syst_tools.add_L1Prefire_unc_hists(results, df, muon_prefiring_helper_stat, muon_prefiring_helper_syst, nominal_axes_thAgn, nominal_cols_thAgn, addhelicity=True)
            df = syst_tools.add_muon_efficiency_unc_hists(results, df, muon_efficiency_helper_stat, muon_efficiency_helper_syst, nominal_axes_thAgn, nominal_cols_thAgn, what_analysis=thisAnalysis, addhelicity=True)
        df = syst_tools.add_theory_hists(results, df, args, dataset.name, corr_helpers, qcdScaleByHelicity_helper, nominal_axes_thAgn, nominal_cols_thAgn, for_wmass=True, addhelicity=True)
    else:
        #FIXME: hardcoded to keep mass weights, this would be done in add_theory_hists
        df = syst_tools.define_mass_weights(df, dataset.name)
        syst_tools.add_massweights_hist(results, df, nominal_axes_thAgn, nominal_cols_thAgn, proc=dataset.name, addhelicity=True)
######################################################
######################################################
######################################################

smearing_weights_procs = []

def build_graph(df, dataset):
    logger.info(f"build graph for dataset: {dataset.name}")
    results = []
    isW = dataset.name in common.wprocs
    isZ = dataset.name in common.zprocs
    isWorZ = isW or isZ

    cvh_helper = data_calibration_helper if dataset.is_data else mc_calibration_helper
    jpsi_helper = data_jpsi_crctn_helper if dataset.is_data else mc_jpsi_crctn_helper
    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        df = df.Define("weight", "std::copysign(1.0, genWeight)")

    weightsum = df.SumAndCount("weight")

    axes = nominal_axes
    cols = nominal_cols

    if args.unfolding and dataset.name == "ZmumuPostVFP":
        df = unfolding_tools.define_gen_level(df, args.genLevel, dataset.name, mode="dilepton")

        if hasattr(dataset, "out_of_acceptance"):
            logger.debug("Reject events in fiducial phase space")
            df = unfolding_tools.select_fiducial_space(df, mode="dilepton", pt_min=args.pt[1], pt_max=args.pt[2], 
                mass_min=mass_min, mass_max=mass_max, selections=unfolding_selections, accept=False)
        else:
            logger.debug("Select events in fiducial phase space")
            df = unfolding_tools.select_fiducial_space(df, mode="dilepton", pt_min=args.pt[1], pt_max=args.pt[2], 
                mass_min=mass_min, mass_max=mass_max, selections=unfolding_selections, accept=True)

            unfolding_tools.add_xnorm_histograms(results, df, args, dataset.name, corr_helpers, qcdScaleByHelicity_helper, unfolding_axes, unfolding_cols)
            axes = [*nominal_axes, *unfolding_axes] 
            cols = [*nominal_cols, *unfolding_cols]

    if not args.noAuxiliaryHistograms and isZ:
        # gen level variables before selection
        for obs in auxiliary_gen_axes:
            df_gen = df
            df_gen = df_gen.DefinePerSample("exp_weight", "1.0")

            df_gen = theory_tools.define_theory_weights_and_corrs(df_gen, dataset.name, corr_helpers, args)

            results.append(df_gen.HistoBoost(f"gen_{obs}", [all_axes[obs]], [obs, "nominal_weight"]))
            df_gen = syst_tools.add_theory_hists(results, df_gen, args, dataset.name, corr_helpers, qcdScaleByHelicity_helper, [all_axes[obs]], [obs], base_name=f"gen_{obs}", for_wmass=False)
    if args.addHelicityHistos and isWorZ: # should be isW to do also Wtaunu
        df = theory_tools.define_prefsr_vars(df)
        if hasattr(dataset, "out_of_acceptance"):
            logger.debug("Reject events in fiducial phase space")
            df = select_fiducial_space_theoryAgnostic(df, axis_ptVgen.edges[-1], axis_absYVgen.edges[-1], accept=False)
        else:
            logger.debug("Select events in fiducial phase space for theory agnostic analysis")
            df = select_fiducial_space_theoryAgnostic(df, axis_ptVgen.edges[-1], axis_absYVgen.edges[-1], accept=True)
            add_xnorm_histograms_theoryAgnostic(results, df, args, dataset.name, corr_helpers, qcdScaleByHelicity_helper, theoryAgnostic_axes, theoryAgnostic_cols)
            # helicity axis is special, defined through a tensor later, theoryAgnostic_ only includes W rapidity and pt for now
            axes = [*nominal_axes, *theoryAgnostic_axes]
            cols = [*nominal_cols, *theoryAgnostic_cols]

    df = df.Filter("HLT_IsoTkMu24 || HLT_IsoMu24")

    df = muon_selections.veto_electrons(df)
    df = muon_selections.apply_met_filters(df)

    df = muon_calibration.define_corrected_muons(df, cvh_helper, jpsi_helper, args, dataset, smearing_helper, bias_helper)

    df = muon_selections.select_veto_muons(df, nMuons=2)
    df = muon_selections.select_good_muons(df, args.pt[1], args.pt[2], nMuons=2, use_trackerMuons=args.trackerMuons, use_isolation=True)

    # for dilepton analysis we will call trigMuons (nonTrigMuons) those with charge plus (minus). In fact both might be triggering, naming scheme might be improved
    df = muon_selections.define_trigger_muons(df, what_analysis=thisAnalysis)

    df = muon_selections.select_z_candidate(df, mass_min, mass_max)

    df = muon_selections.select_standalone_muons(df, dataset, args.trackerMuons, "trigMuons")
    df = muon_selections.select_standalone_muons(df, dataset, args.trackerMuons, "nonTrigMuons")

    df = muon_selections.apply_triggermatching_muon(df, dataset, "trigMuons_eta0", "trigMuons_phi0", "nonTrigMuons_eta0", "nonTrigMuons_phi0")

    df = df.Define("ptll", "ll_mom4.pt()")
    df = df.Define("yll", "ll_mom4.Rapidity()")
    df = df.Define("absYll", "std::fabs(yll)")
    # "renaming" to write out corresponding axis
    df = df.Define("etaPlus", "trigMuons_eta0")
    df = df.Define("etaMinus", "nonTrigMuons_eta0")
    df = df.Define("ptPlus", "trigMuons_pt0")
    df = df.Define("ptMinus", "nonTrigMuons_pt0")

    df = df.Define("etaSum", "nonTrigMuons_eta0 + trigMuons_eta0") 
    df = df.Define("etaDiff", "trigMuons_eta0-nonTrigMuons_eta0") # plus - minus 

    if args.csVarsHist:
        df = df.Define("csSineCosThetaPhill", "wrem::csSineCosThetaPhi(nonTrigMuons_mom4, trigMuons_mom4)")
        df = df.Define("cosThetaStarll", "csSineCosThetaPhill.costheta")
        df = df.Define("phiStarll", "std::atan2(csSineCosThetaPhill.sinphi, csSineCosThetaPhill.cosphi)")

    logger.debug(f"Define weights and store nominal histograms")

    if dataset.is_data:
        results.append(df.HistoBoost("nominal", axes, cols))
    else:
        df = df.Define("weight_pu", pileup_helper, ["Pileup_nTrueInt"])
        columnsForSF = ["trigMuons_pt0", "trigMuons_eta0", "trigMuons_SApt0", "trigMuons_SAeta0", "trigMuons_uT0", "trigMuons_charge0", "trigMuons_passTrigger0",
                        "nonTrigMuons_pt0", "nonTrigMuons_eta0", "nonTrigMuons_SApt0", "nonTrigMuons_SAeta0", "nonTrigMuons_uT0", "nonTrigMuons_charge0", "nonTrigMuons_passTrigger0"]
        df = muon_selections.define_muon_uT_variable(df, isWorZ, smooth3dsf=args.smooth3dsf, colNamePrefix="trigMuons")
        df = muon_selections.define_muon_uT_variable(df, isWorZ, smooth3dsf=args.smooth3dsf, colNamePrefix="nonTrigMuons")
        if not args.smooth3dsf:
            columnsForSF.remove("trigMuons_uT0")
            columnsForSF.remove("nonTrigMuons_uT0")

        # FIXME: add flags for pass_trigger for both leptons
        df = df.Define("weight_fullMuonSF_withTrackingReco", muon_efficiency_helper, columnsForSF)
        df = df.Define("weight_newMuonPrefiringSF", muon_prefiring_helper, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId"])

        df = df.Define("exp_weight", "weight_pu*weight_fullMuonSF_withTrackingReco*weight_newMuonPrefiringSF*L1PreFiringWeight_ECAL_Nom")
        df = theory_tools.define_theory_weights_and_corrs(df, dataset.name, corr_helpers, args)

        results.append(df.HistoBoost("weight", [hist.axis.Regular(100, -2, 2)], ["nominal_weight"], storage=hist.storage.Double()))
        results.append(df.HistoBoost("nominal", axes, [*cols, "nominal_weight"]))

    for obs in ["ptll", "mll", "yll", "etaPlus", "etaMinus", "ptPlus", "ptMinus"]:
        if dataset.is_data:
            results.append(df.HistoBoost(f"nominal_{obs}", [all_axes[obs]], [obs]))
        else:
            results.append(df.HistoBoost(f"nominal_{obs}", [all_axes[obs]], [obs, "nominal_weight"]))
            if isWorZ:
                df = syst_tools.add_theory_hists(results, df, args, dataset.name, corr_helpers, qcdScaleByHelicity_helper, [all_axes[obs]], [obs], base_name=f"nominal_{obs}", for_wmass=False)

    if not args.noAuxiliaryHistograms and isZ:
        # gen level variables
        for obs in auxiliary_gen_axes:
            results.append(df.HistoBoost(f"nominal_{obs}", [all_axes[obs]], [obs, "nominal_weight"]))
            df = syst_tools.add_theory_hists(results, df, args, dataset.name, corr_helpers, qcdScaleByHelicity_helper, [all_axes[obs]], [obs], base_name=f"nominal_{obs}", for_wmass=False)

    # test plots
    if args.validationHists:
        df_plusTrig = df.Filter("trigMuons_passTrigger0")
        df_minusTrig = df.Filter("nonTrigMuons_passTrigger0")
        df_bothTrig = df.Filter("trigMuons_passTrigger0 && nonTrigMuons_passTrigger0")
        df_plusTrigOnly = df.Filter("trigMuons_passTrigger0 && !nonTrigMuons_passTrigger0")
        df_minusTrigOnly = df.Filter("nonTrigMuons_passTrigger0 && !trigMuons_passTrigger0")
        for obs in ["etaPlus", "etaMinus", "ptPlus", "ptMinus"]:
            if dataset.is_data:
                results.append(df_plusTrig.HistoBoost(f"nominal_{obs}_plusTrig", [all_axes[obs]], [obs]))
                results.append(df_minusTrig.HistoBoost(f"nominal_{obs}_minusTrig", [all_axes[obs]], [obs]))
                results.append(df_bothTrig.HistoBoost(f"nominal_{obs}_bothTrig", [all_axes[obs]], [obs]))
                results.append(df_plusTrigOnly.HistoBoost(f"nominal_{obs}_plusTrigOnly", [all_axes[obs]], [obs]))
                results.append(df_minusTrigOnly.HistoBoost(f"nominal_{obs}_minusTrigOnly", [all_axes[obs]], [obs]))
            else:
                results.append(df_plusTrig.HistoBoost(f"nominal_{obs}_plusTrig", [all_axes[obs]], [obs, "nominal_weight"]))
                results.append(df_minusTrig.HistoBoost(f"nominal_{obs}_minusTrig", [all_axes[obs]], [obs, "nominal_weight"]))
                results.append(df_bothTrig.HistoBoost(f"nominal_{obs}_bothTrig", [all_axes[obs]], [obs, "nominal_weight"]))
                results.append(df_plusTrigOnly.HistoBoost(f"nominal_{obs}_plusTrigOnly", [all_axes[obs]], [obs, "nominal_weight"]))
                results.append(df_minusTrigOnly.HistoBoost(f"nominal_{obs}_minusTrigOnly", [all_axes[obs]], [obs, "nominal_weight"]))

    if not dataset.is_data and not args.onlyMainHistograms:

        df = syst_tools.add_muon_efficiency_unc_hists(results, df, muon_efficiency_helper_stat, muon_efficiency_helper_syst, axes, cols, what_analysis=thisAnalysis, smooth3D=args.smooth3dsf)
        df = syst_tools.add_L1Prefire_unc_hists(results, df, muon_prefiring_helper_stat, muon_prefiring_helper_syst, axes, cols)

        # n.b. this is the W analysis so mass weights shouldn't be propagated
        # on the Z samples (but can still use it for dummy muon scale)
        if isWorZ:

            df = syst_tools.add_theory_hists(results, df, args, dataset.name, corr_helpers, qcdScaleByHelicity_helper, axes, cols, for_wmass=False)

            reco_sel = "vetoMuonsPre"
            require_prompt = "tau" not in dataset.name
            df = muon_calibration.define_genFiltered_recoMuonSel(df, reco_sel, require_prompt)
            reco_sel_GF = muon_calibration.getColName_genFiltered_recoMuonSel(reco_sel, require_prompt)
            df = muon_calibration.define_matched_gen_muons_kinematics(df, reco_sel_GF)
            df = muon_calibration.calculate_matched_gen_muon_kinematics(df, reco_sel_GF)
            df = muon_calibration.define_matched_reco_muon_kinematics(df, reco_sel_GF)

            ####################################################
            # nuisances from the muon momemtum scale calibration 
            if (args.muonCorrData in ["massfit", "lbl_massfit"]):
                input_kinematics = [
                    f"{reco_sel_GF}_recoPt",
                    f"{reco_sel_GF}_recoEta",
                    f"{reco_sel_GF}_recoCharge",
                    f"{reco_sel_GF}_genPt",
                    f"{reco_sel_GF}_genEta",
                    f"{reco_sel_GF}_genCharge"
                ]
                if diff_weights_helper:
                    df = df.Define(f'{reco_sel_GF}_response_weight', diff_weights_helper, [*input_kinematics])
                    input_kinematics.append(f'{reco_sel_GF}_response_weight')

                # muon scale variation from stats. uncertainty on the jpsi massfit
                df = df.Define(
                    "nominal_muonScaleSyst_responseWeights_tensor", data_jpsi_crctn_unc_helper,
                    [*input_kinematics, "nominal_weight"]
                )
                muonScaleSyst_responseWeights = df.HistoBoost(
                    "nominal_muonScaleSyst_responseWeights", axes,
                    [*cols, "nominal_muonScaleSyst_responseWeights_tensor"],
                    tensor_axes = data_jpsi_crctn_unc_helper.tensor_axes, storage=hist.storage.Double()
                )
                results.append(muonScaleSyst_responseWeights)

                df = muon_calibration.add_resolution_uncertainty(df, axes, results, cols, smearing_uncertainty_helper, reco_sel_GF)

                # add the ad-hoc Z non-closure nuisances from the jpsi massfit to muon scale unc
                df = df.DefinePerSample("AFlag", "0x01")
                df = df.Define(
                    "Z_non_closure_parametrized_A", z_non_closure_parametrized_helper,
                    [*input_kinematics, "nominal_weight", "AFlag"]
                )
                hist_Z_non_closure_parametrized_A = df.HistoBoost(
                    "nominal_Z_non_closure_parametrized_A",
                    axes, [*cols, "Z_non_closure_parametrized_A"],
                    tensor_axes = z_non_closure_parametrized_helper.tensor_axes,
                    storage=hist.storage.Double()
                )
                results.append(hist_Z_non_closure_parametrized_A)

                df = df.DefinePerSample("MFlag", "0x04")
                df = df.Define(
                    "Z_non_closure_parametrized_M", z_non_closure_parametrized_helper,
                    [*input_kinematics, "nominal_weight", "MFlag"]
                )
                hist_Z_non_closure_parametrized_M = df.HistoBoost(
                    "nominal_Z_non_closure_parametrized_M",
                    axes, [*cols, "Z_non_closure_parametrized_M"],
                    tensor_axes = z_non_closure_parametrized_helper.tensor_axes,
                    storage=hist.storage.Double()
                )
                results.append(hist_Z_non_closure_parametrized_M)
            ####################################################

            # Don't think it makes sense to apply the mass weights to scale leptons from tau decays
            if not "tau" in dataset.name:
                syst_tools.add_muonscale_hist(
                    results, df, args.muonCorrEtaBins, args.muonCorrMag, isW, axes, cols,
                    muon_eta="trigMuons_eta0") ## FIXME: what muon to choose ?

    if hasattr(dataset, "out_of_acceptance"):
        # Rename dataset to not overwrite the original one
        dataset.name = "Bkg"+dataset.name

    if isWorZ and args.addHelicityHistos:
        whistosbyHelicity(df, results, dataset, axes, cols, args)
        ## TODO: this part should be better melted in the rest of the code, there is too much duplication of what could happen later in the loop
        if hasattr(dataset, "out_of_acceptance"):
            # Rename dataset to not overwrite the original one
            dataset.name = "Bkg"+dataset.name
    
    
    return results, weightsum

resultdict = narf.build_and_run(datasets, build_graph)

if not args.noScaleToData:
    scale_to_data(resultdict)
    aggregate_groups(datasets, resultdict, args.aggregateGroups)

output_tools.write_analysis_output(resultdict, f"{os.path.basename(__file__).replace('py', 'hdf5')}", args, update_name=not args.forceDefaultName)
