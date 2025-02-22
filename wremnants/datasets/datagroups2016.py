from utilities import boostHistHelpers as hh, logging
from wremnants import histselections as sel

logger = logging.child_logger(__name__)
    
def make_datagroups_2016(dg, combine=False, pseudodata_pdfset = None, applySelection=True, excludeGroups=None, filterGroups=None):
    # reset datagroups
    dg.groups = {}

    if dg.mode == "wmass" and applySelection:
        sigOp = sel.signalHistWmass
        fakeOp = sel.fakeHistABCD
    else:
        sigOp = None
        fakeOp = None

    dg.addGroup("Data",
        members = dg.getSafeListFromDataset(["dataPostVFP"]),
        selectOp = sigOp,
    )
    dg.addGroup("Zmumu",
        members = list(filter(lambda y: y.group == "Zmumu", dg.datasets.values())),
        selectOp = sigOp,
    ) 
    dg.addGroup("Ztautau",
        members = list(filter(lambda y: y.group == "Ztautau", dg.datasets.values())),
        selectOp = sigOp,
    )

    if pseudodata_pdfset and dg.combine:
        dg.addGroup(f"pdf{pseudodata_pdfset.upper()}_sum",
            label = f"pdf{pseudodata_pdfset.upper()}",
            color = "dimgray"
        )
    if dg.mode in ["vgen", "wmass"]:
        dg.addGroup("Wmunu",
            members = list(filter(lambda y: y.group == "Wmunu", dg.datasets.values())),
            selectOp = sigOp,
        )
        dg.addGroup("Wtaunu",
            members = list(filter(lambda y: y.group == "Wtaunu", dg.datasets.values())),
            selectOp = sigOp,
        )
        dg.addGroup("DYlowMass",
            members = list(filter(lambda y: y.group == "DYlowMass", dg.datasets.values())),
            selectOp = sigOp,
        )
        dg.addGroup("Top",
            members = list(filter(lambda y: y.group == "Top", dg.datasets.values())),
            selectOp = sigOp,
        )
        dg.addGroup("Diboson",
            members = list(filter(lambda y: y.group == "Diboson", dg.datasets.values())),
            selectOp = sigOp,
        )
        dg.addGroup("QCD",
            members = list(filter(lambda y: y.group == "QCD", dg.datasets.values())),
            selectOp = sigOp,
        )   
    else:
        dg.addGroup("Other",
            members = [x for x in dg.datasets.values() if not x.is_data and x.group not in ["Zmumu", "Ztautau"] and x.group != "QCD"],
        )

    dg.filterGroups(filterGroups)
    dg.excludeGroups(excludeGroups)

    if dg.mode == "wmass":
        # add all processes to the fake contributions after filtered and excluded groups
        dg.addGroup("Fake",
            members = [member for sublist in [v.members for k, v in dg.groups.items() if k != "QCD"] for member in sublist],
            scale = lambda x: 1. if x.is_data else -1,
            selectOp = fakeOp,
        )
        dg.filterGroups(filterGroups)
        dg.excludeGroups(excludeGroups)

    return dg
