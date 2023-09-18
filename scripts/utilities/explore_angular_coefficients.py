#!/usr/bin/env python3

import os
import sys
from utilities import boostHistHelpers as hh, logging
from utilities import input_tools
from wremnants import theory_tools
import h5py
from narf import ioutils
import ROOT
from wremnants import plot_tools
import matplotlib.pyplot as plt
import mplhep as hep
import pdb

def load_helicity_moments_for_sample_from_file(sampleName = "ZmumuPostVFP", filePath = None):
    if not os.path.exists(filePath):
        raise Exception(f"{filePath} doesn't exits")

    file = h5py.File(filePath, "r")
    res = narf.ioutils.pickle_load_h5py(file["results"])
    h = input_tools.load_and_scale(res, sampleName, "helicity_moments_scale")[{"muRfact" : 1.j, "muFfact" : 1.j}].project('massVgen', 'y', 'ptVgen', 'chargeVgen', 'helicity')

    return h

def make_Ais_for_observable(obs, h, rebin = None):
    if rebin is not None:
        h = hh.rebinHist(h, obs, rebin)
    h_coeff = theory_tools.moments_to_angular_coeffs(h.project(obs, "helicity"))

    return h_coeff

def plot_hist(h, xlabel = "", ylabel = "", xrange = None, yrange = None, corr = 1):
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (8, 8), sharex = True)

    xValues = h.axes.centers[0]
    yValues = h.values() * corr
#    pdb.set_trace()
    plt.plot(xValues, yValues, marker = '.', linestyle = 'solid')

    ax.set_xlabel(xlabel, fontsize = 22)
    ax.set_ylabel(ylabel, fontsize = 22)
    if xrange is not None:
        ax.set_xlim(xrange)
    if yrange is not None:
        ax.set_ylim(yrange)

    return fig

def save_fig(fig, outname):
    fig.set_tight_layout(True)
    fig.savefig(outname)
    print(f"{outname} created")

def plot_Ai_from_hist_for_observable(Ai, h, obs, outname):
    hToPlot = h[{"helicity" : complex(Ai)}]
    corr = 1
    if Ai == 4: corr = -1
    save_fig(plot_hist(hToPlot, xlabel = obs, ylabel = f"$A_{Ai}$", xrange = None, yrange = None, corr = corr), outname)

def main():
    filePath = "/eos/home-f/fvazzole/TheoryAgnostic/AngularCoefficients/w_z_gen_dists_NonClosureCorl.hdf5" #FIXME

    # Load the helicity moments
    minnloZHelMom = load_helicity_moments_for_sample_from_file(sampleName = "ZmumuPostVFP", filePath = filePath)
    minnloWmHelMom = load_helicity_moments_for_sample_from_file(sampleName = "WminusmunuPostVFP", filePath = filePath)
    minnloWpHelMom = load_helicity_moments_for_sample_from_file(sampleName = "WplusmunuPostVFP", filePath = filePath)

    # Project the helicity moments to the diffent axes
    minnloZcoeff_ptVgen = make_Ais_for_observable("ptVgen", minnloZHelMom, [i for i in range(0, 100, 10)])
    minnloZcoeff_massVgen = make_Ais_for_observable("massVgen", minnloZHelMom)



    for i in range(-1,8):
        plot_Ai_from_hist_for_observable(i, minnloZcoeff_ptVgen, "$p_{T, V}$", f"/eos/home-f/fvazzole/TheoryAgnostic/AngularCoefficients/A_{i}_vs_pT_minnlo.png")
        plot_Ai_from_hist_for_observable(i, minnloZcoeff_massVgen, "$m_V$", f"/eos/home-f/fvazzole/TheoryAgnostic/AngularCoefficients/A_{i}_vs_mV_minnlo.png")

if __name__ == "__main__":
    main()
