#!/usr/bin/env python3

import os
import sys
from utilities import boostHistHelpers as hh, logging
from utilities.io_tools import input_tools
from wremnants import theory_tools
import h5py
import narf
from narf import ioutils
import ROOT
from wremnants import plot_tools
import matplotlib.pyplot as plt
import mplhep as hep
import pdb
import argparse
import uproot

def load_helicity_moments_for_sample_from_file(sampleName = "ZmumuPostVFP", filePath = None):
    if not os.path.exists(filePath):
        raise Exception(f"{filePath} doesn't exits")

    file = h5py.File(filePath, "r")
    res = narf.ioutils.pickle_load_h5py(file["results"])
    h = input_tools.load_and_scale(res, sampleName, "nominal_gen_helicity_moments_scale")[{"muRfact" : 1.j, "muFfact" : 1.j}].project('massVgen', 'y', 'ptVgen', 'chargeVgen', 'helicity')

    return h

def make_Ais_for_observable(obs, h, rebin = None):
    if rebin is not None:
        h = hh.rebinHist(h, obs, rebin)
    h_coeff = theory_tools.moments_to_angular_coeffs(h.project(obs, "helicity"), sumW2 = False)

    return h_coeff

def plot_hist(h, xlabel = None, ylabel = None, label=None, xrange = None, yrange = None, corr = 1, ax=None):

    if ax is None:
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (8, 8), sharex = True)

    xValues = h.axes.centers[0]
    yValues = h.values() * corr
#    yVariances = h.variances()

    ax.plot(xValues, yValues, marker = '.', linestyle = 'dotted', label=label)
#    plt.errorbar(xValues, yValues, yerr=yVariances, marker = ".", linestyle = "solid")

    if xlabel:
        ax.set_xlabel(xlabel, fontsize = 22)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize = 22)
    ax.set_xscale("log")
    if xrange is not None:
        ax.set_xlim(xrange)
    if yrange is not None:
        ax.set_ylim(yrange)

    try:
        return fig, ax
    except NameError:
        return ax

def save_fig(fig, outname):
    fig.set_tight_layout(True)
    fig.savefig(outname)
    print(f"{outname} created")

def get_dyturbo_Ai(Ai, dyturbo_fname):
    file = uproot.open(dyturbo_fname)
    h = file[f"wgt_a{Ai}_y_qt"].to_hist().project("xaxis")
    hs = file['s_qt'].to_hist()
    h = hh.divideHists(h, hs)
    return h

def get_atlas_Ai(Ai, atlas_dname):
    # from https://www.hepdata.net/record/76986
    fname = f"HEPData-ins1466778-v1-Table_{Ai+2}.root"
    file = uproot.open(atlas_dname+"/"+fname)
    h = file['Table '+str(Ai + 2)]['Hist1D_y1'].to_hist()
    return h

def plot_Ai_from_hist_for_observable(Ai, h, obs, outname, dyturbo_fname=None, atlas_dname=None):
    hToPlot = h[{"helicity" : complex(Ai)}]
    corr = 1
    # if Ai == 4 or Ai == 1 or Ai == 3: corr = -1
    fig, ax = plot_hist(hToPlot, xlabel = obs, ylabel = f"$A_{Ai}$", label="minlo", xrange = None, yrange = None, corr = corr)

    if dyturbo_fname and os.path.exists(dyturbo_fname) and Ai != -1:
        if obs != "$p_{T, V}$":
            raise Exception("Dyturbo only implemented for pT")
        h_dyturbo_Ai = get_dyturbo_Ai(Ai, dyturbo_fname)
        ax = plot_hist(h_dyturbo_Ai, label="dyturbo", ax=ax)

    if atlas_dname and os.path.isdir(atlas_dname) and Ai != -1:
        if obs != "$p_{T, V}$":
            raise Exception("ATLAS only implemented for pT")
        h_atlas_Ai = get_atlas_Ai(Ai, atlas_dname)
        ax = plot_hist(h_atlas_Ai, label="ATLAS", ax=ax)

    ax.legend()
    save_fig(fig, outname)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help = "Input file name")
    parser.add_argument("-o", "--outdir", help = "Output directory")
    parser.add_argument("-d", "--dyturbo", help = "Dyturbo file name", default=None, required=False)
    parser.add_argument("-a", "--atlas", help = "ATLAS directory name", default=os.path.join(os.getenv("WREM_HELPERS"), "data"), required=False)
    args = parser.parse_args()

    filePath = f"{args.file}"

    # Load the helicity moments
    minnloZHelMom = load_helicity_moments_for_sample_from_file(sampleName = "ZmumuPostVFP", filePath = filePath)
    minnloWmHelMom = load_helicity_moments_for_sample_from_file(sampleName = "WminusmunuPostVFP", filePath = filePath)
    minnloWpHelMom = load_helicity_moments_for_sample_from_file(sampleName = "WplusmunuPostVFP", filePath = filePath)

    # Project the helicity moments to the diffent axes
    minnloZcoeff_ptVgen = make_Ais_for_observable("ptVgen", minnloZHelMom, [0.00e+00, 1.00e+00, 2.00e+00, 3.00e+00, 4.00e+00, 5.00e+00, 6.00e+00, 7.00e+00, 8.00e+00, 9.00e+00, 1.00e+01, 1.10e+01, 1.20e+01, 1.30e+01, 1.40e+01, 1.50e+01, 1.60e+01, 1.70e+01, 1.80e+01, 1.90e+01, 2.00e+01, 2.10e+01, 2.20e+01, 2.30e+01, 2.40e+01, 2.50e+01, 2.60e+01, 2.70e+01, 2.80e+01, 2.90e+01, 3.00e+01, 1.00e+02, 2.20e+02, 1.30e+04])
    minnloZcoeff_massVgen = make_Ais_for_observable("massVgen", minnloZHelMom)


    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    for i in range(-1,8):
        plot_Ai_from_hist_for_observable(i, minnloZcoeff_ptVgen, "$p_{T, V}$", f"{args.outdir}/A_{i}_vs_pT_minnlo.png", dyturbo_fname=args.dyturbo, atlas_dname=args.atlas)
        plot_Ai_from_hist_for_observable(i, minnloZcoeff_massVgen, "$m_V$", f"{args.outdir}/A_{i}_vs_mV_minnlo.png")

if __name__ == "__main__":
    main()