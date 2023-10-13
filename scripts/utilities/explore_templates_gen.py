import sys
import argparse
import os
import lz4.frame
import pickle
import hist
from utilities import boostHistHelpers as hh, logging
import hdf5plugin
import h5py
import narf
from narf import ioutils
from utilities import input_tools as it
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
hep.style.use("CMS")
plt.style.use(hep.style.CMS)
plt.rcParams['figure.dpi'] = 50


def main():
    # argument parser for input file and output directory
    parser = argparse.ArgumentParser(description='Process input and output files.')
    parser.add_argument('-f', '--file', type=str, required=True, help='Input file path')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output directory path')
    parser.add_argument('-lb', '--lower_bound', type=float, default=None, required=False, help="Integrating lower bound window of pT")
    parser.add_argument('-ub', '--upper_bound', type=float, default=None, required=False, help="Integrating upper bound window of pT")
    parser.add_argument('-it', '--integrateTheta', type=float, default=False, required=False, help='Integrate Theta')
    parser.add_argument('-ip', '--integratePhi', type=float, default=False, required=False, help='Integrate Phi')
    options = parser.parse_args()

    if options.integrateTheta and options.integratePhi:
        raise ValueError("Cannot integrate both theta and phi!")

    LUMI = '16.8'
    INPUT_FILE = options.file
    OUTPUT_LABEL = options.output + '/{}.png'

    if not os.path.isdir(options.output):
        os.mkdir(options.output)
        print("Made output directory", options.output)

    # load input file
    h5file = h5py.File(INPUT_FILE, "r")
    results = narf.ioutils.pickle_load_h5py(h5file["results"])

    # grab the template
    h = results['ZmumuPostVFP']['output']['nominal_gen'].get()
#    print(h.axes)
    h2 = results['ZmumuPostVFP']['output']['helicity_moments_scale'].get()
#    print(h2.axes)
        
    h = h[::sum, 5j:8j:sum, :, :]
    hep.hist2dplot(h)
    plt.savefig(OUTPUT_LABEL.format('cosTheta_phi'))
    plt.savefig(OUTPUT_LABEL.format('cosTheta_phi').replace(".pdf", ".png"))
    plt.clf()

    hPhi = h.project('phiStarll')
    hep.histplot(hPhi)
    plt.savefig(OUTPUT_LABEL.format('phi'))
    plt.savefig(OUTPUT_LABEL.format('phi').replace(".pdf", ".png"))
    plt.clf()

    hTheta = h.project('cosThetaStarll')
    hep.histplot(hTheta)
    plt.savefig(OUTPUT_LABEL.format('cosTheta'))
    plt.savefig(OUTPUT_LABEL.format('cosTheta').replace(".pdf", ".png"))
    plt.clf()

    # plot the P_i templates
    for hel in range(len(h2.axes[-3])):

        # integrate template across other axes --> yields (cosTheta, phi, pTVgen) for chosen helicity
        h_pi = h2[::sum, 5j:8:sum, :, :, hel, 1j, 1j]
        h_pi_new = hist.Hist(h_pi.axes[1], hist.axis.Regular(10, 0, 2*np.pi, circular = True, name = "phiStarll"), storage = hist.storage.Double())
#        h_pi_new = hist.Hist(h_pi.axes[1], h_pi.axes[0], storage = hist.storage.Double())
        h_pi_new.values()[...] = h_pi.values().T

        original_values = h_pi_new.values()
        new_values = np.copy(original_values)
        N = new_values.shape[1]
        new_values[:, 0:int(N/2)] = original_values[:, int(N/2):int(N)]
        new_values[:, int(N/2):int(N)] = original_values[:, 0:int(N/2)]
        h_pi_new.values()[...] = new_values

        h_pi = h_pi_new
        # if requested, integrate over a particular kinematic angle
        if options.integrateTheta: h_pi = h_pi[::sum, :]
        elif options.integratePhi: h_pi = h_pi[:, ::sum]

        # now plot
        fig = plt.figure()
        ax = fig.subplots()
        if options.integrateTheta or options.integratePhi: _ = hep.histplot(h_pi, ax=ax)
        else: _ = hep.hist2dplot(h_pi, ax=ax)
        hep.cms.label(llabel='Preliminary',data=False, lumi=LUMI, ax=ax)
        label = "P_{}".format(hel-1) if hel != 0 else "UL"
        var_labels = ''
        if not options.integrateTheta: var_labels += r'\cos \theta_{CS}'
        if not options.integrateTheta and not options.integratePhi: var_labels += ', '
        if not options.integratePhi: var_labels += r'\phi_{CS}'
        fig.suptitle("Gen. ${label}$ ( ${var_labels}$ ), $p_T^Z = {lb}-{ub}$ GeV".format(label = label, var_labels = var_labels, lb = options.lower_bound, ub = options.upper_bound), fontsize=24)
        if options.integrateTheta: ax.set_xlabel(r"$\phi_{CS}$")
        elif options.integratePhi: ax.set_xlabel(r"$\cos\theta_{CS}$")
        else:
            ax.set_xlabel(r"$\cos\theta_{CS}$")
            ax.set_ylabel(r"$\phi_{CS}$")
        if options.integrateTheta: label += "_phi"
        elif options.integratePhi: label += "_cosTheta"
        fig.savefig(OUTPUT_LABEL.format(label))
#        fig.savefig(OUTPUT_LABEL.format(label).replace("pdf", "png"))
        print("Saving to", OUTPUT_LABEL.format(label))
#        fig.show()



if __name__ == '__main__':
    main()
