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
	options = parser.parse_args()

	LUMI = 'xx'
	INPUT_FILE = options.file
	OUTPUT_LABEL = options.output + '/{}.pdf'

	if not os.path.isdir(options.output):
		os.mkdir(options.output)
		print("Made output directory", options.output)

	# load input file
	h5file = h5py.File(INPUT_FILE, "r")
	results = narf.ioutils.pickle_load_h5py(h5file["results"])

	# grab the template
	h = results['ZmumuPostVFP']['output']['nominal'].get()

	print("Found the following axes:")
	for ax in h.axes:
		print("\t", ax.name)
		
	# plot the P_i templates
	for hel in range(len(h.axes[-1])):

		# integrate template across other axes --> yields (cosTheta, phi, pTVgen, helicity)
		h_pi = h[::sum,::sum,:,:,::sum,:, hel]

		# and over specific pTVgen window
		if options.lower_bound is None: options.lower_bound = h_pi.axes[-2].edges[0]
		if options.upper_bound is None: options.upper_bound = h_pi.axes[-2].edges[-1]
		h_pi = h_pi[:,:,options.lower_bound*1.0j:options.upper_bound*1.0j:sum]

		# now plot
		fig = plt.figure()
		ax = fig.subplots()
		_ = hep.hist2dplot(h_pi, ax=ax)
		hep.cms.label(llabel='Preliminary',data=False, lumi=LUMI, ax=ax)
		label = "P_{}".format(hel-1) if hel != 0 else "UL"
		fig.suptitle("Templated ${}$".format(label) + \
			r"($\cos \theta_{CS}, \phi_{CS}$), " + \
			"$p_T^Z = {}-{}$ GeV".format(options.lower_bound, options.upper_bound),
			fontsize=24)
		ax.set_xlabel(r"$\cos\theta_{CS}$")
		ax.set_ylabel(r"$\phi_{CS}$")
		fig.savefig(OUTPUT_LABEL.format(label))
		fig.savefig(OUTPUT_LABEL.format(label).replace("pdf", "png"))
		print("Saving to", OUTPUT_LABEL.format(label))
		fig.show()


if __name__ == '__main__':
	main()