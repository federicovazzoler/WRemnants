import narf
import h5py
import argparse
import hist

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', nargs='+', required=True, type=str, help='List of input files to compare')
    parser.add_argument('-b', default='Z', required=False, help='Which boson')
    return parser.parse_args()

def load_results_from_path(path):
    results = None
    try:
        file = h5py.File(path)
        results = narf.ioutils.pickle_load_h5py(file['results'])
    except:
        print(f'Unable to load hist from {path}')
    return results

def get_helicity_weights_hist(results):
    helWeightHist = hist.Hist(*results.axes, storage=hist.storage.Double())
    helWeightHist.values(flow=True)[...] = helWeightHist.values(flow=True)
    return helWeightHist

def main():
    args = parse_args()
    paths = args.i

    helWeightHists = []
    for path in paths:
        res = load_results_from_path(path)
        helWeightHists.append(get_helicity_weights_hist(res))

    import pdb
    pdb.set_trace()    

if __name__ == "__main__":
    main()