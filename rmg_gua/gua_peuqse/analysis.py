import os
import sys
import yaml
from multiprocessing import Pool
import glob
import matplotlib.pyplot as plt
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath("")))))
repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(""))))

def plot_peuquse_progress(results_path):
    """
    plots the rmse 
    """
    proj_path = os.path.dirname(results_path)
    output_files = glob.glob(os.path.join(proj_path, "peuqse_out", "*out.yaml"))
    max_i = -1
    max_val = 0
    max_proc = ""
    
    # load the experimental data  
    expt_path = os.path.join(proj_path, "config", "ct_expt_list.yaml")
    with open(expt_path, "r") as f: 
        expt_data = yaml.safe_load(f)
    
    meohlist = []
    colist = []
    co2list = []
    h2list = []
    h2olist = []
    for data in expt_data:
        meohlist.append(data["species_out"]["CH3OH"])
        colist.append(data["species_out"]["CO"])
        co2list.append(data["species_out"]["CO2"])
        h2list.append(data["species_out"]["H2"])
        h2olist.append(data["species_out"]["H2O"])

    def get_max_methanol(file):
        max_i = -1
        max_val = 0
        min_err = 1e6
        max_proc = ""
        proc = None
        loaded = False
        while not loaded:
            with open(file, "r") as f:
                proc = yaml.load(f, Loader=yaml.CLoader)
            if proc and len(proc) > 0:
                loaded = True
                print(f"iterations = {len(proc)}")


        for idx, entry in enumerate(proc):
            try:
                # evaluate error instead
                error_val = abs(entry[1][0][0] - meohlist[0])
                
                if error_val < min_err:
                    # print(entry[1][0])
                    max_i = idx
                    max_val = entry[1][0][0]
                    err = error_val
                    max_proc = file

            except IndexError:
                print("entry not completed yet")
        print(err)
        return {"max_index":max_i, 
                "max_value":max_val, 
                "error":err,
                "file":max_proc,}

    pool = Pool()
    results = pool.map(get_max_methanol, output_files)
    best_run = min(results, key=lambda x: list(x.values())[2])
    best_err = max(results, key=lambda x: list(x.values())[1])

    with open(best_run["file"], "r") as f: 
        best_data = yaml.load(f, Loader=yaml.CLoader)

    x = []
    y = []

    for idx, dat in enumerate(best_data):
        # get x as run # and y as sse
        meohout = np.array(dat[1][0])
        rss = np.sum(np.array(meohlist) - meohout)**2
        y.append(rss)
        x.append(idx)

    plt.scatter(x,y)


def load_perts_results(pert_list):
    """ 
    load the perturbations to a dictionary
    """
    ck_rule_dict_path = os.path.join(results_path, "ck_rule_dict.yaml")
    with open(ck_rule_dict_path, 'r') as f:
        ck_rule_dict = yaml.load(f, Loader=yaml.FullLoader)

    rule_config_file = os.path.join(results_path, "rule_config.yaml")
    with open(rule_config_file, "r") as f:
        rule_dict_orig = yaml.safe_load(f)

    # # load in the perturbations. haven't thought of a more intelligent way to do this, 
    # but order is preserved in yaml for dict for python 3.6+ if we safe load 
    # and save with sort_keys = False
    rule_dict = deepcopy(rule_dict_orig)
    ck_matches = {}
    ck_no_match = []
    count = 0

    if len(pert_list) > 0:
        assert len(pert_list) == len(rule_dict_orig)*3 + 5, "number of reactions in rxn_list does not match number of reactions in rule_dict_orig"
        print("changing reactions")
        for rule, vals in rule_dict_orig.items():
            rule_dict[rule]["A"] = pert_list[count]
            count += 1
            rule_dict[rule]["E0"] = pert_list[count]
            count += 1
            rule_dict[rule]["alpha"] = pert_list[count]
            count += 1
        # load the binding energy perturbations last values in dict
        thermo_pert_dict = {'C': 0., 'O': 0., 'N': 0., 'H': 0., 'vdw': 0.}
        for num, (key, value) in enumerate(thermo_pert_dict.items()):
            thermo_pert_dict[key] = pert_list[count]
            count += 1


    return rule_dict, thermo_pert_dict