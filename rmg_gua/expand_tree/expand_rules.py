# get the specific node for the sensitive reaction. find out the highest populated parent node 
import pickle
import os
from copy import deepcopy
from rmgpy.data.base import Entry
from rmgpy.reaction import Reaction
from rmgpy.data.kinetics import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.kinetics import StickingCoefficientBEP, StickingCoefficient, SurfaceArrheniusBEP, SurfaceArrhenius


def make_new_rule(model_path, rmg_db_folder, sens_rxn=None):
    """
    make a new rule for the most sensitive reaction.
    model_path: path to the model folder
    rmg_db_folder: path to the rmg database folder
    sens_rxn: the most sensitive reaction, if known 

    if sens_rxn is not specified, then it will search for a pickle file in the model folder. 
    """
    if not sens_rxn:
        cmkn_pickle_path = os.path.join(model_path, "sens_cmkn_dict.pickle")
        with open(cmkn_pickle_path, "rb") as f:
            cmkn_dict = pickle.load(f)
        
        # for now, select the most sensitive reaction. sort dict.
        not_found = True
        while not_found:
            # min_key = min(cmkn_dict, key=cmkn_dict.get)
            min_key = min(cmkn_dict.items(), key=lambda k: k[0][0])[0]
            most_sens_rxn = cmkn_dict[min_key][1]

            # remove the most sensitive reaction from the dict
            cmkn_dict.pop(min_key)

            # verify that the reaction is an estimate
            if most_sens_rxn.kinetics.comment.startswith("Estimated"):
                print("found most sensitive reaction is an estimate")
                not_found = False 
            
            if len(cmkn_dict) == 0:
                print("all sens rxns are from library or training")
                break
    else: 
        most_sens_rxn = sens_rxn

    # load the kinetics database
    # rmg_db_folder= "/Users/blais.ch/Documents/_01_code/RMG_env_1/RMG-database/"
    # Specify the path to the families
    families_dir = os.path.join(rmg_db_folder,"input","kinetics","families")
    if not os.path.exists(families_dir):
        raise OSError(f'Path to rules does not exist:\n{families_dir}')

    # Specify the path to the libraries
    kinetic_libraries_dir = os.path.join(rmg_db_folder,"input","kinetics","libraries","Surface")
    if not os.path.exists(kinetic_libraries_dir):
        raise OSError(f'Path to kinetic libraries does not exist:\n{kinetic_libraries_dir}')

    # do not load training, this will create more rules. in future, may want to 
    # create those rules for perturbation, not sure. 
    kinetics_families = ['all']
    kinetics_database = KineticsDatabase()
    kinetics_database.load_families(
            path=families_dir,
            families=kinetics_families,)


    # breadcrumb cjb need to determine which sensitive reactions we should make rules for. 
    # start with only one. add a test for most sensitive reaction, if a rule already exists, go for second most sensitive.

    # test for reaction: get the nodes that would generate this rxn
    rxn = most_sens_rxn
    family = rxn.get_source()
    template = rxn.template

    # extract info from reaction object for new rule entry
    kinetics_database.families[family].add_atom_labels_for_reaction(rxn)
    template = kinetics_database.families[family].get_reaction_template(rxn)
    source_info = kinetics_database.families[family].extract_source_from_comments(rxn)
    parent = source_info[1][1]["rules"][0][0]
    rank = parent.rank + 1

    # get index for new rule
    last_index = 0
    for entry in kinetics_database.families[family].rules.entries: 
        index = kinetics_database.families[family].rules.entries[entry][0].index
        if index > last_index:
            last_index = index
    new_index = last_index + 1

    data = most_sens_rxn.kinetics

    if isinstance(data, StickingCoefficient):
        data = StickingCoefficientBEP(
            # todo: perhaps make a method StickingCoefficient.StickingCoefficientBEP
            #  analogous to Arrhenius.to_arrhenius_ep
            A=deepcopy(data.A),
            n=deepcopy(data.n),
            alpha=0.5,
            E0=deepcopy(data.Ea),
            Tmin=deepcopy(data.Tmin),
            Tmax=deepcopy(data.Tmax),
            coverage_dependence=deepcopy(data.coverage_dependence),
        )
    elif isinstance(data, SurfaceArrhenius):
        data = SurfaceArrheniusBEP(
            # todo: perhaps make a method SurfaceArrhenius.toSurfaceArrheniusBEP
            #  analogous to Arrhenius.to_arrhenius_ep
            A=deepcopy(data.A),
            n=deepcopy(data.n),
            alpha=0.5,
            E0=deepcopy(data.Ea),
            Tmin=deepcopy(data.Tmin),
            Tmax=deepcopy(data.Tmax),
            coverage_dependence=deepcopy(data.coverage_dependence),
        )
    else:
        logging.error(f"Cannot convert {data} to BEP.")

    new_entry = Entry(
        index=new_index,
        label=';'.join([g.label for g in template]),
        item=Reaction(reactants=[g.item for g in template], products=[]),
        data=data,
        rank=rank,
        short_desc="Rate rule generated for uncertainty",
        long_desc="Rate rule generated for uncertainty",
    )

    try:
        print(f"key already exists for {new_entry.label}, updating entry")
        kinetics_database.families[family].rules.entries[new_entry.label].append(new_entry)
    except KeyError:
        kinetics_database.families[family].rules.entries[new_entry.label] = [new_entry]

    # save the new rules
    kinetics_database.families[family].rules.save(os.path.join(families_dir, family, 'rules' + '.py'))

