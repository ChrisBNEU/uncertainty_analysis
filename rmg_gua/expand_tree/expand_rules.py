# get the specific node for the sensitive reaction. find out the highest populated parent node 
import pickle
import os
from copy import deepcopy
from rmgpy.quantity import ScalarQuantity
from rmgpy.data.base import Entry
from rmgpy.reaction import Reaction
from rmgpy.data.kinetics import KineticsDatabase
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.kinetics import StickingCoefficientBEP, StickingCoefficient, SurfaceArrheniusBEP, SurfaceArrhenius

def generate_thermo(rxn):
    for reac in rxn.reactants:
        reac.thermo = thermo_database.get_thermo_data(reac,metal_to_scale_to="Cu")
    for prod in rxn.products:
        prod.thermo = thermo_database.get_thermo_data(prod,metal_to_scale_to="Cu")

    return rxn

def get_e0(rxn, alpha=0.5):

    dhrxn = rxn.get_enthalpy_of_reaction(298)
    Ea = rxn.kinetics.Ea.value_si
    E0 = Ea - 0.5*dhrxn
    
    E0_quantity = ScalarQuantity(value=E0, units='J/mol',)
    return E0_quantity

def make_new_rule(model_path, rmg_db_folder, sens_rxn=None, update_old=False):
    """
    make a new rule for the most sensitive reaction.
    model_path: path to the model folder
    rmg_db_folder: path to the rmg database folder
    sens_rxn: the most sensitive reaction, if known 

    if sens_rxn is not specified, then it will search for a pickle file in the model folder. 
    """
    # do not make a new rule unless we have to. 
    new_rule = False

    if not sens_rxn:
        cmkn_pickle_path = os.path.join(model_path, "sens_cmkn_dict.pickle")
        with open(cmkn_pickle_path, "rb") as f:
            cmkn_dict = pickle.load(f)
        
        # for now, select the most sensitive reaction. sort dict.
        not_found = True
        while not_found:
            # min_key = min(cmkn_dict, key=cmkn_dict.get)
            min_key = min(cmkn_dict.items(), key=lambda k: k[0][0])[0]
            rxn = cmkn_dict[min_key][1]

            # remove the most sensitive reaction from the dict
            cmkn_dict.pop(min_key)

            # verify that the reaction is an estimate
            if rxn.kinetics.comment.startswith("Estimated"):
                print("found most sensitive reaction is an estimate")
                not_found = False 
            
            if len(cmkn_dict) == 0:
                print("all sens rxns are from library or training")
                break
    else: 
        rxn = sens_rxn

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
    global kinetics_database
    kinetics_families = ['all']
    kinetics_database = KineticsDatabase()
    kinetics_database.load_families(
            path=families_dir,
            families=kinetics_families,)

    # only do this if there is no thermo for rxn
    if rxn.reactants[0].thermo is None:
        global thermo_database
        thermo_database = ThermoDatabase()
        tdb_path = os.path.join(rmg_db_folder, "input", "thermo")
        thermo_database.load(path = tdb_path, libraries = "surfaceThermoPt111", surface=True)

        # match all of the reactants and products for the reaction
        rxn = generate_thermo(rxn)

    # breadcrumb cjb need to determine which sensitive reactions we should make rules for. 
    # start with only one. add a test for most sensitive reaction, if a rule already exists, go for second most sensitive.

    # test for reaction: get the nodes that would generate this rxn


    # allow user to supply template if it is a library reaction
    if isinstance(rxn, TemplateReaction):
        template = rxn.template
        family = rxn.get_source()

        kinetics_database.families[family].add_atom_labels_for_reaction(rxn)
        template = [item.label for item in kinetics_database.families[family].get_reaction_template(rxn)]

        try: 
            rank = kinetics_database.families[family].rules.entries[";".join(template)][0].rank
            source_info = kinetics_database.families[family].extract_source_from_comments(rxn)
            parent = source_info[1][1]["rules"][0][0]
            rank = parent.rank
        except KeyError:
            new_rule = True
            rank = kinetics_database.families[family].rules.entries[";".join(template)][0].rank
            

        # get E0 for reaction
        data_E0 = get_e0(rxn, alpha=rxn.kinetics.alpha.value_si)

    # if its a library reaction, generate a family reaction from it. 
    elif isinstance(rxn, LibraryReaction): 

        try: 
            fam_rxn_list = kinetics_database.generate_reactions_from_families(
                reactants=rxn.reactants,
                products=rxn.products,
                only_families=None,
                resonance=True,
            )
            if len(fam_rxn_list) == 1:
                fam_rxn = fam_rxn_list[0]
                family = fam_rxn_list[0].get_source()
                # template = fam_rxn_list[0].template
                kinetics_database.families[family].add_atom_labels_for_reaction(fam_rxn)
                template_obj = kinetics_database.families[family].get_reaction_template(fam_rxn)
            else: 
                raise Exception(f"{len(fam_rxn_list)} matches for library reaction, requires manual selection")
            
            kinetics_database.families[family].add_atom_labels_for_reaction(fam_rxn)
            template = [item.label for item in template_obj]
            
            try: 
                print((template))
                rank = kinetics_database.families[family].rules.entries[";".join(template)][0].rank
            except KeyError:
                new_rule = True
                rank = 0
            E0_data = get_e0(rxn, alpha=0.5)

        except Exception as e:
            raise Exception("could not generate family reaction from library reaction, error: ", str(e))
    else: 
        raise Exception("must be either template or library reaction")


    # get index for new rule, or use same index for old rule
    if new_rule:
        last_index = 0
        for entry in kinetics_database.families[family].rules.entries: 
            index = kinetics_database.families[family].rules.entries[entry][0].index
            if index > last_index:
                last_index = index
        new_index = last_index + 1
    else: 
        new_index = kinetics_database.families[family].rules.entries[";".join(template)][0].index

    data = rxn.kinetics
    if len(rxn.comment) > 0:
        desc = rxn.comment

    if isinstance(data, StickingCoefficient):
        data = StickingCoefficientBEP(
            # todo: perhaps make a method StickingCoefficient.StickingCoefficientBEP
            #  analogous to Arrhenius.to_arrhenius_ep
            A=deepcopy(data.A),
            n=deepcopy(data.n),
            alpha=0.5,
            E0=deepcopy(E0_data),
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
            E0=deepcopy(E0_data),
            Tmin=deepcopy(data.Tmin),
            Tmax=deepcopy(data.Tmax),
            coverage_dependence=deepcopy(data.coverage_dependence),
        )
    else:
        logging.error(f"Cannot convert {data} to BEP.")

    new_entry = Entry(
        index=new_index,
        label=';'.join(template),
        item=Reaction(reactants=[g.item for g in template_obj], products=[]),
        data=data,
        rank=rank,
        short_desc="Rate rule generated for uncertainty",
        long_desc=f"{desc}",
    )


    try:
        print(f"key already exists for {new_entry.label}, updating entry")
        if update_old:
            kinetics_database.families[family].rules.entries[new_entry.label] = [new_entry]
        else: 
            print("leaving entry as is. to update, set update_old=True")
    except KeyError:
        kinetics_database.families[family].rules.entries[new_entry.label].append(new_entry)

    # save the new rules
    kinetics_database.families[family].rules.save(os.path.join(families_dir, family, 'rules' + '.py'))

