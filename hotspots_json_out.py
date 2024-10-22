from ccdc.protein import Protein
from ccdc.io import MoleculeWriter
from hotspots.hs_io import HotspotWriter
from hotspots.calculation import Runner
import pandas as pd
import re
import json
from time import perf_counter
import os
from pathlib import Path
import argparse
# import glob 

WORKING_DIR = Path(__file__).parent.resolve()
NUM_ROTS = 3000


def find_hotspots(file, args, profile=False):
    pdb_id, file_type = file.split(".")
    
    out_dir = f"{args.input_dir}/results/{NUM_ROTS}/{pdb_id}_{file_type}"

    if not os.path.exists(f"{args.input_dir}/results"):
        os.mkdir(f"{args.input_dir}/results")
    if not os.path.exists(f"{args.input_dir}/results/{NUM_ROTS}"):
        os.mkdir(f"{args.input_dir}/results/{NUM_ROTS}")
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    timings = {"pdb_id": pdb_id,
               "file_type": file_type}
    
    t0 = perf_counter()
    print(f"Begin processing for {pdb_id}...")
    filepath = args.input_dir + f"/{file}"
    prot = Protein.from_file(filepath)

    prot.remove_all_waters()
    prot.add_hydrogens()
    for c in prot.cofactors:
        prot.remove_cofactor(c.identifier)
    for l in prot.ligands:
        prot.remove_ligand(l.identifier)

    # Create lookup dictionaries of chain and residue IDs
    chain_id_lookup = {}
    residue_id_lookup = {}
    for c in prot.chains:
       chain_id_lookup[c.identifier] = c.author_identifier
       residue_id_lookup[c.identifier] = {}
       for r in c.residues:
           # residue identifier is returned as "C:RES123" where "C" is chain ID
           rid = r.identifier.split(":")[-1]
           rid_author = r.author_identifier.split(":")[-1]
           residue_id_lookup[c.identifier][rid] = rid_author

    t1 = perf_counter()
    timings['pre-processing'] = t1-t0
    if profile == True:
        print(f"Pre-processing protein = {t1-t0:.2f} seconds")
        
    save_cavity = "cavity" in args.retain
    
    runner = Runner(pdb_id=pdb_id, output_path=out_dir, save_cavity=save_cavity)

    settings = Runner.Settings()
    settings.nrotations = NUM_ROTS
   
    # Only SuperStar jobs are parallelised (one job per processor). By default there are 3 jobs, when calculating charged interactions there are 5.
    print("file_ID =", pdb_id)
    results = runner.from_protein(prot,
                                  nprocesses=1,
                                  charged_probes=False,
                                  buriedness_method='ghecom',
                                  settings=settings)
    
    t2 = perf_counter()
    timings['calculation'] = t2-t1
    if profile:
        print(f"Whole calculation process = {t2-t1:.1f} seconds")

        # Creates "results/pdb1/out.zip"
    # with HotspotWriter(out_dir) as writer:
    #     writer.write(results)

    # print(f"Calculations for {file_ID} saved to {out_dir}")

    # t4 = perf_counter()
    # if profile == True:
    #     print(f"Save initial results = {t4-t2:.2f} seconds")
    p = results.score(results.protein, as_dict=True)
    t3 = perf_counter()
    timings['scoring'] = t3-t2
    if profile:
        print(f"Scoring = {t3-t2:.1f} seconds")

    df = (pd.DataFrame.from_dict(p, orient='index', columns= ['type','chain_ID','residue_ID','score'])
                      .reset_index()
                      .rename(columns={'index':"atom_ID"})
                      .groupby(['chain_ID','residue_ID'])[['atom_ID','type','score']]
                      .agg(list)                      
                      )

    df.to_csv(f"{out_dir}/{pdb_id}_hotspots.csv")

    t4 = perf_counter()
    timings['create_dataframe'] = t4-t3
    if profile:
        print(f"Create dataframe = {t4-t3:.2f} seconds")

    final_dict = {
        "data_resource": "CSD/Fragment Hotspots",
        "resource_version": "CSD v5.45",
        "software_version": "Hotspots v1.06",
        "resource_entry_url": "https://www.ccdc.cam.ac.uk/open-source-products/fragment-hotspots/",
        "release_date": "20/02/2024",
        "pdb_id": pdb_id,
        "chains": [],
        "sites" : [],
        "evidence_code_ontology": [
        {
        "eco_term": "automatically integrated combinatorial computational and experimental evidence used in automatic assertion",
        "eco_code": "ECO_0007667"
        }
        ]
    }

    for chain, single_chain_df in df.groupby(level=0):
        single_chain_df = single_chain_df.droplevel(0)
        
        author_chain = chain_id_lookup[chain]
        chain_dict = {'chain_label' : author_chain,
                    'residues' : []
                    }
        
        for residue in single_chain_df.index:
            residue_dict = {}
            if residue not in residue_id_lookup[chain]:
                print(f"WARNING: Skipping residue {residue} reported by hotspots as not in lookup table for {chain}")
                continue
            author_residue = residue_id_lookup[chain][residue]
            try:
                residue_dict = {'pdb_res_label' : re.findall("\d+", author_residue)[0],
                                'aa_type' : re.findall("[A-Z]+", author_residue)[0],
                            'site_data' : []
                            }
            except:
                print(f"Error found for {residue}")

            num_sites_in_residue = len(single_chain_df.loc[residue, "atom_ID"])

            try:
                for i in range(num_sites_in_residue):
                    atom_dict = {'site_id_ref' : single_chain_df.loc[residue, "atom_ID"][i],
                                'raw_score' : single_chain_df.loc[residue, "score"][i],
                                'raw_score_unit' : single_chain_df.loc[residue, "type"][i][0],
                                'confidence_classification' : 'null'
                                }
                    
                    site_dict = {'site_id' : single_chain_df.loc[residue, "atom_ID"][i],
                                'label' : single_chain_df.loc[residue, "type"][i][0]
                                }
                    
                    residue_dict['site_data'].append(atom_dict)
                    final_dict['sites'].append(site_dict)
            except:
                print(f"Error found for {residue}")

            chain_dict['residues'].append(residue_dict)
            
        final_dict['chains'].append(chain_dict)
                

    # Serializing json
    json_object = json.dumps(final_dict, indent=4)
    
    # Writing to sample.json
    with open(f"{out_dir}/{pdb_id}_hotspots.json", "w") as outfile:
        outfile.write(json_object)

    print(f'Hotspots for {pdb_id} written to "{out_dir}/{pdb_id}_hotspots.json"')

    t5 = perf_counter()
    if profile:
        print(f"Process JSON - {t5-t4:.2f} seconds")
    timings['total'] = t5-t0

    return timings


if __name__ == '__main__':

    retainable_files = ['mol2', 'cavity']
    
    parser = argparse.ArgumentParser(description='Run Hotspots on all suitable files in a subdirectory')
    parser.add_argument('input_dir',
                        help='the subdirectory containing the input files',
                        type=str)
    parser.add_argument('-s', '--sort',
                        action='store_true',
                        help='sort all files by size and run smallest first')
    parser.add_argument('--retain', 
                        choices=retainable_files,
			default=[],
			nargs='*',
			help='retain intermediate files')
    args = parser.parse_args()

    input_dir = WORKING_DIR / args.input_dir  

    # Storing list of all files (file paths) 
    # in the given directory in list_of_files 

    # Sort list of file names by size  
    list_of_files = filter(lambda x: os.path.isfile 
                        (os.path.join(input_dir, x)), 
                            os.listdir(input_dir) )

    list_of_files = sorted(list_of_files,
                        key =  lambda x: os.stat
                        (os.path.join(input_dir, x)).st_size) 
    
    # Iterate over sorted list of file  
    # names and print them along with size one by one 
     
    for name_of_file in list_of_files: 
        path_of_file = os.path.join(input_dir, name_of_file) 
        size_of_file  = os.stat(path_of_file).st_size  
        print(size_of_file, ' -->', name_of_file)

    try:    
        ref_df = pd.read_csv("processed_files_and_timings.csv")
        print("ref_df loaded")
        id_set = set(ref_df['pdb_id'])

    except FileNotFoundError:
        print("No record of processed files found")
        print("Creating a new file instead")
        ref_df = pd.DataFrame()
        id_set = set()
    
    print(f"{(list_of_files)} to process")
    for file in list_of_files:
        pdb_id, file_type = file.split(".")

        if file_type not in {'cif', 'ent', 'pdb'}:
            print(f"{file} is not of a recognised file type")
            continue

        if pdb_id in id_set:
            print(f"{pdb_id} already processed")
            continue
        
        if 1:#try:
            timings = find_hotspots(file, args, profile=False)
            print(f"\n{file} processed in {timings['total']:.1f} seconds")

            # find_hotspots returns a dict that must be converted to df
            results_list = []
            results_list.append(timings)
            new_df = pd.DataFrame.from_dict(results_list)

            # add results to existing and save latest version
            ref_df = pd.concat([ref_df, new_df], axis=0)
            print(f"{len(ref_df)} of {len(list_of_files)} files complete!\n")
            ref_df.to_csv("processed_files_and_timings.csv", index=False)

        else:#except Exception as error:
            print(f"An error occurred for {file}:", type(error).__name__, "–", error)
