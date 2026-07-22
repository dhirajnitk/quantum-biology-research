#!/usr/bin/env python3
"""
This script wil generate the gitlab CI yaml
"""

import argparse
import sys

import json
import yaml
import os
import glob
import copy

def fname_without_extension(fname):
    return ".".join(fname.split(".")[0:-1])

def fname_without_path(fname):
    return stripped_fname.split("/")[-1]

def fname_without_path_or_extension(fname):
    stripped_fname = ".".join(fname.split(".")[0:-1])
    without_path = stripped_fname.split("/")[-1]
    return without_path

def ci_diagonalize_matrix_job_template():
    return {
        "extends": [".diagonalize-matrix"],
        "variables": {},
    }

def ci_create_matrix_job_template():
    return {
        "extends": [".create-matrix"],
        "variables": {},
    }

def ci_combine_realizations_job_template():
    return {
        "extends": [".combine-realizations"],
        "variables": {},
    }

def ci_create_plot_data_job_template():
    return {
        "extends": [".create-plot-data"],
        "variables": {},
    }

def ci_prepare_pipelines_job_template():
    return {
        "extends": [".superradiance"],
        "variables": {},
    }

def get_base_yaml(filename):
    with open(filename, "r") as f:
        return yaml.safe_load(f)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base_yaml",
        required=True,
        help="A YAML file containing the base structure for the pipeline",
    )
    parser.add_argument(
        "--noise_values",
        nargs='+',
        required=True,
        help="The noise values",
    )
    parser.add_argument(
        "--num_realizations",
        type=int,
        required=True,
        help="A YAML file containing the base structure for the pipeline",
    )
    return parser.parse_args()

def main():
    args = parse_args()

    base_yaml = get_base_yaml(args.base_yaml)
    base_yaml["variables"] = {}
    base_yaml["variables"]["PLOT_ONLY"] = os.environ.get("PLOT_ONLY")
    # If this variable is not "", then it overrides the previous 2
    dir_to_scan = os.environ.get("DIR_TO_SCAN")
    base_yaml["variables"]["DIR_TO_SCAN"] = dir_to_scan
    base_yaml["variables"]["LOCAL_MAX_DECAY_RATE_DIRECTORY"] = os.environ.get("PARENT_LOCAL_MAX_DECAY_RATE_DIRECTORY")
    base_yaml["variables"]["LOCAL_NPY_MATRICES_DIRECTORY"] = os.environ.get("PARENT_LOCAL_NPY_MATRICES_DIRECTORY")
    base_yaml["variables"]["LOCAL_QUANTUM_YIELD_DIRECTORY"] = os.environ.get("PARENT_LOCAL_QUANTUM_YIELD_DIRECTORY")

    all_create_matrix_yamls = []

    if dir_to_scan != "":

        # Gets all files in directory
        files_in_dir_to_scan = glob.glob(f"{dir_to_scan}/*")

        # Loops through all files in directory, and creates a create-matrix job for each one
        for filename in files_in_dir_to_scan:

            create_matrix_yaml = {}
            create_matrix_job_name = f"create-matrix:{fname_without_path_or_extension(filename)}"

            create_matrix_yaml[create_matrix_job_name] = ci_create_matrix_job_template()
            create_matrix_yaml[create_matrix_job_name]["variables"]["MAX_DECAY_RATE_DIR"] = "${LOCAL_MAX_DECAY_RATE_DIRECTORY}/" + fname_without_path_or_extension(filename)
            create_matrix_yaml[create_matrix_job_name]["variables"]["QUANTUM_YIELD_DIR"] = "${LOCAL_QUANTUM_YIELD_DIRECTORY}/" + fname_without_path_or_extension(filename)
            create_matrix_yaml[create_matrix_job_name]["variables"]["POS_DIPOLE_FILENAME"] = filename
            create_matrix_yaml[create_matrix_job_name]["variables"]["NPY_MATRIX_FILE"] = "${LOCAL_NPY_MATRICES_DIRECTORY}/" + fname_without_path_or_extension(filename)

            all_create_matrix_yamls.append(create_matrix_yaml)

    else:
        base_yaml["variables"]["POS_DIPOLE_NAME"] = os.environ.get("POS_DIPOLE_FILE")
        base_yaml["variables"]["POS_DIPOLE_FILENAME"] = f"{os.environ.get('POS_DIPOLE_DIR')}/{os.environ.get('POS_DIPOLE_FILE')}.txt"
        base_yaml["variables"]["NPY_MATRIX_FILE"] = f"{os.environ.get('LOCAL_NPY_MATRICES_DIRECTORY')}/{os.environ.get('POS_DIPOLE_FILE')}"
        base_yaml["variables"]["MAX_DECAY_RATE_DIR"] = f"{os.environ.get('LOCAL_MAX_DECAY_RATE_DIRECTORY')}/{os.environ.get('POS_DIPOLE_FILE')}"
        base_yaml["variables"]["QUANTUM_YIELD_DIR"] = f"{os.environ.get('LOCAL_QUANTUM_YIELD_DIRECTORY')}/{os.environ.get('POS_DIPOLE_FILE')}"
        # This code was already there
        pos_dipole_type = os.environ.get("POS_DIPOLE_FILE")
        create_matrix_yaml = {}
        create_matrix_job_name = f"create-matrix:{pos_dipole_type}"
        create_matrix_yaml[create_matrix_job_name] = ci_create_matrix_job_template()
        all_create_matrix_yamls.append(create_matrix_yaml)

    all_diagonalize_matrix_yamls = []
    all_combine_realizations_yamls = []
    all_create_plot_data_yamls = []

    if base_yaml["variables"]["DIR_TO_SCAN"] != "":
        for filename in files_in_dir_to_scan:
            diagonalize_matrix_yaml = {}
            combine_realizations_yaml = {}
            create_plot_data_yaml = {}
            previous_job_name = None
            for noise in args.noise_values:
                plot_rule = {}
                plot_rule["if"] = '$PLOT_ONLY == "false"'

                diagonalize_job_name = f"diagonalize-matrix:{fname_without_path_or_extension(filename)}:noise-{noise}"
                diagonalize_matrix_yaml[diagonalize_job_name] = ci_diagonalize_matrix_job_template()
                diagonalize_matrix_yaml[diagonalize_job_name]["variables"]["NOISE"] = noise
                diagonalize_matrix_yaml[diagonalize_job_name]["variables"]["MAX_DECAY_RATE_DIR"] = "${LOCAL_MAX_DECAY_RATE_DIRECTORY}/" + fname_without_path_or_extension(filename)
                diagonalize_matrix_yaml[diagonalize_job_name]["variables"]["QUANTUM_YIELD_DIR"] = "${LOCAL_QUANTUM_YIELD_DIRECTORY}/" + fname_without_path_or_extension(filename)
                diagonalize_matrix_yaml[diagonalize_job_name]["variables"]["NPY_MATRIX_FILE"] = "${LOCAL_NPY_MATRICES_DIRECTORY}/" + fname_without_path_or_extension(filename)

                diagonalize_matrix_yaml[diagonalize_job_name]["rules"] = []
                diagonalize_matrix_yaml[diagonalize_job_name]["rules"].append(plot_rule)

                if noise == "0":
                    diagonalize_matrix_yaml[diagonalize_job_name]["variables"]["CI_NODE_INDEX"] = 0
                else:
                    diagonalize_matrix_yaml[diagonalize_job_name]["parallel"] = args.num_realizations

                # The combine realization jobs depende on the associated noise diagonalize jobs
                
                realizations_job_name = f"combine-realizations:{fname_without_path_or_extension(filename)}:noise-{noise}"
                combine_realizations_yaml[realizations_job_name] = ci_combine_realizations_job_template()
                combine_realizations_yaml[realizations_job_name]["variables"]["NOISE"] = noise
                combine_realizations_yaml[realizations_job_name]["variables"]["MAX_DECAY_RATE_DIR"] = "${LOCAL_MAX_DECAY_RATE_DIRECTORY}/" + fname_without_path_or_extension(filename)
                combine_realizations_yaml[realizations_job_name]["variables"]["QUANTUM_YIELD_DIR"] = "${LOCAL_QUANTUM_YIELD_DIRECTORY}/" + fname_without_path_or_extension(filename)

                combine_realizations_yaml[realizations_job_name]["needs"] = [diagonalize_job_name]
                combine_realizations_yaml[realizations_job_name]["rules"] = []
                combine_realizations_yaml[realizations_job_name]["rules"].append(plot_rule)

            create_plot_data_job_name = f"create-plot-data:{fname_without_path_or_extension(filename)}"
            create_plot_data_yaml[create_plot_data_job_name] = ci_create_plot_data_job_template()
            create_plot_data_yaml[create_plot_data_job_name]["variables"]["MAX_DECAY_RATE_DIR"] = "${LOCAL_MAX_DECAY_RATE_DIRECTORY}/" + fname_without_path_or_extension(filename)
            create_plot_data_yaml[create_plot_data_job_name]["variables"]["QUANTUM_YIELD_DIR"] = "${LOCAL_QUANTUM_YIELD_DIRECTORY}/" + fname_without_path_or_extension(filename)
            create_plot_data_yaml[create_plot_data_job_name]["variables"]["POS_DIPOLE_NAME"] = fname_without_path_or_extension(filename)

            all_diagonalize_matrix_yamls.append(diagonalize_matrix_yaml)
            all_combine_realizations_yamls.append(combine_realizations_yaml)
            all_create_plot_data_yamls.append(create_plot_data_yaml)

        prepare_pipelines_yaml = get_base_yaml(".gitlab-ci-prepare-pipelines.yml")
        prepare_pipelines_yaml["variables"] = {}
        prepare_pipelines_yaml["variables"]["PLOT_ONLY"] = os.environ.get("PLOT_ONLY")
        prepare_pipelines_yaml["variables"]["DIR_TO_SCAN"] = os.environ.get("DIR_TO_SCAN")
        prepare_pipelines_yaml["variables"]["LOCAL_MAX_DECAY_RATE_DIRECTORY"] = os.environ.get("PARENT_LOCAL_MAX_DECAY_RATE_DIRECTORY")
        prepare_pipelines_yaml["variables"]["LOCAL_QUANTUM_YIELD_DIRECTORY"] = os.environ.get("PARENT_LOCAL_QUANTUM_YIELD_DIRECTORY")
        prepare_pipelines_yaml["variables"]["LOCAL_NPY_MATRICES_DIRECTORY"] = os.environ.get("PARENT_LOCAL_NPY_MATRICES_DIRECTORY")

        for filename in files_in_dir_to_scan:
            prepare_pipelines_job_name = f"superradiance:{fname_without_path_or_extension(filename)}"
            prepare_pipelines_job_artifact_name = f"gitlab-ci-calculate-superradiance-{fname_without_path_or_extension(filename)}.yml"

            prepare_pipelines_yaml[prepare_pipelines_job_name] = copy.deepcopy(prepare_pipelines_yaml[".superradiance"])
            prepare_pipelines_yaml[prepare_pipelines_job_name]["trigger"]["include"][0]["artifact"] = prepare_pipelines_job_artifact_name

            calculate_superradiance_yaml_name = f"gitlab-ci-calculate-superradiance-{fname_without_path_or_extension(filename)}.yml"
            calculate_superradiance_yaml = get_base_yaml(".gitlab-ci-calculate-superradiance.yml")
            calculate_superradiance_yaml["variables"] = {}
            calculate_superradiance_yaml["variables"]["LOCAL_MAX_DECAY_RATE_DIRECTORY"] = os.environ.get("LOCAL_MAX_DECAY_RATE_DIRECTORY")
            calculate_superradiance_yaml["variables"]["LOCAL_QUANTUM_YIELD_DIRECTORY"] = os.environ.get("LOCAL_QUANTUM_YIELD_DIRECTORY")
            calculate_superradiance_yaml["variables"]["LOCAL_NPY_MATRICES_DIRECTORY"] = os.environ.get("LOCAL_NPY_MATRICES_DIRECTORY")

            with open(calculate_superradiance_yaml_name, 'x') as yaml_name:
                ind = files_in_dir_to_scan.index(filename)
                yaml.dump(calculate_superradiance_yaml, yaml_name)
                yaml.dump(all_create_matrix_yamls[ind], yaml_name)
                yaml.dump(all_diagonalize_matrix_yamls[ind], yaml_name)
                yaml.dump(all_combine_realizations_yamls[ind], yaml_name)
                yaml.dump(all_create_plot_data_yamls[ind], yaml_name)

        with open("gitlab-ci-prepare-pipelines.yml", 'x') as yaml_name:
            yaml.dump(prepare_pipelines_yaml, yaml_name)

    else:
        diagonalize_matrix_yaml = {}
        combine_realizations_yaml = {}
        create_plot_data_yaml = {}
        previous_job_name = None
        for noise in args.noise_values:
            plot_rule = {}
            plot_rule["if"] = '$PLOT_ONLY == "false"'

            diagonalize_job_name = f"diagonalize-matrix:noise-{noise}"
            diagonalize_matrix_yaml[diagonalize_job_name] = ci_diagonalize_matrix_job_template()
            diagonalize_matrix_yaml[diagonalize_job_name]["variables"]["NOISE"] = noise
            diagonalize_matrix_yaml[diagonalize_job_name]["rules"] = []
            diagonalize_matrix_yaml[diagonalize_job_name]["rules"].append(plot_rule)
            if noise == "0":
                diagonalize_matrix_yaml[diagonalize_job_name]["variables"]["CI_NODE_INDEX"] = 0
            else:
                diagonalize_matrix_yaml[diagonalize_job_name]["parallel"] = args.num_realizations
            # Depend on the previous noise job.
            if previous_job_name:
                diagonalize_matrix_yaml[diagonalize_job_name]["needs"] = [previous_job_name]
            previous_job_name = diagonalize_job_name
            #
            # The combine realization jobs depende on the associated noise diagonalize jobs
            #
            realizations_job_name = f"combine-realizations:noise-{noise}"
            combine_realizations_yaml[realizations_job_name] = ci_combine_realizations_job_template()
            combine_realizations_yaml[realizations_job_name]["variables"]["NOISE"] = noise
            combine_realizations_yaml[realizations_job_name]["needs"] = [diagonalize_job_name]
            combine_realizations_yaml[realizations_job_name]["rules"] = []
            combine_realizations_yaml[realizations_job_name]["rules"].append(plot_rule)

        create_plot_data_job_name = f"create-plot-data"
        create_plot_data_yaml[create_plot_data_job_name] = ci_create_plot_data_job_template()

        all_diagonalize_matrix_yamls.append(diagonalize_matrix_yaml)
        all_combine_realizations_yamls.append(combine_realizations_yaml)
        all_create_plot_data_yamls.append(create_plot_data_yaml)

if __name__ == "__main__":
    main()

