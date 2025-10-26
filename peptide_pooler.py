import os
import sys
from datetime import datetime

import yaml
import numpy as np
import pandas as pd

import robotools
from loguru import logger


def arg_parser():
    parser = argparse.ArgumentParser(description="Peptide pooler")
    parser.add_argument("--config", type=str, default="config.yaml", help="Path to config file")
    parser.add_argument("--input", type=str, help="Path to input data file")
    parser.add_argument("--outdir", type=str, help="Path to output directory")
    return parser.parse_args()

def parse_trays(labware):
    """
    Define the positions of indexed peptides in trays
    Args:
        labware: dict
            trays: list of dicts
                name: tray name
                samples: list of peptide indices
    Returns:
        peptides_positions: dict
            key: peptide index
            value: dict
                tray: tray name
                well: well name
    """
    peptides_positions = dict()
    for tray in labware["trays"]:
        for well, sample in enumerate(tray["samples"]):
            peptides_positions[sample] = {"tray": tray["name"], "well": "A" + str(well + 1).zfill(2)}
    return peptides_positions

def parse_delution_plate(labware):
    """
    Define the positions of indexed peptides in delution plate
    Args:
        labware: dict
            delution_plate: dict
                name: delution plate name
                samples: list of peptide indices
    Returns:
        delution_plate_positions: dict
        key: peptide index
        value: dict
            delution_plate: delution plate name
            well: well name
    """
    delution_plate_positions = dict()
    wells_96 = [let + str(num) for num in range(1, 13) for let in ["A", "B", "C", "D", "E", "F", "G", "H"]]
    for well, sample in enumerate(labware["delution_plate"]["samples"]):
        delution_plate_positions[sample] = {"tray": labware["delution_plate"]["name"], "well": wells_96[well]}
    return delution_plate_positions

def split_volume(volume, factor):
    

def main(args):
    # read the config file
    with open(args.config, "r") as fr:
        # if the config file was wrongly configured, stop the script
        try:
            config = yaml.load(fr, Loader=yaml.SafeLoader)
        except yaml.YAMLError as e:
            logger.error(f"Error loading config: {e}")
            exit(1)
    # define the launch timestamp according to the rule of the config file
    timestamp = datetime.now().strftime(config["additional"]["date_preset"])
    logger.remove()
    logger.add(sys.stdout, level="INFO")
    # define the output directory
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    logger.add(f"{args.outdir}/{timestamp}_peptide_pooler.log", level="DEBUG")
    logger.debug(f"Config loaded: {config}")
    
    # check the trays and the delution plate
    peptides_positions = parse_trays(config["labware"])
    delution_plate_positions = parse_delution_plate(config["labware"])
    # if the peptides are not found in the delution plate, stop the script
    for element in peptides_positions.keys():
        if element not in delution_plate_positions.keys():
            logger.error(f"Peptide {element} not found in delution plate")
            exit(1)
    # read the input data file
    df = pd.read_excel(args.input, names=["pep_number", "sequence", "conc"])
    logger.debug(f"Input data loaded: {df.head(3)}")
    # iterate over the input data file and check if the indexed peptides are found in the trays and the delution plate
    for index, row in df.iterrows():
        if df['pep_number'] not in peptides_positions.keys() or df['pep_number'] not in delution_plate_positions.keys():
            logger.error(f"Peptide {df['pep_number']}:{df['sequence']} not found in trays or delution plate, check the config file")
            exit(1)
    
    min_volume = config["volumes"]["min_pipette_volume"]
    max_volume = config["volumes"]["max_peptide_usage_volume"]
    max_result_volume = config["volumes"]["max_pool_volume"]

    df.sort_values(by="conc", inplace=True)
    # количество пептидов
    num_of_peptides = df.shape[0]
    # определяем индекс середины списка пептидов
    half_index = num_of_peptides // 2
    # получае половинки DataFrame
    half_df = df[:half_index]
    outer_df = df[half_index:]
    # определяем максимальную концентрацию пептидов на нижней половине
    half_max_conc = half_df['conc'].max()
    # количество пептида в максимальной концентрации при минимальном объеме на нижней половине
    half_amount = half_max_conc * min_volume

    # определяем такие объемы пептидов, которые получатся при минимальном использовании максимального пептида
    delute_volumes = (half_amount / half_df['conc']).round()
    # определяем суммарный объем пептидов на нижней половине
    inner_volume = delute_volumes.sum()

    # определяем суммарный объем пептидов на верхней половине
    outer_volume = outer_df.shape[0] * min_volume

    # определяем, подходит ли суммарный объем пептидов на нижней и верхней половине для пула
    criteria_fitting = (inner_volume + outer_volume) <= max_result_volume
    # определяем, помещается ли максимально дозволенный для отбора объем
    criteria_lowest_fits = delute_volumes.max() < max_volume

    logger.debug(f"Half index: {half_index}, Criteria fitting: {criteria_fitting}, Criteria lowest fits: {criteria_lowest_fits}, Max delute volume: {delute_volumes.max()}")

    # определяем последний успешный индекс
    last_success = half_index

    # если оба критерия выполняются, то станем искать более выгодные варианты
    if criteria_fitting and criteria_lowest_fits:
        # до тех пор, пока оба критерия выполняются, ищем более выгодные варианты
        while criteria_fitting and criteria_lowest_fits:
            if half_index <= df.shape[0] - 1:
                last_success = half_index
                half_index += 1
                half_df = df[:half_index]
                outer_df = df[half_index:]
                half_max_conc = half_df['conc'].max()
                
                half_amount = half_max_conc * min_volume
                
                delute_volumes = (half_amount / half_df['conc']).round()
                inner_volume = delute_volumes.sum()
                
                outer_volume = outer_df.shape[0] * min_volume
                
                criteria_fitting = (inner_volume + outer_volume) <= max_result_volume
                criteria_lowest_fits = delute_volumes.max() < max_volume
                logger.debug(f"Half index: {half_index}, Criteria fitting: {criteria_fitting}, Criteria lowest fits: {criteria_lowest_fits}, Max delute volume: {delute_volumes.max()}")
            else:
                logger.warning(f"Half index is out of range, last success index: {last_success}")
                last_success = -1
                break
    # если оба критария не выполняются, то станем смотреть в сторону менее выгодных вариантов
    else:
        # до тех пор, пока оба критерия не выполняются, ищем менее выгодные варианты
        while not (criteria_fitting and criteria_lowest_fits):
            if half_index >= 0:
                half_index -= 1
                half_df = df[:half_index]
                outer_df = df[half_index:]
                half_max_conc = half_df['conc'].max()
                
                half_amount = half_max_conc * min_volume
                
                delute_volumes = (half_amount / half_df['conc']).round()
                inner_volume = delute_volumes.sum()
                
                outer_volume = outer_df.shape[0] * min_volume
                
                criteria_fitting = (inner_volume + outer_volume) <= max_result_volume
                criteria_lowest_fits = delute_volumes.max() < max_volume
                logger.debug(f"Half index: {half_index}, Criteria fitting: {criteria_fitting}, Criteria lowest fits: {criteria_lowest_fits}, Max delute volume: {delute_volumes.max()}")
                last_success = half_index
            else:
                logger.warning(f"Half index is out of range, last success index: {last_success}")
                last_success = -1
                break

    logger.info(f"Last success index: {last_success}")

    df['water_for_delution'] = 0.0
    df['peptide_for_delution'] = 0.0

    df['deluted_volume'] = 0.0
    df['deluted_concentration'] = 0.0

    if last_success != -1:
        half_index = last_success
    else:
        half_index = 0

    

    half_df = df[:half_index]
    outer_df = df[half_index:]
    half_max_conc = half_df['conc'].max()

    half_amount = half_max_conc * min_volume
    
    logger.success(f"Half index: {half_index}, Half amount: {half_amount}")

    df.loc[outer_df.index, 'water_for_delution'] = min_volume * (outer_df['conc'] - half_max_conc) / half_max_conc
    df.loc[outer_df.index, 'peptide_for_delution'] = min_volume
    df.loc[outer_df.index, 'deluted_volume'] = min_volume
    df.loc[outer_df.index, 'deluted_concentration'] = outer_df['conc'] * min_volume / (min_volume + df.loc[outer_df.index, 'water_for_delution'])

    for idx, row in df.loc[outer_df.index].iterrows():
        if row['water_for_delution'] <= config["additional"]["delution_factors"]["low"]["volume"]:
            df.loc[idx, 'water_for_delution'] = round(row['water_for_delution'] * config["additional"]["delution_factors"]["low"]["factor"])
            df.loc[idx, 'peptide_for_delution'] = round(row['peptide_for_delution'] * config["additional"]["delution_factors"]["low"]["factor"])
        elif row['water_for_delution'] <= config["additional"]["delution_factors"]["medium"]["volume"]:
            df.loc[idx, 'water_for_delution'] = round(row['water_for_delution'] * config["additional"]["delution_factors"]["medium"]["factor"])
            df.loc[idx, 'peptide_for_delution'] = round(row['peptide_for_delution'] * config["additional"]["delution_factors"]["medium"]["factor"])
            df.loc[idx, 'deluted_concentration'] = row['conc'] * row['peptide_for_delution'] / (row['peptide_for_delution'] + row['water_for_delution'])
        elif row['water_for_delution'] <= config["additional"]["delution_factors"]["high"]["volume"]:
            df.loc[idx, 'water_for_delution'] = round(row['water_for_delution'] * config["additional"]["delution_factors"]["high"]["factor"])
            df.loc[idx, 'peptide_for_delution'] = round(row['peptide_for_delution'] * config["additional"]["delution_factors"]["high"]["factor"])
            df.loc[idx, 'deluted_concentration'] = row['conc'] * row['peptide_for_delution'] / (row['peptide_for_delution'] + row['water_for_delution'])
        else:
            df.loc[idx, 'water_for_delution'] = round(row['water_for_delution'])
            df.loc[idx, 'peptide_for_delution'] = round(row['peptide_for_delution'])

    df.loc[half_df.index, 'deluted_concentration'] = half_df['conc']
    df.loc[half_df.index, 'deluted_volume'] = (half_amount / half_df['conc']).round()

if __name__ == "__main__":
    args = arg_parser()
    main(args)
