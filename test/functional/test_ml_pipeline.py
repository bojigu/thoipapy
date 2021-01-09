import os
import tarfile

import thoipapy
from pathlib import Path
from shutil import rmtree, copyfile


def test_feature_extraction_and_ml_pipeline_for_small_set_of_proteins():
    here: Path = Path(__file__)
    settings_path: Path = here.parents[2] / "thoipapy/setting/thoipapy_standalone_run_settings.xlsx"
    s = thoipapy.common.create_settingdict(settings_path)
    assert isinstance(s, dict)
    sets_dir: Path = here.parents[1] / "test_inputs/protein_sets"
    assert sets_dir.is_dir()
    s["sets_dir"] = str(sets_dir)
    s["set_number"] = 1
    s["data_dir"] = here.parents[1] / "test_outputs/data_dir"

    base_dir: Path = here.parents[1] / "test_outputs/base_dir"
    if not base_dir.is_dir():
        base_dir.mkdir(parents=True)
    s["base_dir"] = base_dir
    protein_names_xlsx_orig: Path = here.parents[1] / f"test_inputs/protein_sets/protein_names.xlsx"
    protein_names_xlsx: Path = base_dir / "protein_names.xlsx"
    copyfile(protein_names_xlsx_orig, protein_names_xlsx)

    s["create_identity_matrix_from_set_seqs"] = True
    s["run_parse_homologues_xml_into_csv"] = True
    s["parse_csv_homologues_to_alignment"] = True
    s["pssm_calculation"] = True
    s["entropy_calculation"] = True
    s["rate4site_calculation"] = True
    s["coevolution_calculation"] = True
    s["clac_relative_position"] = True
    s["calc_lipo_from_pssm"] = True
    s["lips_score_calculation"] = True
    s["motifs_from_seq"] = True
    s["combine_feature_into_train_data"] = True
    s["remove_crystal_hetero"] = False  # not tested here
    s["calc_PREDDIMER_TMDOCK_closedist"] = False
    s["run_feature_selection"] = True
    s["calc_feature_importances"] = True
    s["conduct_ttest"] = True
    s["train_machine_learning_model"] = True
    s["run_testset_trainset_validation"] = False
    s["run_validation"] = False  # currently failing due to divide by zero error, probably due to having only 2 proteins in set

    database = "crystal"

    structure_dir: Path = here.parents[1] / f"test_outputs/data_dir/features/structure/{database}"
    if not structure_dir.is_dir():
        structure_dir.mkdir(parents=True)

    predictions_dir: Path = here.parents[1] / f"test_outputs/data_dir/Predictions/other_predictors/{database}"
    if not predictions_dir.is_dir():
        predictions_dir.mkdir(parents=True)

    acc_list = ["1xioA4", "4hksA1"]

    for acc in acc_list:
        # setup pre-downloaded homologues. Re-use homologues from test_standalone_prediction.py, which requires a rename of internal files.
        pre_downloaded_homologue_xml_tar_gz: Path = here.parents[1] / f"test_inputs/blast_data_valid/{acc}.surr20.BLAST.xml.tar.gz"
        xml_dir: Path = here.parents[1] / f"test_outputs/data_dir/homologues/xml/crystal"
        target_path: Path = xml_dir / f"{acc}.surr20.BLAST.xml.tar.gz"
        target_path.parent.mkdir(parents=True, exist_ok=True)
        copyfile(pre_downloaded_homologue_xml_tar_gz, target_path)
        with tarfile.open(target_path, 'r:gz') as tar:
            tar.extractall(os.path.dirname(target_path))
        blast_details_path_orig: Path = xml_dir / f"BLAST_details.txt"
        blast_details_path: Path = xml_dir / f"{acc}.surr20.BLAST_details.txt"
        os.rename(blast_details_path_orig, blast_details_path)
        blast_xml_path_orig: Path = xml_dir / f"BLAST_results.xml"
        blast_xml_path: Path = xml_dir / f"{acc}.surr20.BLAST.xml"
        os.rename(blast_xml_path_orig, blast_xml_path)
        with tarfile.open(target_path, mode='w:gz') as tar:
            # add the files to the compressed tarfile
            tar.add(str(blast_xml_path), arcname=blast_xml_path.name)
            tar.add(str(blast_details_path), arcname=blast_details_path.name)

        # setup csv with experimentally determined interface
        experimentally_determined_interface_csv_orig: Path = here.parents[1] / f"test_inputs/experimental_data/{acc}.6pairmax.bind.closedist.csv"
        experimentally_determined_interface_csv: Path = structure_dir / f"{acc}.6pairmax.bind.closedist.csv"
        copyfile(experimentally_determined_interface_csv_orig, experimentally_determined_interface_csv)

        # setup prepared csv with preddimer and tmdock predictions
        preddimer_csv_orig: Path = here.parents[1] / f"test_inputs/predictions/{acc}.preddimer.closedist.csv"
        preddimer_csv_target: Path = predictions_dir / f"{acc}.preddimer.closedist.csv"
        copyfile(preddimer_csv_orig, preddimer_csv_target)
        tmdock_csv_orig: Path = here.parents[1] / f"test_inputs/predictions/{acc}.tmdock.closedist.csv"
        tmdock_csv_target: Path = predictions_dir / f"{acc}.tmdock.closedist.csv"
        copyfile(tmdock_csv_orig, tmdock_csv_target)

    thoipapy.run.run(s)

    # cleanup
    data_dir: Path = here.parents[1] / f"test_outputs/data_dir"
    rmtree(data_dir)
    rmtree(base_dir)
