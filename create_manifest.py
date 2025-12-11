"""Create manifests for every combination of Array Type * Genome version from Illumina or SeSAMe annotation"""
import os
from pylluminator.annotations import ArrayType, GenomeVersion
from pylluminator.utils import get_logger

from sesame_annotations import SesameAnnotations
from updated_annotations import UpdatedAnnotations
import pandas as pd

LOGGER = get_logger()

root_dir = '_generated_data/annotations'
gi_root_dir = '_generated_data/genome_info'
os.makedirs(root_dir, exist_ok=True)

pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

for at in ArrayType:
    for gv in GenomeVersion:
        if at.is_human() != gv.is_human():
            continue

        current_dir = f'{root_dir}/{gv}/{at}/'
        os.makedirs(current_dir, exist_ok=True)

        LOGGER.info(f'\n------------------- {gv} {at} - from sesame')
        anno_sesame = SesameAnnotations(at, gv)
        if anno_sesame.manifest is not None and anno_sesame.probe_infos is not None:
            LOGGER.info(f'saving {current_dir}/probe_infos.csv')
            anno_sesame.probe_infos.to_csv(f'{current_dir}/probe_infos.csv')

        if at == ArrayType.HUMAN_EPIC_V2 and gv == GenomeVersion.HG38:
            LOGGER.info(f'\n------------------- {gv} {at} - updated annotations')
            updated_anno = UpdatedAnnotations(at, gv)
            if updated_anno.manifest is not None and updated_anno.probe_infos is not None:
                LOGGER.info(f'saving {current_dir}/probe_infos.updated.csv')
                updated_anno.probe_infos.to_csv(f'{current_dir}/probe_infos.updated.csv')

                # save genome info
                current_gi_dir = f'{gi_root_dir}/{gv}/'
                os.makedirs(gi_root_dir, exist_ok=True)
                os.makedirs(current_gi_dir, exist_ok=True)

                LOGGER.info(f'saving {current_gi_dir}/transcripts_exons.updated.csv')
                updated_anno.genome_info.transcripts_exons.to_csv(f'{current_gi_dir}/transcripts_exons.updated.csv', index=False)

                LOGGER.info(f'saving {current_gi_dir}/transcripts_list.updated.csv')
                updated_anno.genome_info.transcripts_list.to_csv(f'{current_gi_dir}/transcripts_list.updated.csv', index=False)

                LOGGER.info(f'saving {current_gi_dir}/gap_info.updated.csv')
                gi = updated_anno.genome_info.gap_info
                gi.columns = [str.lower(c) for c in gi.columns]
                gi.to_csv(f'{current_gi_dir}/gap_info.updated.csv', index=False)

                LOGGER.info(f'saving {current_gi_dir}/seq_length.updated.csv')
                seq_length = pd.DataFrame(updated_anno.genome_info.seq_length.values(), updated_anno.genome_info.seq_length.keys()).reset_index()
                seq_length.columns=['chromosome', 'seq_length']
                seq_length.to_csv(f'{current_gi_dir}/seq_length.updated.csv', index=False)

                LOGGER.info(f'saving {current_gi_dir}/chromosome_regions.updated.csv')
                updated_anno.genome_info.chromosome_regions.to_csv(f'{current_gi_dir}/chromosome_regions.updated.csv')