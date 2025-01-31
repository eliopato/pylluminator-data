"""Create manifests for every combination of Array Type * Genome version from Illumina or SeSAMe annotation"""

import os.path

from pylluminator.annotations import ArrayType, GenomeVersion
from pylluminator.utils import get_logger

from sesame_annotations import SesameAnnotations
import pandas as pd

LOGGER = get_logger()

root_dir = '_generated_data/annotations'
os.makedirs(root_dir, exist_ok=True)

pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

for at in ArrayType:
    for gv in GenomeVersion:
        if at.is_human() != gv.is_human():
            continue
        LOGGER.info(f'\n------------------- {gv} {at}')
        anno_sesame = SesameAnnotations(at, gv)

        if anno_sesame.manifest is None:
            continue
        current_dir = f'{root_dir}/{gv}/{at}/'
        os.makedirs(current_dir, exist_ok=True)

        if anno_sesame.probe_infos is not None:
            LOGGER.info(f'saving {current_dir}/probe_infos.csv')
            anno_sesame.probe_infos.to_csv(f'{current_dir}/probe_infos.csv')

        if anno_sesame.island_relation is not None:
            LOGGER.info(f'saving {current_dir}/island_relation.csv')
            anno_sesame.island_relation.to_csv(f'{current_dir}/island_relation.csv')
