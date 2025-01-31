"""Create manifests for every combination of Array Type * Genome version from Illumina or SeSAMe annotation"""

import os.path

from pylluminator.annotations import ArrayType, GenomeVersion
from pylluminator.utils import get_logger

from illumina_annotations import IlluminaAnnotations
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

        # create manifest
        LOGGER.info(f'\n------------------- {gv} {at}')
        anno_illu = IlluminaAnnotations(at, gv)
        anno_sesame = SesameAnnotations(at, gv)

        # compare manifests
        print('--------------------------------------------------------------------')
        if anno_sesame.manifest is not None and anno_illu.manifest is not None:
            if anno_sesame.manifest.equals(anno_illu.manifest):
                print(f'similar for {at} {gv}')
            else:
                print(f'different for {at} {gv}')
                joined = anno_sesame.probe_infos.join(anno_illu.probe_infos, lsuffix='_s', rsuffix='_i')
                print(f'missing {joined['probe_id_s'].isna().sum()} probes in sesame')
                print(f'missing {joined['probe_id_i'].isna().sum()} probes in illumina')
                common_probes = joined.dropna(subset=['probe_id_s', 'probe_id_i'])
                print(f'{len(common_probes)} common probes out of {len(joined)}')
                for col in ['chromosome', 'end', 'start', 'type', 'channel', 'address_a', 'address_b', 'probe_type', 'strand']:
                    common_probes[f'{col}_i'] = common_probes[f'{col}_i'].astype('str')
                    common_probes[f'{col}_s'] = common_probes[f'{col}_s'].astype('str')
                    diff = common_probes[common_probes[f'{col}_i'] != common_probes[f'{col}_s']]
                    if len(diff) > 0:
                        print(f'{len(diff)} different rows in col {col}')
                        print(diff[[f'{col}_i', f'{col}_s']])
        else:
            if anno_sesame.manifest is None:
                print(f'missing sesame manifest for {at} {gv}')
            if anno_illu.manifest is None:
                print(f'missing illumina manifest for {at} {gv}')
