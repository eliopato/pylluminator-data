"""Brute force script to try to create manifests for every combination of Array Type * Genome version from SeSAMe annotation"""

import os.path
from illuminator.annotations import ArrayType, GenomeVersion
from sesame_annotations import SesameAnnotations

root_dir = 'generated_data/annotations'
os.makedirs(root_dir, exist_ok=True)

for at in ArrayType:
    for gv in GenomeVersion:
        print(f'------------------- {gv} {at}')
        sesame_anno = SesameAnnotations(at, gv)
        if sesame_anno.probe_infos is None:
            continue
        current_dir = f'{root_dir}/{gv}/{at}/'
        os.makedirs(current_dir, exist_ok=True)

        if sesame_anno.probe_infos is not None:
            print(f'saving {current_dir}/probe_infos.csv')
            sesame_anno.probe_infos.to_csv(f'{current_dir}/probe_infos.csv')