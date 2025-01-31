"""Read data downloaded from SeSAMe annotations and restructure them to match pylluminator data structure"""
import os

from pylluminator.annotations import ArrayType, GenomeVersion, GenomeInfo
from pylluminator.utils import get_resource_folder, download_from_link, column_names_to_snake_case, concatenate_non_na
from pylluminator.utils import get_logger

import pandas as pd

LOGGER = get_logger()

LINKS = {
    'mask': {
        GenomeVersion.HG38: {
            ArrayType.HUMAN_MSA: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/MSA/MSA.hg38.mask.tsv.gz',
            ArrayType.HUMAN_EPIC_V2: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg38.mask.tsv.gz',
            ArrayType.HUMAN_EPIC_PLUS: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC+/EPIC+.hg38.mask.tsv.gz',
            ArrayType.HUMAN_EPIC: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/refs/heads/main/Anno/EPIC/EPIC.hg38.mask.tsv.gz',
            # todo different formats
            # ArrayType.HUMAN_450K: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/archive/202209/HM450.hg38.mask.tsv.gz',
            # ArrayType.HUMAN_27K: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM27/archive/202209/HM27.hg38.mask.tsv.gz',
        },
        GenomeVersion.HG19: {},
        GenomeVersion.MM10: {
            ArrayType.MOUSE_MM285: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/MM285/N296070/MM285.mm10.mask.tsv.gz'
        }
    },
    'manifest': {
        GenomeVersion.HG38: {
            ArrayType.HUMAN_MSA: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/MSA/MSA.hg38.manifest.tsv.gz',
            ArrayType.HUMAN_EPIC_V2: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg38.manifest.tsv.gz',
            ArrayType.HUMAN_EPIC_PLUS: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC+/EPIC+.hg38.manifest.tsv.gz',
            ArrayType.HUMAN_EPIC: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/EPIC.hg38.manifest.tsv.gz',
            ArrayType.HUMAN_450K: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg38.manifest.tsv.gz',
            ArrayType.HUMAN_27K: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM27/HM27.hg38.manifest.tsv.gz',
            ArrayType.MAMMAL_40: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/Mammal40/Mammal40.hg38.manifest.tsv.gz'
        },
        GenomeVersion.HG19: {
            ArrayType.HUMAN_EPIC_V2: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg19.manifest.tsv.gz',
            ArrayType.HUMAN_EPIC: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/EPIC.hg19.manifest.tsv.gz',
            ArrayType.HUMAN_450K: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg19.manifest.tsv.gz',
            ArrayType.HUMAN_27K: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM27/HM27.hg19.manifest.tsv.gz'
        },
        GenomeVersion.MM10: {
            ArrayType.MOUSE_MM285: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/MM285/N296070/MM285.mm10.manifest.tsv.gz'
        },
        GenomeVersion.MM39: {
            ArrayType.MOUSE_MM285: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/MM285/MM285.mm39.manifest.tsv.gz'
        }
    },
    'gene': {
        GenomeVersion.HG38: {
            ArrayType.HUMAN_MSA: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/MSA/MSA.hg38.manifest.gencode.v41.tsv.gz',
            ArrayType.HUMAN_EPIC_V2: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg38.manifest.gencode.v41.tsv.gz',
            ArrayType.HUMAN_EPIC: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/EPIC.hg38.manifest.gencode.v36.tsv.gz',
            ArrayType.HUMAN_450K: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg38.manifest.gencode.v36.tsv.gz',
            ArrayType.HUMAN_27K: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM27/HM27.hg38.manifest.gencode.v36.tsv.gz'
        },
        GenomeVersion.HG19: {
            ArrayType.HUMAN_EPIC_V2: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg19.manifest.gencode.v26lift37.tsv.gz',
            ArrayType.HUMAN_EPIC: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/EPIC.hg19.manifest.gencode.v26lift37.tsv.gz',
            ArrayType.HUMAN_450K: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg19.manifest.gencode.v26lift37.tsv.gz',
            ArrayType.HUMAN_27K: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM27/HM27.hg19.manifest.gencode.v26lift37.tsv.gz'
        },
        GenomeVersion.MM10: {
            ArrayType.MOUSE_MM285: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/MM285/N296070/MM285.mm10.manifest.gencode.vM25.tsv.gz'
        },
    },
    'island_relation': {
        GenomeVersion.HG38: {
            ArrayType.HUMAN_MSA: 'https://github.com/zhou-lab/KYCG_knowledgebase_MSA/raw/refs/heads/main/hg38/CGI.20220904.gz',
            ArrayType.HUMAN_EPIC_V2: 'https://github.com/zhou-lab/KYCG_knowledgebase_EPICv2/raw/refs/heads/main/hg38/CGI.20220904.gz',
            ArrayType.HUMAN_EPIC: 'https://github.com/zhou-lab/KYCG_knowledgebase_EPIC/raw/refs/heads/main/hg38/CGI.20220904.gz',
            ArrayType.HUMAN_450K: 'https://github.com/zhou-lab/KYCG_knowledgebase_HM450/raw/refs/heads/main/hg38/CGI.20220904.gz',
            ArrayType.HUMAN_27K: 'https://github.com/zhou-lab/KYCG_knowledgebase_HM27/raw/refs/heads/main/hg38/CGI.20220904.gz',
            ArrayType.MAMMAL_40: 'https://github.com/zhou-lab/KYCG_knowledgebase_Mammal40/raw/refs/heads/main/hg38/CGI.20220904.gz'
        },
        GenomeVersion.MM10: {
            ArrayType.MOUSE_MM285: 'https://github.com/zhou-lab/KYCG_knowledgebase_MM285/raw/refs/heads/main/mm10/CGI.20220904.gz'
        }
    }
}

class SesameAnnotations:
    """Extract meaningful information from Sesame data files, and create dataframes with pylluminator format"""

    def __init__(self, array_type: ArrayType, genome_version: GenomeVersion, load_all=True):
        if load_all:
            LOGGER.info('Loading SeSAMe annotations')

        self.array_type = array_type
        self.genome_version = genome_version
        if load_all:
            self.mask = self.load_annotation('mask')
            self.manifest = self.load_annotation('manifest')
            self.genome_info = self.load_annotation('genome_info')
            self.gene = self.load_annotation('gene')
            self.island_relation = self.load_annotation('island_relation')
            if self.manifest is not None:
                self.probe_infos = self.make_pylluminator_probe_info()

    def load_annotation(self, kind: str) -> pd.DataFrame | None:
        """Download or read an annotation file. Kind must be 'mask', 'manifest', 'genome_info' or 'gene'"""

        LOGGER.debug(f'>> loading {kind} for {self.array_type} {self.genome_version} from Sesame')

        # genome info files are handled separately
        if kind == 'genome_info':
            return GenomeInfo('default', self.genome_version)

        # now we can handle mask and manifest files, check that the parameter is not something else
        if kind not in LINKS.keys():
            LOGGER.error(f'Unknown annotation {kind}, must be one of {LINKS.keys()}')
            return None

        if self.genome_version not in LINKS[kind].keys():
            LOGGER.error(f'SeSAMe - Unsupported {kind} genome version {self.genome_version}')
            return None

        if self.array_type not in LINKS[kind][self.genome_version]:
            LOGGER.error(f'SeSAMe - Unsupported {kind} array type {self.array_type} for genome version {self.genome_version}')
            return None

        # get the annotation resource folder
        data_folder = get_resource_folder('tmp')
        filelink = LINKS[kind][self.genome_version][self.array_type]
        local_filepath = data_folder.joinpath(filelink.split('/')[-1])

        # all the 'cgi' files are called the same no matter the array type/genome version so we need to delete them each time
        if kind == 'island_relation' and os.path.exists(local_filepath):
            os.remove(local_filepath)

        # if the csv manifest file doesn't exist, download it from sesame repository
        return_status = download_from_link(filelink, data_folder)
        if return_status == -1:
            return None

        # now read the downloaded manifest file
        df = pd.read_csv(str(local_filepath), delimiter='\t')

        # uniformization - who likes camel case ?
        df = column_names_to_snake_case(df)

        # extract probe type from probe id (first letters, identifies control probes, snp...)
        df['probe_type'] = df['probe_id'].str.extract(r'^([a-zA-Z]+)')

        # set dataframes index + specific processing for manifest file
        if kind == 'manifest':
            # for type I probes that have both address A and address B set, split them in two rows
            df['illumina_id'] = df.apply(lambda x: concatenate_non_na(x, ['address_a', 'address_b']), axis=1)
            df = df.explode('illumina_id', ignore_index=True)
            df = df.rename(columns={'design_type': 'type', 'cpg_chrm': 'chromosome',
                                    'cpg_beg': 'start', 'cpg_end': 'end', 'probe_strand': 'strand'}, errors='ignore')
            df['chromosome'] = df['chromosome'].str.lower().str.replace('chr', '').str.upper()
            # turn some columns into categories as it speeds up further processing
            ini_size = len(df)
            df = df.dropna(subset=['illumina_id'])
            LOGGER.info(f'dropped {ini_size - len(df)} probes with missing illumina ID')
            df = df.astype({'illumina_id': 'int', 'type': 'category', 'probe_type': 'category',
                            'channel' :'category', 'chromosome': 'category', 'start': 'Int64', 'end': 'Int64'})
            df = df.set_index('illumina_id')
            df['probe_type'] = df.probe_type.cat.rename_categories({'rs': 'snp'})  # to improve readability
            if 'strand' not in df.columns:
                LOGGER.info('creating probe strand column')
                df['strand'] = '*'
        else:
            df = df.set_index('probe_id')
            if kind == 'mask':
                df = df.rename(columns={'mask': 'mask_info'})

        return df

    def make_pylluminator_probe_info(self) -> pd.DataFrame | None:
        """Extract useful information from Sesame Manifest, Masks and Genes annotation and merge it in one dataframe
        :return: a pd.DataFrame with IlluminaID as indexes, probes as rows and probes info as columns"""

        if self.manifest is None:
            LOGGER.warning('Make pylluminator probe info : provide a manifest first')
            return None

        # select manifest column
        manifest = self.manifest[['probe_id', 'type', 'probe_type', 'channel', 'address_a', 'address_b', 'start', 'end',
                                  'chromosome', 'strand']]

        if self.mask is not None:
            if 'mask_uniq' not in self.mask.columns:
                LOGGER.warning('Make pylluminator probe info : mask missing mask_uniq column, resetting it')
                # todo handle old mask versions
                self.mask = None
            else:
                # select mask column (`mask_uniq` or column `mask_info` to get all  the information)
                mask = self.mask[['mask_uniq']].rename(columns={'mask_uniq': 'mask_info'})
                mask.mask_info = mask.mask_info.str.replace(',', ';').replace('"', '')
                manifest = manifest.join(mask, on='probe_id')

        if self.gene is not None:
            # select genes columns
            genes = self.gene[['genes_uniq']].rename(columns={'genes_uniq': 'genes'})
            manifest = manifest.join(genes, on='probe_id')
            # manifest.transcript_types = manifest.transcript_types.apply(lambda x: ';'.join(set(str(x).replace('nan', '').split(';'))))

        if self.island_relation is not None:
            self.island_relation['rel'] = self.island_relation['knowledgebase'].str.replace('CGI;', '')
            self.island_relation = self.island_relation['rel']  # keep only this column
            # print(self.island_relation)

        return manifest.sort_index()