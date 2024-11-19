"""Read data downloaded from SeSAMe annotations and restructure them to match illuminator data structure"""

from illuminator.annotations import ArrayType, GenomeVersion, GenomeInfo
from illuminator.utils import get_resource_folder, download_from_link, column_names_to_snake_case, concatenate_non_na
from illuminator.utils import get_logger

import pandas as pd

LOGGER = get_logger()

class SesameAnnotations:
    """Extract meaningful information from Sesame data files, and create dataframes with Illuminator format"""

    def __init__(self, array_type: ArrayType, genome_version: GenomeVersion):
        self.array_type = array_type
        self.genome_version = genome_version
        self.mask = self.load_annotation('mask')
        self.manifest = self.load_annotation('manifest')
        self.genome_info = self.load_annotation('genome_info')
        self.gene = self.load_annotation('gene')
        self.probe_infos = self.make_illuminator_probe_info()

    def load_annotation(self, kind: str) -> pd.DataFrame | None:
        """Download or read an annotation file. Kind must be 'mask', 'manifest', 'genome_info' or 'gene'"""

        LOGGER.debug(f'>> loading {kind} for {self.array_type} {self.genome_version} from Sesame')

        # genome info files are handled separately
        if kind == 'genome_info':
            return GenomeInfo('default', self.genome_version)

        # now we can handle mask and manifest files, check that the parameter is not something else
        if kind not in ['mask', 'manifest', 'gene']:
            LOGGER.warning(f'Unknown annotation {kind}, must be one of `mask`, `manifest`, `gene`')
            return None

        # get the annotation resource folder
        data_folder = get_resource_folder('tmp')

        # build the filenames depending on the array type and genome version
        if kind == 'gene':
            # gene files have sub-versions ..
            gene_file_version = ''
            if self.genome_version == GenomeVersion.HG38:
                if self.array_type in [ArrayType.HUMAN_MSA, ArrayType.HUMAN_EPIC_V2]:
                    gene_file_version = 'v41'
                elif self.array_type in [ArrayType.HUMAN_450K, ArrayType.HUMAN_27K, ArrayType.HUMAN_EPIC]:
                    gene_file_version = 'v36'
            elif self.genome_version == GenomeVersion.HG19:
                if self.array_type in [ArrayType.HUMAN_EPIC_V2, ArrayType.HUMAN_450K, ArrayType.HUMAN_27K, ArrayType.HUMAN_EPIC]:
                    gene_file_version = 'v26lift37'
            elif self.genome_version == GenomeVersion.MM10 and self.array_type == ArrayType.MOUSE_MM285:
                gene_file_version = 'vM25'
            if gene_file_version == '':
                LOGGER.warning(f'Gene information : unsupported version {self.genome_version} / {self.array_type}')
                return None
            filename = f'{self.array_type}.{self.genome_version}.manifest.gencode.{gene_file_version}.tsv.gz'
        else:
            filename = f'{self.array_type}.{self.genome_version}.{kind}.tsv.gz'

        # file path specificity for mouse array version MM10
        if self.genome_version == 'mm10':
            filename = f'N296070/{filename}'

        filepath = data_folder.joinpath(filename)

        # if the csv manifest file doesn't exist, download it from sesame repository
        if not filepath.exists():
            print(f'file not found in {filepath}')
            dl_link = f'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/{self.array_type}/{filename}'
            return_status = download_from_link(dl_link, data_folder, filename)
            if return_status == -1:
                return None

        # now read the downloaded manifest file
        df = pd.read_csv(str(filepath), delimiter='\t')

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
            print(f'dropped {ini_size - len(df)} probes with missing illumina ID')
            df = df.astype({'illumina_id': 'int', 'type': 'category', 'probe_type': 'category',
                            'channel' :'category', 'chromosome': 'category', 'start': 'Int64', 'end': 'Int64'})
            df = df.set_index('illumina_id')
            df['probe_type'] = df.probe_type.cat.rename_categories({'rs': 'snp'})  # to improve readability
            if 'strand' not in df.columns:
                print('creating probe strand column')
                df['strand'] = '*'
        else:
            df = df.set_index('probe_id')
            if kind == 'mask':
                df = df.rename(columns={'mask': 'mask_info'})

        LOGGER.info('loading done\n')
        return df

    def make_illuminator_probe_info(self) -> pd.DataFrame | None:
        """Extract useful information from Sesame Manifest, Masks and Genes annotation and merge it in one dataframe
        :return: a pd.DataFrame with IlluminaID as indexes, probes as rows and probes info as columns"""

        if self.manifest is None:
            LOGGER.warning('Make illuminator probe info : provide a manifest first')
            return None

        # select manifest column
        manifest = self.manifest[['probe_id', 'type', 'probe_type', 'channel', 'address_a', 'address_b', 'start', 'end',
                                  'chromosome', 'strand']]

        if self.mask is not None:
            # select mask column (`mask_uniq` or column `mask_info` to get all  the information)
            mask = self.mask[['mask_uniq']].rename(columns={'mask_uniq': 'mask_info'})
            mask.mask_info = mask.mask_info.str.replace(',', ';').replace('"', '')
            manifest = manifest.join(mask, on='probe_id')

        if self.gene is not None:
            # select genes columns
            genes = self.gene[['genes_uniq', 'transcript_types']].rename(columns={'genes_uniq': 'genes'})
            manifest = manifest.join(genes, on='probe_id')
            manifest.transcript_types = manifest.transcript_types.apply(lambda x: ';'.join(set(str(x).replace('nan', '').split(';'))))
        return manifest