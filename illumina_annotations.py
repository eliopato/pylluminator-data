from pylluminator.annotations import ArrayType, GenomeVersion, GenomeInfo
from pylluminator.utils import get_logger
from pylluminator.utils import get_resource_folder, download_from_link, column_names_to_snake_case, concatenate_non_na
from pyliftover import LiftOver

import pandas as pd

from sesame_annotations import SesameAnnotations

LOGGER = get_logger()


def find_tokens_in_file(file_path, tokens) -> dict:
    line_numbers = {}  # dictionary to store line numbers for each token

    with open(file_path, 'r') as file:
        for line_num, line in enumerate(file, start=1):
            for token in tokens:
                if token in line:
                    line_numbers[token] = line_num
                    break

    return line_numbers


class IlluminaAnnotations:
    """Extract meaningful information from Illumina data files, and create dataframes with pylluminator format"""

    def __init__(self, array_type: ArrayType, genome_version: GenomeVersion):
        LOGGER.info('Loading Illumina annotations')
        self.array_type = array_type
        self.genome_version = genome_version
        self.manifest = self.load_manifest()
        if self.manifest is None:
            self.genome_info = None
            self.mask = None
            self.probe_infos = None
            return
        self.genome_info = GenomeInfo('default', genome_version)
        self.mask = self.load_mask()
        self.probe_infos = self.make_pylluminator_probe_info()

    def load_manifest(self) -> pd.DataFrame | None:
        """Download or read a manifest file

        :return: the manifest as a dataframe
        :rtype: pandas.DataFrame | None"""

        # get the annotation resource folder
        data_folder = get_resource_folder('tmp')

        if self.array_type == ArrayType.HUMAN_450K:
            filename = 'humanmethylation450_15017482_v1-2.csv'
        elif self.array_type == ArrayType.HUMAN_MSA:
            filename = 'MSA-48v1-0_20102838_A1.csv'
        elif self.array_type == ArrayType.HUMAN_EPIC:
            filename = 'infinium-methylationepic-v-1-0-b5-manifest-file.csv'
        elif self.array_type == ArrayType.HUMAN_EPIC_V2:
            filename = 'MethylationEPIC v2.0 Files/EPIC-8v2-0_A2.csv'
        elif self.array_type == ArrayType.MOUSE_MM285:
            filename = 'MouseMethylation-12v1-0_A2.csv'
        else:
            LOGGER.warning(f'Illumina annotation : unsupported array type {self.array_type}')
            return None

        filepath = data_folder.joinpath(filename)

        # if the csv manifest file doesn't exist, download it from illumina
        if not filepath.exists():

            if self.array_type == ArrayType.HUMAN_450K:
                # version May 23, 2013
                dl_link = 'https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv'
            elif self.array_type == ArrayType.HUMAN_MSA:
                # version March 4, 2024
                dl_link = 'https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/infiniummethylationscreening/MSA-48v1-0_20102838_A1.csv'
            elif self.array_type == ArrayType.HUMAN_EPIC:
                # version March 13, 2020
                dl_link = 'https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip'
            elif self.array_type == ArrayType.HUMAN_EPIC_V2:
                # version August 16, 2024
                dl_link = 'https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/InfiniumMethylationEPICv2.0ProductFiles(ZIPFormat).zip'
            elif self.array_type == ArrayType.MOUSE_MM285:
                # version May 25, 2021
                dl_link = 'https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/mouse-methylation/infinium-mouse-methylation-manifest-file-csv.zip'
            else:
                LOGGER.warning(f'Unsupported array type {self.array_type}')
                return None

            return_status = download_from_link(dl_link, data_folder, decompress=True)
            if return_status == -1:
                return None

        tokens = ['[Controls]', '[Assay]']
        line_numbers = find_tokens_in_file(filepath, tokens)
        start_row = line_numbers['[Assay]'] if '[Assay]' in line_numbers else 0
        total_nrow = line_numbers['[Controls]'] - start_row - len(line_numbers) if '[Controls]' in line_numbers else None

        # now read the downloaded manifest file
        df = pd.read_csv(str(filepath), delimiter=',', skiprows=start_row, nrows=total_nrow, low_memory=False)

        # uniformization - who likes camel case ?
        df = column_names_to_snake_case(df)
        LOGGER.info(df.columns)

        # get probe ID column - if both ilmn_id and name columns exist, name column is usually without the 4 letters suffix
        df['probe_id'] = df['ilmn_id'] if 'ilmn_id' in df.columns else df['name']

        # extract probe type from probe id (first letters, identifies control probes, snp...)
        df['probe_type'] = df['probe_id'].str.extract(r'^([a-zA-Z]+)')

        df = df.rename(columns={'design_type': 'type', 'infinium_design_type': 'type',
                                'cpg_chrm': 'chromosome', 'chr': 'chromosome',
                                'color_channel': 'channel',
                                'cpg_beg': 'start', 'cpg_end': 'end', 'mapinfo': 'start',
                                'probe_strand': 'strand',
                                'UCSC_RefGene_Name': 'genes',
                                'strand_fr': 'strand', # todo check
                                #'Relation_to_UCSC_CpG_Island': 'island_relation',
                                'address_a_id': 'address_a', 'address_b_id': 'address_b'}, errors='ignore')

        df['start'] = df['start'].astype('Int64')
        df['address_a'] = df['address_a'].astype('Int64')
        df['genome_build'] = df['genome_build'].astype('str')
        if 'end' in df.columns:
            df['end'] = df['end'].astype('Int64')

        # Function to perform the LiftOver conversion
        def convert_position(row):
            result = lo.convert_coordinate(row['chromosome'], row['start'])#, row['strand'])
            if result:
                new_chromosome, new_position = result[0]
                return pd.Series([new_chromosome, new_position])

            return pd.Series([row['chromosome'], row['start']])

        # check if we need to lift over the genome
        # human genome
        print(df.columns)
        if self.genome_version == GenomeVersion.HG19:
            idxs_hg38 = df[df.genome_build.str.contains('38')].index
            if  len(idxs_hg38) > 0:
                print(f'lift over to hg19 for {len(idxs_hg38)} probes')
                lo = LiftOver('hg38', 'hg19')
                # todo : assignment doesnt work here, results in NAs
                df.loc[idxs_hg38, ['chromosome', 'start']] = df.loc[idxs_hg38].apply(convert_position, axis=1)
        elif self.genome_version == GenomeVersion.HG38:
            idxs_hg19 = df[df.genome_build.str.contains('37')].index
            if len(idxs_hg19) > 0:
                print(f'lift over to hg38 for {len(idxs_hg19)} probes')
                lo = LiftOver('hg19', 'hg38')
                df.loc[idxs_hg19, ['chromosome', 'start']] = df.loc[idxs_hg19].apply(convert_position, axis=1)
        # mouse genome
        elif self.genome_version == GenomeVersion.MM10:
            idxs_mm39 = df[df.genome_build.str.contains('mm39')].index
            if len(idxs_mm39) > 0:
                print(f'lift over to mm10 for {len(idxs_mm39)} probes')
                lo = LiftOver('mm39', 'mm10')
                df.loc[idxs_mm39, ['chromosome', 'start']] = df.loc[idxs_mm39].apply(convert_position, axis=1)
        elif self.genome_version == GenomeVersion.MM39:
            idxs_mm10 = df[df.genome_build.str.contains('10')].index
            if len(idxs_mm10) > 0:
                print(f'lift over to mm39 for {len(idxs_mm10)} probes')
                lo = LiftOver('mm10', 'mm39')
                df.loc[idxs_mm10, ['chromosome', 'start']] = df.loc[idxs_mm10].apply(convert_position, axis=1)

        # length 2 for CpG, length 1 for SNP and CpH. beg is 0-based and end is 1-based like in bed files.
        df['end'] = df['start'] + 2
        df.loc[df.probe_type != 'cg', 'end'] -= 1

        # set dataframes index + specific processing for manifest file
        # for type I probes that have both address A and address B set, split them in two rows
        df['illumina_id'] = df.apply(lambda x: concatenate_non_na(x, ['address_a', 'address_b']), axis=1)
        df = df.explode('illumina_id', ignore_index=True)
        df['chromosome'] = df['chromosome'].str.lower().str.replace('chr', '').str.upper()
        # ensure channels are 'R' and 'G', not 'Red' and 'Grn'
        df.loc[~df.channel.isna(), 'channel'] = df.loc[~df.channel.isna()].channel.str[0]
        # turn some columns into categories as it speeds up further processing
        ini_size = len(df)
        df = df.dropna(subset=['illumina_id'])
        if ini_size - len(df) > 0:
            LOGGER.info(f'dropped {ini_size - len(df)} probes with missing illumina ID')
        df = df.astype({'illumina_id': 'int', 'type': 'category', 'probe_type': 'category',
                        'channel' :'category', 'chromosome': 'category'})#, 'start': 'Int64', 'end': 'Int64'}, errors='ignore')
        df = df.set_index('illumina_id')
        df['probe_type'] = df.probe_type.cat.rename_categories({'rs': 'snp'})  # to improve readability
        if 'strand' not in df.columns:
            LOGGER.info('creating probe strand column')
            df['strand'] = '*'

        return df

    def load_mask(self):
        """Load the mask for the corresponding array type and genome version, from SeSAMe data"""
        sesame_annot = SesameAnnotations(self.array_type, self.genome_version, load_all=False)
        return sesame_annot.load_annotation('mask')


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

        return manifest.sort_index()

