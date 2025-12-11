"""Updated EPICv2 annotation as defined by https://www.biorxiv.org/content/10.1101/2025.03.12.642895v2 """
import posixpath

from pylluminator.annotations import ArrayType, GenomeVersion, GenomeInfo
from pylluminator.utils import get_resource_folder, download_from_link, column_names_to_snake_case, concatenate_non_na
from pylluminator.utils import get_logger

import pandas as pd

LOGGER = get_logger(level=0)

LINKS = {
    'mask': {
        GenomeVersion.HG38: {
            ArrayType.HUMAN_EPIC_V2: 'https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg38.mask.tsv.gz',
        }
    },
    'manifest': {
        GenomeVersion.HG38: {
            ArrayType.HUMAN_EPIC_V2: 'https://zenodo.org/records/14933469/files/EPICv2_reannotated_manifest_v1.0.csv.zip'
        }
    },
    #genome_info also needs the manifest file
    'genome_info': {
        GenomeVersion.HG38: {
            ArrayType.HUMAN_EPIC_V2: 'https://zenodo.org/records/14933469/files/EPICv2_reannotated_manifest_v1.0.csv.zip'
        }
    }
}

class UpdatedAnnotations:
    """Extract meaningful information from annotation files, and create dataframes with pylluminator format"""

    def __init__(self, array_type: ArrayType, genome_version: GenomeVersion, load_all=True):
        if load_all:
            LOGGER.info('Loading updated annotations')

        self.array_type = array_type
        self.genome_version = genome_version
        if load_all:
            self.mask = self.load_annotation('mask')
            self.manifest = self.load_annotation('manifest')
            self.genome_info = self.load_annotation('genome_info')
            if self.manifest is not None:
                self.probe_infos = self.make_pylluminator_probe_info()

    def load_annotation(self, kind: str) -> pd.DataFrame | None:
        """Download or read an annotation file. Kind must be 'mask', 'manifest', 'genome_info' or 'gene'"""

        LOGGER.debug(f'>> loading {kind} for {self.array_type} {self.genome_version} from updated annotation')

        # now we can handle mask and manifest files, check that the parameter is not something else
        if kind not in LINKS.keys():
            LOGGER.error(f'Unknown annotation {kind}, must be one of {LINKS.keys()}')
            return None

        if self.genome_version not in LINKS[kind].keys():
            LOGGER.error(f'Updated annotations - Unsupported {kind} genome version {self.genome_version}')
            return None

        if self.array_type not in LINKS[kind][self.genome_version]:
            LOGGER.error(f'Updated annotations - Unsupported {kind} array type {self.array_type} for genome version {self.genome_version}')
            return None

        # get the annotation resource folder
        data_folder = get_resource_folder('tmp')
        filelink = LINKS[kind][self.genome_version][self.array_type]
        local_filepath = data_folder.joinpath(filelink.split('/')[-1])

        # if the csv manifest file doesn't exist, download it from sesame repository
        
        if kind !=  'mask':
            # uncompress zip (specific code bc it's a folder with several csv)
            return_status = download_from_link(filelink, data_folder, decompress=True)
            local_filepath = str(local_filepath)[:-4]
            delimiter = ','
        else:
            return_status = download_from_link(filelink, data_folder)
            delimiter = '\t'

        if return_status == -1:
            return None

        # now read the downloaded manifest file
        LOGGER.info(f'read file at {local_filepath}')
        df = pd.read_csv(str(local_filepath), delimiter=delimiter, low_memory=False)

        # genome info files are handled separately
        if kind == 'genome_info':
            
            genome_info = GenomeInfo('illumina', self.genome_version)

            transcript_exon = df[['CHR', 'GENCODEv47_Transcript_Start', 'GENCODEv47_Transcript_End', 'GENCODEv47_Gene_Strand', 'GENCODEv47_Transcript_ID', 'GENCODEv47_Transcript_Name', 
                                   'GENCODEv47_Gene_Name', 'GENCODEv47_Gene_ID', 'GENCODEv47_Gene_Type',  'GENCODEv47_Gene_Start',  'GENCODEv47_Gene_End',
                                   'GENCODEv47_Feature_Type', 'GENCODEv47_Feature_Exon_Number', 'GENCODEv47_Feature_Start', 'GENCODEv47_Feature_End']].dropna(subset='GENCODEv47_Transcript_ID').drop_duplicates()
            transcript_exon = transcript_exon.map(lambda x: x.split(';'))
            transcript_exon = transcript_exon.explode(list(transcript_exon.columns[1:]))
            transcript_exon.columns = [str.lower(c.replace('GENCODEv47_', '')) for c in transcript_exon.columns]
            transcript_exon.chr = transcript_exon.chr.str[0].str.replace('chr', '')
            transcript_exon = transcript_exon.astype({'transcript_start': int, 'transcript_end': int, 'gene_start': int,  'feature_start': int, 'feature_end': int}) #'gene_end': 'Int64',
            transcript_exon = transcript_exon.drop_duplicates()
            transcript_exon = transcript_exon.rename(columns={'chr': 'chromosome'})

            transcript_list = transcript_exon[['transcript_id', 'feature_start', 'feature_end', 'feature_exon_number']].drop_duplicates()
            transcript_list = transcript_list.groupby(['transcript_id', 'feature_exon_number']).agg({'feature_start': 'min', 'feature_end': 'max'}).reset_index()
            transcript_list.columns = ['group_name', 'start', 'end', 'exon_number']

            transcript_exon = transcript_exon.drop(columns=['feature_start', 'feature_end', 'feature_exon_number']).drop_duplicates()

            genome_info.transcripts_exons = transcript_exon
            genome_info.transcripts_list = transcript_list

            return genome_info

        # set dataframes index + specific processing for manifest file
        if kind == 'manifest':
            df = df[['IlmnID', 'Name', 'Infinium_Design_Type', 'Probe_Type', 'col', 'AddressA_ID', 'AddressB_ID','MAPINFO', 'CHR', 'Strand_FR', 'GENCODEv47_Gene_Name', 'Regulatory_Feature_Group', 'Relation_to_UCSC_CpG_Island']].copy()
            df.columns = ['probe_id', 'name', 'type', 'probe_type', 'channel', 'address_a', 'address_b', 'start', 'chromosome', 'strand', 'genes', 'promoter_or_body', 'cgi']
            # clean info of control probes
            df.probe_type = df.probe_type.fillna('ctl')  # NA types are control probes
            ctrl_probes_idxs = df.probe_type == 'ctl'
            df.loc[ctrl_probes_idxs, 'address_a'] = df.loc[df.probe_type == 'ctl', 'probe_id']  # set control probes to type II 
            df.loc[ctrl_probes_idxs, 'type'] = df.loc[ctrl_probes_idxs, 'type'].fillna('II')  # set control probes to type II 
            df.loc[ctrl_probes_idxs, 'probe_id'] = 'ctl_' + df.loc[ctrl_probes_idxs, 'address_a'] + '_' +df.loc[ctrl_probes_idxs, 'name']   # set control probes to type II 
            df.loc[ctrl_probes_idxs, 'strand'] = df.loc[ctrl_probes_idxs, 'strand'].fillna('*')
            df = df.drop(columns='name')
            # clean address b
            df.loc[(df.address_b == 'AVG'), 'address_b'] = None
            # add illumina ID
            df['illumina_id'] = df['address_a']
            # add address_b illumina ID
            add_b = df[~pd.isna(df.address_b)].copy()
            add_b.illumina_id = add_b.address_b
            df = pd.concat([df.set_index('illumina_id'), add_b.set_index('illumina_id')])
            # uniformize promoter
            # df.loc[~pd.isna(df.promoter_or_body) & df.promoter_or_body.str.startswith('Promoter'), 'promoter_or_body'] = 'p'
            # fix start and end
            df.loc[~pd.isna(df.start), 'start'] = df.loc[~pd.isna(df.start), 'start'] - 1
            df.loc[~pd.isna(df.start), 'end'] = df.loc[~pd.isna(df.start), 'start'] + 2
            # remove redundant "chr" from chromosome values
            df.loc[~pd.isna(df.chromosome), 'chromosome'] = df.loc[~pd.isna(df.chromosome), 'chromosome'].str.replace('chr', '')
            # keep only unique genes in genes list
            df.loc[~pd.isna(df.genes), 'genes'] = df.loc[~pd.isna(df.genes), 'genes'].str.split(';').map(lambda x: ';'.join(set(x)))
            df.loc[~pd.isna(df.cgi), 'cgi'] = df.loc[~pd.isna(df.cgi), 'cgi'].str.split(';').map(lambda x: ';'.join(set(x)))
            # rename snp category
            df.loc[df.probe_type == 'rs', 'probe_type'] = 'snp'
            # remove any NA probe ID (like "Controls" header)
            df = df.dropna(subset='probe_id')
        else:
            df = column_names_to_snake_case(df).set_index('probe_id')
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
        manifest = self.manifest

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

        return manifest.sort_index().sort_values('probe_id')