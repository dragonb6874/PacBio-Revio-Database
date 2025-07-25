#!/bin/python
import subprocess
import argparse
import sys
import os
import re
import pandas as pd
import xml.etree.ElementTree as ET
import logging
import glob
from functools import lru_cache
from typing import Dict, List, Tuple, Optional, Any, Union

# ==============================================================================
# 로깅 설정
# ==============================================================================
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    filename='xml_parser.log',
    filemode='w',
    datefmt='%m/%d/%Y %I:%M:%S %p'
)
logger = logging.getLogger(__name__)

# ==============================================================================
# 전역 상수 및 설정
# ==============================================================================
NAMESPACES = {
    'pbds': 'http://pacificbiosciences.com/PacBioDatasets.xsd',
    'pbmeta': 'http://pacificbiosciences.com/PacBioCollectionMetadata.xsd',
    'pbsample': 'http://pacificbiosciences.com/PacBioSampleInfo.xsd',
    'pbbase': 'http://pacificbiosciences.com/PacBioBaseDataModel.xsd'
}

DEFAULT_CONFIG = {
    'columns_to_drop': ['RunDetails_StartedBy', 'WellInfo_CreatedAt', 'WellInfo_ModifiedAt']
}

SALESFORCE_EXCEL_PATH = "/BI-NFS/RunLog/PacBio/Revio/Revio log book)-START 20230420.xlsx"
ISSUE_SAMPLE_PATH = "/dlst/raw_data/revio-issue-samples"
DEFAULT_OUTPUT_FILENAME = "revio_report.tsv"
PROCESSED_FILES_SUFFIX = "_processed_files.txt"
FAILED_FILES_SUFFIX = "_failed_files.txt"
ALL_XML_GLOB_PATTERN = "/dlst/tmp-Revio/r*/*/pb_formats/*bc*xml"
EXCLUDE_XML_PATTERN_REGEX = re.compile(r'(unass|fail)')

# ==============================================================================
# 유틸리티 함수
# ==============================================================================
def get_options() -> argparse.Namespace:
    """명령줄 옵션을 파싱합니다."""
    parser = argparse.ArgumentParser(description="Revio Sequencing Report Generation")
    parser.add_argument("-r", "--revioxml", dest="revioxml", metavar="REVIOXML",
                      help="Revio xml file full path", default="None")
    parser.add_argument("-l", "--xmllist", dest="xmllist", metavar="XMLLIST",
                      help="File containing list of XML file paths", default=None)
    parser.add_argument("-o", "--output", dest="output", metavar="OUTPUT",
                      help="Output CSV file name", default=DEFAULT_OUTPUT_FILENAME)
    return parser.parse_args()

def format_number_with_commas(number: Union[int, float, str]) -> str:
    """숫자에 천 단위 구분 쉼표를 추가합니다."""
    try:
        return f"{int(float(str(number).replace(',', ''))):,}"
    except (ValueError, TypeError):
        return str(number)

def safe_extract_element_text(element: Optional[ET.Element], default: str = "") -> str:
    """요소에서 안전하게 텍스트를 추출합니다."""
    return element.text if element is not None and element.text else default

def safe_extract_element_attrib(element: Optional[ET.Element], attr: str, default: str = "") -> str:
    """요소에서 안전하게 속성을 추출합니다."""
    return element.attrib.get(attr, default) if element is not None else default

# ==============================================================================
# XML 파싱 클래스
# ==============================================================================
class XMLParser:
    """XML 파일을 파싱하고 필요한 데이터를 추출하는 클래스"""
    
    def __init__(self):
        self.namespaces = NAMESPACES
    
    @staticmethod
    @lru_cache(maxsize=32)
    def parse_xml(xml_file: str) -> Optional[ET.Element]:
        """XML 파일을 파싱하고 루트 요소를 반환합니다."""
        try:
            return ET.parse(xml_file).getroot()
        except Exception as e:
            logger.error(f"Error parsing XML file {xml_file}: {e}")
            return None
    
    def extract_biosample_name(self, root: ET.Element) -> str:
        """XML 루트 요소에서 BioSample 이름을 추출합니다."""
        # 방법 1: pbsample:BioSample 요소 직접 찾기
        element = root.find('.//pbsample:BioSample', self.namespaces)
        name = safe_extract_element_attrib(element, 'Name')
        if name:
            return name

        # 방법 2: pbmeta:WellSample 내의 pbsample:BioSample 찾기
        wellsample = root.find('.//pbmeta:WellSample', self.namespaces)
        if wellsample is not None:
            element = wellsample.find('.//pbsample:BioSample', self.namespaces)
            name = safe_extract_element_attrib(element, 'Name')
            if name:
                return name

        # 방법 3: ConsensusReadSet의 Name 속성에서 추출
        if root.tag.endswith('ConsensusReadSet'):
            name_attribute = safe_extract_element_attrib(root, 'Name')
            if name_attribute:
                # 괄호 안의 내용 추출
                match = re.search(r'\((.*?)\)', name_attribute)
                if match:
                    return match.group(1)
                # 괄호 앞부분 추출
                parts = name_attribute.split('(')
                if len(parts) > 1:
                    return parts[0].strip()
                return name_attribute

        return "NoName"
    
    def extract_bam_path(self, root: ET.Element, xml_file: str) -> Optional[pd.DataFrame]:
        """XML에서 BAM 파일 경로와 UniqueId를 추출합니다."""
        try:
            xml_dir = os.path.dirname(os.path.abspath(xml_file))
            bam_data = []

            for elem in root.findall('.//pbbase:ExternalResource', self.namespaces):
                resource_id = safe_extract_element_attrib(elem, 'ResourceId')
                unique_id = safe_extract_element_attrib(elem, 'UniqueId')
                
                if resource_id.endswith('.bam') and unique_id:
                    path = os.path.normpath(os.path.join(xml_dir, resource_id)) if not os.path.isabs(resource_id) else resource_id
                    bam_data.append({'Path': path, 'UniqueId': unique_id})

            return pd.DataFrame(bam_data) if bam_data else None
        except Exception as e:
            logger.warning(f"BAM 경로 추출 중 오류 발생: {e}")
            return None

    def extract_cell_statistics(self, xml_file: str) -> Optional[pd.DataFrame]:
        """CellStatistics 정보를 추출합니다."""
        try:
            xml_abs_path = os.path.abspath(xml_file)
            stat_dir = os.path.join(os.path.dirname(os.path.dirname(xml_abs_path)), 'statistics')

            if not (os.path.exists(stat_dir) and os.path.isdir(stat_dir)):
                logger.warning(f"Statistics 디렉토리를 찾을 수 없음: {stat_dir}")
                return None

            report_files = glob.glob(os.path.join(stat_dir, '*ccs_report.txt'))
            if not report_files:
                logger.warning(f"Statistics 디렉토리에 ccs_report.txt 파일을 찾을 수 없음: {stat_dir}")
                return None

            # 최신 파일 선택
            if len(report_files) > 1:
                report_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
                logger.info(f"여러 개의 ccs_report.txt 파일 중 가장 최근 파일 선택: {os.path.basename(report_files[0])}")

            return self._parse_statistics_file(report_files[0])
        except Exception as e:
            logger.error(f"CellStatistics 추출 중 오류 발생: {e}")
            return None

    def _parse_statistics_file(self, report_path: str) -> Optional[pd.DataFrame]:
        """통계 파일을 파싱합니다."""
        try:
            stats_data = {}
            with open(report_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and ':' in line:
                        key, value = [item.strip() for item in line.split(':', 1)]
                        if key and not key.startswith((" ", "-")):
                            stats_data[key.replace(' ', '_')] = value

            return pd.DataFrame([stats_data]) if stats_data else None
        except Exception as e:
            logger.warning(f"보고서 파일 파싱 실패 {report_path}: {e}")
            return None

    def extract_run_details(self, root: ET.Element) -> Optional[pd.DataFrame]:
        """RunDetails 및 CollectionMetadata 정보를 추출합니다."""
        try:
            data = {}

            # ConsensusReadSet 정보
            if root.tag.endswith('ConsensusReadSet'):
                data['WellSampleName'] = safe_extract_element_attrib(root, 'Name')
                data['BioSampleUniqueId'] = safe_extract_element_attrib(root, 'UniqueId')

            # CollectionMetadata 정보
            metadata = root.find('.//pbmeta:CollectionMetadata', self.namespaces)
            if metadata is not None:
                data.update({
                    'InstrumentId': safe_extract_element_attrib(metadata, 'InstrumentId'),
                    'InstrumentName': safe_extract_element_attrib(metadata, 'InstrumentName'),
                    'WellUniqueId': safe_extract_element_attrib(metadata, 'UniqueId')
                })

            # RunDetails 정보
            for run_details in root.findall('.//pbmeta:RunDetails', self.namespaces):
                for child in run_details:
                    tag = child.tag.split('}')[-1]
                    data['RunName' if tag == "Name" else tag] = safe_extract_element_text(child)

            return pd.DataFrame([data]) if data else None
        except Exception as e:
            logger.warning(f"RunDetails 추출 중 오류 발생: {e}")
            return None

    def extract_dataset_metadata(self, root: ET.Element) -> Optional[pd.DataFrame]:
        """DataSetMetadata 정보를 추출합니다."""
        try:
            data = {}

            # TotalLength, NumRecords 추출
            for tag_name in ['TotalLength', 'NumRecords']:
                element = root.find(f'.//pbds:{tag_name}', self.namespaces)
                if element is not None and element.text:
                    data[tag_name] = format_number_with_commas(int(element.text))

            # Context 정보
            metadata = root.find('.//pbmeta:CollectionMetadata', self.namespaces)
            data['Context'] = safe_extract_element_attrib(metadata, 'Context')

            # DNABarcode 정보
            barcode = root.find('.//pbsample:DNABarcode', self.namespaces)
            data['DNABarcodeName'] = safe_extract_element_attrib(barcode, 'Name')

            return pd.DataFrame([data]) if data else None
        except Exception as e:
            logger.warning(f"DataSetMetadata 추출 중 오류 발생: {e}")
            return None

    def extract_wellinfo(self, root: ET.Element) -> Optional[pd.DataFrame]:
        """WellSample 정보를 추출합니다."""
        try:
            data = {}
            well_sample = root.find('.//pbmeta:WellSample', self.namespaces)
            
            if well_sample is not None:
                data.update(well_sample.attrib)
                for child in well_sample:
                    tag = child.tag.split('}')[-1]
                    if tag != 'BioSamples':
                        data[tag] = safe_extract_element_text(child)

            return pd.DataFrame([data]) if data else None
        except Exception as e:
            logger.warning(f"WellSample 정보 추출 중 오류 발생: {e}")
            return None

# ==============================================================================
# 데이터 처리 클래스
# ==============================================================================
class DataProcessor:
    """XML에서 추출된 데이터를 처리하고 정리하는 클래스"""
    
    def __init__(self):
        self.xml_parser = XMLParser()
    
    @staticmethod
    def validate_and_clean_data(data_dict: Dict[str, Any]) -> Dict[str, List[Any]]:
        """데이터 딕셔너리를 검증하고 정리합니다."""
        # 모든 값을 리스트로 변환
        for key, value in list(data_dict.items()):
            if not isinstance(value, list):
                data_dict[key] = [value]

        # 모든 리스트 길이를 통일
        lengths = {len(v) for v in data_dict.values()}
        if len(lengths) > 1:
            max_length = max(lengths)
            for key, value in list(data_dict.items()):
                if len(value) < max_length:
                    data_dict[key] = value * max_length if len(value) == 1 else value + [None] * (max_length - len(value))
        
        return data_dict

    def process_xml_file(self, xml_file: str) -> Tuple[pd.DataFrame, str]:
        """단일 XML 파일을 처리합니다."""
        result_dict = {'XMLFilePath': [xml_file], 'BioSampleName': ["NoName"]}

        root = self.xml_parser.parse_xml(xml_file)
        if not root:
            logger.error(f"XML 파싱 실패: {xml_file}")
            return pd.DataFrame({'Error': ['XML 파싱 실패'], 'XMLFilePath': [xml_file]}), "Error"

        # BioSample 이름 추출
        biosample_name = self.xml_parser.extract_biosample_name(root)
        if biosample_name:
            result_dict['BioSampleName'] = [biosample_name]

        # 데이터 추출
        extractors = [
            ('SampleInfo_', lambda: self.xml_parser.extract_dataset_metadata(root)),
            ('RunDetails_', lambda: self.xml_parser.extract_run_details(root)),
            ('WellInfo_', lambda: self.xml_parser.extract_wellinfo(root)),
            ('CellStat_', lambda: self.xml_parser.extract_cell_statistics(xml_file)),
            ('BamFile_', lambda: self.xml_parser.extract_bam_path(root, xml_file)),
        ]

        for prefix, extract_func in extractors:
            self._process_extraction_result(extract_func(), prefix, result_dict)

        result_dict = self.validate_and_clean_data(result_dict)

        try:
            return pd.DataFrame(result_dict), biosample_name
        except Exception as e:
            logger.error(f"DataFrame 생성 중 오류 ({xml_file}): {e}")
            return pd.DataFrame({'Error': ['처리 실패'], 'XMLFilePath': [xml_file]}), "Error"

    def _process_extraction_result(self, df_extracted: Optional[pd.DataFrame], prefix: str, result_dict: Dict[str, List[Any]]) -> None:
        """추출 결과를 처리합니다."""
        if df_extracted is None or df_extracted.empty or len(df_extracted) < 1:
            return

        for col in df_extracted.columns:
            value = df_extracted[col].values[0]
            if prefix == "BamFile_":
                key = 'BAMPath' if col == 'Path' else f'BamFile{col}'
            else:
                key = f"{prefix}{col}"
            result_dict[key] = [value]

    @staticmethod
    def clean_dataframe(df: pd.DataFrame) -> pd.DataFrame:
        """DataFrame을 정리합니다."""
        # 불필요한 열 제거
        columns_to_drop = [col for col in DEFAULT_CONFIG['columns_to_drop'] if col in df.columns]
        if columns_to_drop:
            df = df.drop(columns_to_drop, axis=1)

        # XMLFilePath 열을 맨 뒤로 이동
        if 'XMLFilePath' in df.columns:
            columns = [col for col in df.columns if col != 'XMLFilePath'] + ['XMLFilePath']
            df = df[columns]

        return df

    @staticmethod
    def extract_salesforce_data() -> pd.DataFrame:
        """세일즈포스 엑셀 파일에서 데이터를 추출합니다."""
        try:
            lab_run_log = pd.read_excel(SALESFORCE_EXCEL_PATH, sheet_name="Revio Demulti")
            lab_run_log = lab_run_log[['Sample ID', 'pooling ID', 'Sample Name', 'run name', '서비스 ID', '서비스명']]
            lab_run_log.columns = ['Sample ID', 'RunDetails_WellSampleName_tmp', 'BioSampleName', 'RunDetails_RunName', 'ServiceId', 'ServiceName']
            
            # 데이터 정리
            lab_run_log['BioSampleName'] = lab_run_log['BioSampleName'].str.strip().replace('', pd.NA).fillna(lab_run_log['Sample ID'])
            
            # 'ServiceId', 'ServiceName' column을 제외한 나머지 열에 대해 ffill 적용
            #lab_run_log.fillna(method='ffill', inplace=True)
            lab_run_log.loc[:, lab_run_log.columns.difference(['ServiceId', 'ServiceName'])] = lab_run_log.loc[:, lab_run_log.columns.difference(['ServiceId', 'ServiceName'])].ffill()

            
            # 문자열 변환 및 공백 제거
            for col in ['BioSampleName', 'Sample ID']:
                lab_run_log[col] = lab_run_log[col].astype(str).str.strip()
            
            lab_run_log['RunDetails_RunName'] = lab_run_log['RunDetails_RunName'].str.replace(' ', '')
            lab_run_log.drop_duplicates(subset=['BioSampleName', 'RunDetails_RunName', 'ServiceId'], keep='last', inplace=True)
            
            return lab_run_log
        except FileNotFoundError:
            logger.error(f"세일즈포스 엑셀 파일을 찾을 수 없습니다: {SALESFORCE_EXCEL_PATH}")
            return pd.DataFrame()
        except Exception as e:
            logger.error(f"세일즈포스 데이터 추출 중 오류 발생: {e}")
            return pd.DataFrame()

    @staticmethod
    def salesforce_ENAME(df: pd.DataFrame) -> pd.DataFrame:
        """ServiceId에 해당하는 EName 정보를 가져옵니다."""
        salesforce_id_list = df['ServiceId'].unique()
        ename_df = pd.DataFrame({"ServiceId": [], "EName": []})
        
        for i in salesforce_id_list:
            ename = subprocess.check_output(
                f"/ess/dlstibm/home/dragonb/workspace/pipeline/PacDB/salesforce -e {i} | grep EName__c | cut -d'\"' -f4", 
                shell=True
            ).strip().decode('utf-8')
            ename_df = pd.concat([ename_df, pd.DataFrame({"ServiceId": [i], "ServiceEName": [ename]})])
            
        return ename_df

    def process_issue_xml(self) -> pd.DataFrame:
        """이슈 샘플 정보를 처리합니다."""
        if not os.path.exists(ISSUE_SAMPLE_PATH):
            logger.warning(f"이슈 샘플 경로를 찾을 수 없음: {ISSUE_SAMPLE_PATH}")
            return pd.DataFrame()

        issue_folders = [d for d in os.listdir(ISSUE_SAMPLE_PATH) if os.path.isdir(os.path.join(ISSUE_SAMPLE_PATH, d))]
        issue_data = []

        logger.info(f"이슈 샘플 폴더 {len(issue_folders)}개 검색 시작...")

        for sample_name in issue_folders:
            issue_data.append(self._process_single_issue_sample(sample_name))

        return pd.DataFrame(issue_data) if issue_data else pd.DataFrame()

    def _process_single_issue_sample(self, sample_name: str) -> Dict[str, Any]:
        """단일 이슈 샘플을 처리합니다."""
        xml_file_path = os.path.join(ISSUE_SAMPLE_PATH, sample_name, f"{sample_name}.consensusreadset.xml")
        
        data = {
            'BioSampleName': sample_name,
            'XMLFilePath': xml_file_path if os.path.exists(xml_file_path) else None,
            'BAMPath': None,
            'SampleInfo_NumRecords': -1,
            'SampleInfo_TotalLength': -1,
            'RunDetails_WhenStarted': None,
            'RunDetails_WhenCreated': None
        }

        if os.path.exists(xml_file_path):
            root = self.xml_parser.parse_xml(xml_file_path)
            if root:
                self._extract_issue_data(root, xml_file_path, data)
        else:
            logger.warning(f"이슈 샘플 XML 파일을 찾을 수 없음: {xml_file_path}")

        return data

    def _extract_issue_data(self, root: ET.Element, xml_file: str, data: Dict[str, Any]) -> None:
        """이슈 샘플 데이터를 추출합니다."""
        # BAM 경로
        bam_df = self.xml_parser.extract_bam_path(root, xml_file)
        if bam_df is not None and not bam_df.empty and 'Path' in bam_df.columns:
            data['BAMPath'] = bam_df['Path'].iloc[0]

        # NumRecords, TotalLength
        dataset_df = self.xml_parser.extract_dataset_metadata(root)
        if dataset_df is not None and not dataset_df.empty:
            for col in ['NumRecords', 'TotalLength']:
                if col in dataset_df.columns:
                    try:
                        data[f'SampleInfo_{col}'] = int(str(dataset_df[col].iloc[0]).replace(',', ''))
                    except ValueError:
                        logger.warning(f"이슈 샘플 {data['BioSampleName']}의 {col} 변환 실패")

        # CreatedAt
        external_resource = root.find('.//pbbase:ExternalResource', self.xml_parser.namespaces)
        if external_resource is not None and 'CreatedAt' in external_resource.attrib:
            created_at = external_resource.attrib.get('CreatedAt')
            data['RunDetails_WhenStarted'] = created_at
            data['RunDetails_WhenCreated'] = created_at
            logger.info(f"샘플 {data['BioSampleName']}의 CreatedAt 추출: {created_at}")

# ==============================================================================
# 리포트 생성 클래스
# ==============================================================================
class ReportGenerator:
    """XML 파싱 결과를 바탕으로 리포트를 생성하는 클래스"""
    
    def __init__(self):
        self.data_processor = DataProcessor()
    
    def process_xml_list(self, xml_files: List[str], output_file: str) -> None:
        """XML 파일 목록을 처리하고 리포트를 생성합니다."""
        # XML 파일 처리
        result_dfs, processed_files, failed_files = self._process_xml_files(xml_files)
        if not result_dfs:
            logger.warning("처리할 수 있는 XML 파일이 없습니다.")
            return

        # 데이터 결합 및 정리
        combined_df = pd.concat(result_dfs, ignore_index=True)
        combined_df = self.data_processor.clean_dataframe(combined_df)

        # 세일즈포스 데이터와 병합
        combined_df = self._merge_salesforce_data(combined_df)

        # 데이터 필터링
        combined_df, rows_to_remove = self._filter_data(combined_df)

        # 이슈 샘플 처리
        combined_df, re_demux_df_final = self._process_issue_samples(combined_df, rows_to_remove)

        # 최종 정리
        combined_df = self._finalize_dataframe(combined_df)

        # 결과 저장
        self._save_results(combined_df, output_file, processed_files, failed_files, rows_to_remove, re_demux_df_final)

    def _process_xml_files(self, xml_files: List[str]) -> Tuple[List[pd.DataFrame], List[str], List[Tuple[str, str]]]:
        """XML 파일들을 처리합니다."""
        result_dfs = []
        processed_files = []
        failed_files = []
        
        logger.info(f"총 {len(xml_files)}개 XML 파일 처리 시작")

        for i, xml_file in enumerate(xml_files):
            try:
                if os.path.exists(xml_file):
                    logger.info(f"파일 처리 중 ({i+1}/{len(xml_files)}): {xml_file}")
                    df_processed, _ = self.data_processor.process_xml_file(xml_file)
                    result_dfs.append(df_processed)
                    processed_files.append(xml_file)
                else:
                    logger.warning(f"파일을 찾을 수 없음 ({i+1}/{len(xml_files)}): {xml_file}")
                    failed_files.append((xml_file, "파일 없음"))
            except Exception as e:
                logger.error(f"파일 처리 중 오류 발생 ({i+1}/{len(xml_files)}): {xml_file}, 오류: {e}")
                failed_files.append((xml_file, str(e)))

        return result_dfs, processed_files, failed_files

    def _merge_salesforce_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """세일즈포스 데이터와 병합합니다."""
        salesforce_data = self.data_processor.extract_salesforce_data()
        if salesforce_data.empty:
            logger.warning("세일즈포스 데이터가 비어있어 병합을 건너뜁니다.")
            return df
        
        logger.info("=" * 62)
        logger.info("1차 세일즈포스 데이터와 병합 시작... (BioSampleName, RunDetails_RunName 기준)")

        # 1차 병합: BioSampleName과 RunDetails_RunName 기준
        if 'RunDetails_RunName' in df.columns:
            df['RunDetails_RunName'] = df['RunDetails_RunName'].dropna().astype(str).str.strip()
            df = pd.merge(df, salesforce_data, on=['BioSampleName', 'RunDetails_RunName'], how='left')

        # 2차 병합: Sample ID와 RunDetails_RunName 기준으로 빈 ServiceId 채우기
        salesforce_data_tmp = salesforce_data.copy()
        salesforce_data_tmp.drop(columns=['BioSampleName'], inplace=True, errors='ignore') 
        salesforce_data_tmp.rename(columns={'Sample ID': 'BioSampleName'}, inplace=True, errors='ignore')
        
        logger.info("=" * 62)
        logger.info("2차 세일즈포스 데이터와 병합 시작... (Sample ID, RunDetails_RunName 기준)")
        if 'ServiceId' in df.columns and 'BioSampleName' in df.columns:
            empty_serviceid_df = df[df['ServiceId'].isna() | (df['ServiceId'] == '')]
            if not empty_serviceid_df.empty:
                logger.info(f"ServiceId가 비어있는 행 {len(empty_serviceid_df)}개 발견. 세일즈포스 데이터로 채우기 시작.")
                
                for index, row in empty_serviceid_df.iterrows():
                    bio_sample_name = row['BioSampleName']
                    run_details_run_name = row.get('RunDetails_RunName', '')

                    matching_row = salesforce_data_tmp[
                        (salesforce_data_tmp['BioSampleName'] == bio_sample_name) &
                        (salesforce_data_tmp['RunDetails_RunName'] == run_details_run_name)
                    ]

                    if not matching_row.empty:
                        df.at[index, 'ServiceId'] = matching_row['ServiceId'].values[0]
                        df.at[index, 'ServiceName'] = matching_row['ServiceName'].values[0]
                        logger.info(f"{index}행 '{bio_sample_name}' -> {matching_row['ServiceId'].values[0]}로 채움.")
                    else:
                        logger.warning(f"{index}행 '{bio_sample_name}' | '{run_details_run_name}' 과 매치되는 세일즈포스 정보 찾을 수 없음.")

        # 3차 병합: BioSampleName과 RunDetails_WellSampleName_tmp 기준
        df['RunDetails_WellSampleName_tmp'] = df['RunDetails_WellSampleName'].str.split('-Cell').str[0]
        
        logger.info("=" * 62)
        logger.info("3차 세일즈포스 데이터와 병합 시작... (BioSampleName, RunDetails_WellSampleName_tmp 기준)")
        if 'ServiceId' in df.columns and 'BioSampleName' in df.columns:
            empty_serviceid_df = df[df['ServiceId'].isna() | (df['ServiceId'] == '')]

            if not empty_serviceid_df.empty:
                logger.info(f"ServiceId가 비어있는 행 {len(empty_serviceid_df)}개 발견. 세일즈포스 데이터로 채우기 시작.")
                
                for index, row in empty_serviceid_df.iterrows():
                    bio_sample_name = row['BioSampleName']
                    run_details_wellsample_name = row.get('RunDetails_WellSampleName_tmp', '')

                    matching_row = salesforce_data_tmp[
                        (salesforce_data_tmp['BioSampleName'] == bio_sample_name) &
                        (salesforce_data_tmp['RunDetails_WellSampleName_tmp'] == run_details_wellsample_name)
                    ]

                    if not matching_row.empty:
                        df.at[index, 'ServiceId'] = matching_row['ServiceId'].values[0]
                        df.at[index, 'ServiceName'] = matching_row['ServiceName'].values[0]
                        logger.info(f"{index}행 '{bio_sample_name}' -> {matching_row['ServiceId'].values[0]}로 채움.")
                    else:
                        logger.warning(f"{index}행 '{bio_sample_name}' | '{run_details_wellsample_name}' 과 매치되는 세일즈포스 정보 찾을 수 없음.")

        # 마지막으로 BioSampleName을 기준으로 ServiceId 채우기
        logger.info("=" * 62)
        logger.info("4차 세일즈포스 데이터와 병합 시작... (BioSampleName 기준)")
        empty_serviceid_df = df[df['ServiceId'].isna() | (df['ServiceId'] == '')]
        for index, value in enumerate(empty_serviceid_df['BioSampleName']):
            idlist = salesforce_data[salesforce_data['BioSampleName']==value]['ServiceId'].unique()
            if len(idlist) == 1:
                df.at[empty_serviceid_df.index[index], 'ServiceId'] = idlist[0]
                df.at[empty_serviceid_df.index[index], 'ServiceName'] = salesforce_data[salesforce_data['BioSampleName']==value]['ServiceName'].values[0]
            elif len(idlist) > 1:
                logger.warning(f"BioSampleName '{value}'에 대해 여러 개의 ServiceId가 발견되었습니다: {idlist}. ServiceId를 채우지 않습니다.")
            else:
                idlist = salesforce_data[salesforce_data['RunDetails_WellSampleName_tmp']==value]['ServiceId'].unique()
                if len(idlist) == 1:
                    df.at[empty_serviceid_df.index[index], 'ServiceId'] = idlist[0]
                    df.at[empty_serviceid_df.index[index], 'ServiceName'] = salesforce_data[salesforce_data['RunDetails_WellSampleName_tmp']==value]['ServiceName'].values[0]
                else:
                    logger.warning(f"BioSampleName '{value}'에 대해 ServiceId를 찾을 수 없습니다. ServiceId를 채우지 않습니다.")

        # 임시 컬럼 제거
        df.drop(columns=['Sample ID', 'RunDetails_WellSampleName_tmp'], inplace=True, errors='ignore')

        return df

    def _filter_data(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """데이터를 필터링합니다."""
        # BioSampleName 필터링
        df = df[df['BioSampleName'].notna()]
        df = df[df['BioSampleName'] != "NoName"]

        # SampleInfo_NumRecords 필터링
        rows_to_remove = pd.DataFrame()
        if 'SampleInfo_NumRecords' in df.columns:
            df['SampleInfo_NumRecords_Numeric'] = df['SampleInfo_NumRecords'].astype(str).str.replace(',', '')
            df['SampleInfo_NumRecords_Numeric'] = pd.to_numeric(df['SampleInfo_NumRecords_Numeric'], errors='coerce').fillna(0).astype(int)

            rows_to_remove = df[df['SampleInfo_NumRecords_Numeric'] < 1000].drop_duplicates(subset=['BioSampleName', 'BAMPath'])
            if not rows_to_remove.empty:
                logger.info("=" * 62)
                logger.info(f"SampleInfo_NumRecords가 1000 미만인 {len(rows_to_remove)}개의 행이 제거 대상입니다.")
                for _, row in rows_to_remove.iterrows():
                    logger.info(f"제거 대상 행: BioSampleName={row.get('BioSampleName', 'N/A')}, SampleInfo_NumRecords={row.get('SampleInfo_NumRecords', 'N/A')}, BAMPath={row.get('BAMPath', 'N/A')}")

            df = df[df['SampleInfo_NumRecords_Numeric'] >= 1000].copy()
            if not df.empty:
                df['SampleInfo_NumRecords'] = df['SampleInfo_NumRecords_Numeric'].apply(format_number_with_commas)
            df.drop(columns=['SampleInfo_NumRecords_Numeric'], inplace=True, errors='ignore')

        return df, rows_to_remove

    def _process_issue_samples(self, df: pd.DataFrame, rows_to_remove: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """이슈 샘플을 처리합니다."""
        re_demux_df_final = pd.DataFrame()
        
        issue_df = self.data_processor.process_issue_xml()
        if issue_df.empty or rows_to_remove.empty:
            return df, re_demux_df_final

        # 제거 대상이었던 행들 중에서, BioSampleName이 issue_df에 있는 행들을 선택
        re_demux_df = rows_to_remove[rows_to_remove['BioSampleName'].isin(issue_df['BioSampleName'])].copy()
        
        # issue_df에서 re_demux_df로 병합되지 않은 나머지 행들
        rename_demux_df = issue_df[issue_df['BioSampleName'].isin(re_demux_df['BioSampleName']) == False].copy()

        if not re_demux_df.empty:
            logger.info(f"이슈 샘플로 판단되어 복원 및 업데이트 대상인 행 {len(re_demux_df)}개 발견.")
            # issue_df에서 XMLFilePath, BAMPath, SampleInfo_NumRecords 정보 가져오기
            issue_info_to_merge = issue_df[['BioSampleName', 'XMLFilePath', 'BAMPath', 'SampleInfo_NumRecords', 'SampleInfo_TotalLength', 'RunDetails_WhenStarted', 'RunDetails_WhenCreated']].copy()
            issue_info_to_merge.rename(columns={
                'XMLFilePath': 'Issue_XMLFilePath',
                'BAMPath': 'Issue_BAMPath',
                'SampleInfo_NumRecords': 'Issue_SampleInfo_NumRecords',
                'SampleInfo_TotalLength': 'Issue_SampleInfo_TotalLength',
                'RunDetails_WhenCreated': 'Issue_RunDetails_WhenCreated',
                'RunDetails_WhenStarted': 'Issue_RunDetails_WhenStarted'
            }, inplace=True)

            re_demux_df = pd.merge(re_demux_df, issue_info_to_merge, on='BioSampleName', how='left')

            # 기존 XMLFilePath, BAMPath, SampleInfo_NumRecords를 이슈 샘플 정보로 업데이트
            re_demux_df['XMLFilePath'] = re_demux_df['Issue_XMLFilePath'].fillna(re_demux_df['XMLFilePath'])
            re_demux_df['BAMPath'] = re_demux_df['Issue_BAMPath'].fillna(re_demux_df['BAMPath'])
            re_demux_df['RunDetails_WhenStarted'] = re_demux_df['Issue_RunDetails_WhenStarted'].fillna(re_demux_df['RunDetails_WhenStarted'])
            re_demux_df['RunDetails_WhenCreated'] = re_demux_df['Issue_RunDetails_WhenCreated'].fillna(re_demux_df['RunDetails_WhenCreated'])
            
            # SampleInfo_NumRecords 업데이트 (쉼표 포맷 적용)
            re_demux_df['SampleInfo_NumRecords'] = re_demux_df['Issue_SampleInfo_NumRecords'].apply(
                lambda x: format_number_with_commas(x) if pd.notnull(x) and x != -1 else "BarcodeIssue"
            )
            # SampleInfo_TotalLength 업데이트 (쉼표 포맷 적용)
            re_demux_df['SampleInfo_TotalLength'] = re_demux_df['Issue_SampleInfo_TotalLength'].apply(
                lambda x: format_number_with_commas(x) if pd.notnull(x) and x != -1 else "BarcodeIssue"
            )
            
            re_demux_df.drop(columns=['Issue_XMLFilePath', 'Issue_BAMPath', 'Issue_SampleInfo_NumRecords', 
                                     'Issue_SampleInfo_TotalLength', 'Issue_RunDetails_WhenStarted', 'Issue_RunDetails_WhenCreated',
                                     'SampleInfo_NumRecords_Numeric'], inplace=True)
            
            # 특정 컬럼 제외하고 "BarcodeIssue"로 채우기
            preserved_cols = ['BioSampleName', 'XMLFilePath', 'BAMPath', 'ServiceId', 'ServiceName', 
                             'SampleInfo_NumRecords', 'SampleInfo_TotalLength', 'RunDetails_WhenStarted', 'RunDetails_WhenCreated']
            for col in re_demux_df.columns:
                if col not in preserved_cols:
                    re_demux_df[col] = "BarcodeIssue"
            
            logger.info(f"재처리된 이슈 샘플 BioSampleName: {re_demux_df['BioSampleName'].unique()}")
            re_demux_df_final = re_demux_df.copy()
            df = pd.concat([df, re_demux_df_final], ignore_index=True)

        # 재생산된 샘플 처리
        if not rename_demux_df.empty:
            logger.info("=" * 62)
            logger.info(f"이슈 샘플이 아니지만 재생산된 샘플 {len(rename_demux_df)}개 발견.")
            rename_demux_df['SampleInfo_NumRecords'] = rename_demux_df['SampleInfo_NumRecords'].apply(
                lambda x: format_number_with_commas(x) if pd.notnull(x) and x != -1 else "BarcodeIssue"
            )
            rename_demux_df['SampleInfo_TotalLength'] = rename_demux_df['SampleInfo_TotalLength'].apply(
                lambda x: format_number_with_commas(x) if pd.notnull(x) and x != -1 else "BarcodeIssue"
            )
            rename_demux_df.drop(columns=['SampleInfo_NumRecords_Numeric'], inplace=True, errors='ignore')
            df = pd.concat([df, rename_demux_df], ignore_index=True)
            logger.info(f"메뉴얼 생산 & 등록 샘플 BioSampleName: {rename_demux_df['BioSampleName'].unique()}")
            logger.info("=" * 62)

        return df, re_demux_df_final

    def _finalize_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """DataFrame을 최종 정리합니다."""
        # ServiceId와 ServiceName 열을 앞쪽으로 이동 (2, 3번째 위치)
        if 'ServiceId' in df.columns and 'ServiceName' in df.columns:
            cols = list(df.columns)
            try:
                cols.remove('ServiceId')
                cols.insert(1, 'ServiceId')
            except ValueError:
                logger.warning("ServiceId 컬럼이 없어 위치를 조정할 수 없습니다.")
            try:
                cols.remove('ServiceName')
                cols.insert(2, 'ServiceName')
            except ValueError:
                logger.warning("ServiceName 컬럼이 없어 위치를 조정할 수 없습니다.")
            df = df[cols]

        # 시간 순 정렬
        if 'RunDetails_WhenCreated' in df.columns:
            df['RunDetails_WhenCreated'] = pd.to_datetime(df['RunDetails_WhenCreated'], errors='coerce')
            df = df.sort_values(by='RunDetails_WhenCreated')

        # 중복 제거 및 번호 부여
        df = df.drop_duplicates(keep='last')
        df.insert(0, 'Number', range(1, len(df) + 1))

        # ServiceEName 추가
        ename_df = self.data_processor.salesforce_ENAME(df)
        df = pd.merge(df, ename_df, on='ServiceId', how='left')
        df.drop(columns=['EName'], inplace=True, errors='ignore')

        # ServiceEName 열을 ServiceName 옆으로 이동
        if 'ServiceEName' in df.columns:
            cols = list(df.columns)
            cols.remove('ServiceEName')
            service_name_index = cols.index('ServiceName') + 1
            cols.insert(service_name_index, 'ServiceEName')
            df = df[cols]
        else:
            logger.warning("ServiceEName 컬럼이 없어 위치를 조정할 수 없습니다.")

        return df

    def _save_results(self, df: pd.DataFrame, output_file: str, processed_files: List[str], 
                     failed_files: List[Tuple[str, str]], rows_to_remove: pd.DataFrame, 
                     re_demux_df_final: pd.DataFrame) -> None:
        """결과를 저장합니다."""
        try:
            df.to_csv(output_file, index=False, sep="\t")
            logger.info(f"통합 데이터가 {output_file}에 저장되었습니다. 총 {len(df)}개의 레코드가 있습니다.")
        except Exception as e:
            logger.error(f"결과 파일 저장 중 오류 발생 ({output_file}): {e}")

        # 처리된 파일 및 실패한 파일 목록 저장
        base_name = os.path.splitext(output_file)[0]
        with open(f"{base_name}{PROCESSED_FILES_SUFFIX}", 'w') as f:
            for xml_file_path in processed_files:
                f.write(f"{xml_file_path}\n")
        logger.info(f"처리된 파일 목록이 {base_name}{PROCESSED_FILES_SUFFIX}에 저장되었습니다.")

        if failed_files:
            with open(f"{base_name}{FAILED_FILES_SUFFIX}", 'w') as f:
                for xml_file_path, reason in failed_files:
                    f.write(f"{xml_file_path}\t{reason}\n")
            logger.warning(f"{len(failed_files)}개 파일 처리 실패, 목록이 {base_name}{FAILED_FILES_SUFFIX}에 저장되었습니다.")

        # 요약 로그
        logger.info("=" * 62)
        logger.info(f"SUMMARY :: TOTAL PROCESSED XMLs                   :: {len(processed_files)}")
        logger.info(f"SUMMARY :: FAILED XMLs                            :: {len(failed_files)}")
        logger.info(f"SUMMARY :: ROWS REMOVED (<1000 reads initially)   :: {len(rows_to_remove)}")
        logger.info(f"SUMMARY :: ROWS RE-DEMUXED (issue samples)        :: {len(re_demux_df_final)}")
        logger.info(f"SUMMARY :: FINAL RECORDS IN OUTPUT                :: {len(df)}")
        logger.info("=" * 62)

    def process_single_xml(self, xml_file: str, output_file: str) -> Optional[str]:
        """단일 XML 파일을 처리합니다."""
        if not os.path.exists(xml_file):
            logger.error(f"입력 XML 파일을 찾을 수 없음: {xml_file}")
            return None

        logger.info(f"단일 XML 파일 처리 중: {xml_file}")
        result_df, biosample_name = self.data_processor.process_xml_file(xml_file)

        if "Error" in biosample_name or result_df.empty:
            logger.error(f"단일 XML 파일 처리 실패: {xml_file}")
            return None

        result_df = self.data_processor.clean_dataframe(result_df)

        # 출력 파일명 결정
        base_name, ext = os.path.splitext(output_file)
        ext = ext or ".tsv"
        
        if biosample_name not in ["NoName", "Error", None, ""]:
            output_path = f"{base_name}_{biosample_name}_report{ext}"
        else:
            output_path = output_file

        try:
            result_df.to_csv(output_path, index=False, sep="\t")
            logger.info(f"결과가 {output_path}에 저장되었습니다.")
            return output_path
        except Exception as e:
            logger.error(f"단일 XML 결과 파일 저장 중 오류 발생 ({output_path}): {e}")
            return None

# ==============================================================================
# 메인 실행 로직
# ==============================================================================
def main() -> int:
    """메인 실행 함수"""
    options = get_options()
    report_generator = ReportGenerator()

    try:
        if options.xmllist and options.xmllist != "None":
            if not os.path.exists(options.xmllist):
                logger.error(f"XML 목록 파일을 찾을 수 없음: {options.xmllist}")
                return 1
            
            with open(options.xmllist, 'r') as f:
                xml_files = [line.strip() for line in f if line.strip().endswith('.xml')]
            
            if not xml_files:
                logger.error(f"XML 목록 파일이 비어 있거나 유효한 XML 파일이 없음: {options.xmllist}")
                return 1
            
            logger.info(f"XML 목록 파일에서 {len(xml_files)}개의 파일을 찾았습니다.")
            report_generator.process_xml_list(xml_files, options.output)

        elif options.revioxml and options.revioxml != "None":
            if not os.path.exists(options.revioxml):
                logger.error(f"XML 파일을 찾을 수 없음: {options.revioxml}")
                return 1
            report_generator.process_single_xml(options.revioxml, options.output)

        else:
            logger.info(f"기본 경로({ALL_XML_GLOB_PATTERN})에서 모든 XML 파일을 처리합니다.")
            all_xmls = glob.glob(ALL_XML_GLOB_PATTERN)
            xml_files = [f for f in all_xmls if not EXCLUDE_XML_PATTERN_REGEX.search(f)]
            
            if not xml_files:
                logger.warning("처리할 XML 파일이 기본 경로에 없습니다.")
            else:
                logger.info(f"기본 경로에서 {len(xml_files)}개의 XML 파일을 처리합니다.")
                report_generator.process_xml_list(xml_files, options.output)

        return 0
    except Exception as e:
        logger.error(f"메인 실행 중 예상치 못한 오류가 발생했습니다: {e}", exc_info=True)
        return 1

if __name__ == "__main__":
    logger.info("=" * 62)
    logger.info("Revio Sequencing Report Generation 스크립트 시작")
    logger.info(f"실행 인수: {' '.join(sys.argv)}")
    logger.info("=" * 62)
    
    exit_code = main()
    
    logger.info("=" * 62)
    if exit_code == 0:
        logger.info("Revio Sequencing Report Generation 스크립트 성공적으로 완료")
    else:
        logger.error(f"Revio Sequencing Report Generation 스크립트 오류로 종료 (Exit Code: {exit_code})")
    logger.info("=" * 62)
    
    sys.exit(exit_code)