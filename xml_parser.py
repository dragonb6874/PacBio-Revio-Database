#!/bin/python
import subprocess
import argparse # optparse 대신 argparse 임포트
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

# XML 네임스페이스 정의
NAMESPACES = {
    'pbds': 'http://pacificbiosciences.com/PacBioDatasets.xsd',
    'pbmeta': 'http://pacificbiosciences.com/PacBioCollectionMetadata.xsd',
    'pbsample': 'http://pacificbiosciences.com/PacBioSampleInfo.xsd',
    'pbbase': 'http://pacificbiosciences.com/PacBioBaseDataModel.xsd'
}

# 기본 설정
DEFAULT_CONFIG = {
    'columns_to_drop': ['RunDetails_StartedBy', 'WellInfo_CreatedAt', 'WellInfo_ModifiedAt']
}

# 파일 경로 및 패턴 정의
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
        """
        XML 파일을 파싱하고 루트 요소를 반환합니다.
        lru_cache 데코레이터를 사용하여 최근 파싱 결과를 캐싱하여 성능을 개선합니다.
        """
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            return root
        except Exception as e:
            logger.error(f"Error parsing XML file {xml_file}: {e}")
            return None
    
    def extract_biosample_name(self, root: ET.Element) -> str:
        """XML 루트 요소에서 BioSample 이름을 추출합니다. 여러 위치에서 시도합니다."""
        try:
            # 방법 1: pbsample:BioSample 요소 직접 찾기
            bio_sample_element = root.find('.//pbsample:BioSample', self.namespaces)
            if bio_sample_element is not None and 'Name' in bio_sample_element.attrib:
                return bio_sample_element.attrib.get('Name')

            # 방법 2: pbmeta:WellSample 내의 pbsample:BioSample 찾기
            wellsample_element = root.find('.//pbmeta:WellSample', self.namespaces)
            if wellsample_element is not None:
                bio_samples_in_wellsample = wellsample_element.find('.//pbsample:BioSample', self.namespaces)
                if bio_samples_in_wellsample is not None and 'Name' in bio_samples_in_wellsample.attrib:
                    return bio_samples_in_wellsample.attrib.get('Name')

            # 방법 3: ConsensusReadSet의 Name 속성에서 추출 (괄호 안의 내용 또는 괄호 앞부분)
            if root.tag.endswith('ConsensusReadSet') and 'Name' in root.attrib:
                name_attribute = root.attrib.get('Name')
                match = re.search(r'\((.*?)\)', name_attribute)  # 괄호 안의 내용
                if match:
                    return match.group(1)
                parts = name_attribute.split('(')  # 괄호 기준으로 분리
                if len(parts) > 1:
                    return parts[0].strip()  # 괄호 앞부분
                return name_attribute  # 괄호가 없는 경우 전체 이름
        except Exception as e:
            logger.warning(f"BioSample 이름 추출 중 오류 발생: {e}")

        return "NoName"  # 모든 방법 실패 시 기본값
    
    def extract_bam_path(self, root: ET.Element, xml_file: str) -> Optional[pd.DataFrame]:
        """
        XML 루트 요소에서 BAM 파일 경로와 UniqueId를 추출합니다.
        BAM 파일 경로는 XML 파일 위치를 기준으로 상대 경로이거나 절대 경로일 수 있습니다.
        """
        try:
            xml_dir = os.path.dirname(os.path.abspath(xml_file))
            bam_data_list = []

            for elem in root.findall('.//pbbase:ExternalResource', self.namespaces):
                if 'ResourceId' in elem.attrib and 'UniqueId' in elem.attrib:
                    resource_id = elem.attrib.get('ResourceId')
                    if resource_id.endswith('.bam'):
                        # 경로 정규화 (상대 경로인 경우 절대 경로로 변환)
                        path = os.path.normpath(os.path.join(xml_dir, resource_id)) if not os.path.isabs(resource_id) else resource_id
                        bam_data_list.append({
                            'Path': path,
                            'UniqueId': elem.attrib.get('UniqueId')
                        })

            return pd.DataFrame(bam_data_list) if bam_data_list else None
        except Exception as e:
            logger.warning(f"BAM 경로 추출 중 오류 발생: {e}")
            return None

    def extract_cell_statistics(self, xml_file: str) -> Optional[pd.DataFrame]:
        """
        CellStatistics 정보를 담고 있는 ccs_report.txt 파일에서 데이터를 추출합니다.
        XML 파일의 상위 디렉토리 내 'statistics' 폴더에서 해당 파일을 찾습니다.
        """
        try:
            xml_abs_path = os.path.abspath(xml_file)
            parent_dir = os.path.dirname(os.path.dirname(xml_abs_path))  # XML 파일의 상위 디렉토리
            stat_dir = os.path.join(parent_dir, 'statistics')

            if os.path.exists(stat_dir) and os.path.isdir(stat_dir):
                report_files = glob.glob(os.path.join(stat_dir, '*ccs_report.txt'))

                if not report_files:
                    logger.warning(f"Statistics 디렉토리에 ccs_report.txt 파일을 찾을 수 없음: {stat_dir}")
                    return None

                # 여러 파일이 있을 경우 가장 최근 수정된 파일 선택
                if len(report_files) > 1:
                    report_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
                    logger.info(f"여러 개의 ccs_report.txt 파일 중 가장 최근 파일 선택: {os.path.basename(report_files[0])}")

                report_path = report_files[0]
                stats_data = {}
                with open(report_path, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line and ':' in line:
                            key, value = [item.strip() for item in line.split(':', 1)]
                            # 유효한 키만 추출 (공백이나 '-'로 시작하지 않는 키)
                            if key and not key.startswith(" ") and not key.startswith("-"):
                                stats_data[f"{key.replace(' ', '_')}"] = value  # 공백을 '_'로 변경

                if stats_data:
                    return pd.DataFrame([stats_data])
                else:
                    logger.warning(f"보고서 파일에서 데이터를 추출할 수 없음: {report_path}")
                    return None
            else:
                logger.warning(f"Statistics 디렉토리를 찾을 수 없음: {stat_dir}")
                return None
        except Exception as e:
            logger.error(f"CellStatistics 추출 중 오류 발생: {e}")
            return None

    def extract_run_details(self, root: ET.Element) -> Optional[pd.DataFrame]:
        """XML 루트 요소에서 RunDetails 및 CollectionMetadata 관련 정보를 추출합니다."""
        try:
            data = {}

            # ConsensusReadSet 정보 (Name, UniqueId)
            if root.tag.endswith('ConsensusReadSet'):
                for attr, key_name in [('Name', 'WellSampleName'), ('UniqueId', 'BioSampleUniqueId')]:
                    if attr in root.attrib:
                        data[key_name] = root.attrib.get(attr)

            # CollectionMetadata 정보 (InstrumentId, InstrumentName, UniqueId)
            metadata_element = root.find('.//pbmeta:CollectionMetadata', self.namespaces)
            if metadata_element is not None:
                for attr, key_name in [('InstrumentId', 'InstrumentId'), ('InstrumentName', 'InstrumentName'),
                                     ('UniqueId', 'WellUniqueId')]:
                    if attr in metadata_element.attrib:
                        data[key_name] = metadata_element.attrib.get(attr)

            # RunDetails 정보 (하위 모든 요소)
            for run_details_elem in root.findall('.//pbmeta:RunDetails', self.namespaces):
                for child in run_details_elem:
                    tag = child.tag.split('}')[-1]  # 네임스페이스 제거
                    data['RunName' if tag == "Name" else tag] = child.text  # 'Name' 태그는 'RunName'으로 변경

            return pd.DataFrame([data]) if data else None
        except Exception as e:
            logger.warning(f"RunDetails 추출 중 오류 발생: {e}")
            return None

    def extract_dataset_metadata(self, root: ET.Element) -> Optional[pd.DataFrame]:
        """XML 루트 요소에서 DataSetMetadata, Context, DNABarcode 정보를 추출합니다."""
        try:
            data = {}

            # TotalLength, NumRecords 추출 (숫자는 쉼표 포맷 적용)
            for tag_name, key_name in [('TotalLength', 'TotalLength'), ('NumRecords', 'NumRecords')]:
                element = root.find(f'.//pbds:{tag_name}', self.namespaces)
                if element is not None and element.text:
                    data[key_name] = format_number_with_commas(int(element.text))

            # Context 정보 추출
            metadata_element = root.find('.//pbmeta:CollectionMetadata', self.namespaces)
            if metadata_element is not None and 'Context' in metadata_element.attrib:
                data['Context'] = metadata_element.attrib.get('Context')

            # DNABarcode 이름 추출
            barcode_element = root.find('.//pbsample:DNABarcode', self.namespaces)
            if barcode_element is not None and 'Name' in barcode_element.attrib:
                data['DNABarcodeName'] = barcode_element.attrib.get('Name')

            return pd.DataFrame([data]) if data else None
        except Exception as e:
            logger.warning(f"DataSetMetadata 추출 중 오류 발생: {e}")
            return None

    def extract_wellinfo(self, root: ET.Element) -> Optional[pd.DataFrame]:
        """XML 루트 요소에서 WellSample 관련 정보를 추출합니다."""
        try:
            data = {}
            well_sample_element = root.find('.//pbmeta:WellSample', self.namespaces)
            if well_sample_element is not None:
                # WellSample 요소의 모든 속성 추출
                for attr, value in well_sample_element.attrib.items():
                    data[attr] = value

                # WellSample의 하위 요소 추출 (BioSamples 제외)
                for child in well_sample_element:
                    tag = child.tag.split('}')[-1]  # 네임스페이스 제거
                    if tag != 'BioSamples':
                        data[tag] = child.text

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
        """
        데이터 딕셔너리를 검증하고 정리합니다.
        DataFrame 생성을 위해 모든 값 리스트의 길이를 통일시키는 역할을 합니다.
        """
        # 모든 값이 리스트인지 확인하고 아니면 리스트로 변환
        items = list(data_dict.items())
        for key, value in items:
            if not isinstance(value, list):
                data_dict[key] = [value]

        # DataFrame 생성을 위해 모든 리스트의 길이가 같은지 확인
        lengths = {len(v) for v in data_dict.values()}
        if len(lengths) > 1:
            max_length = max(lengths)
            items_for_length_check = list(data_dict.items())
            for key, value in items_for_length_check:
                if len(value) < max_length:
                    # 리스트 길이를 맞추기: 단일 값이면 반복, 아니면 None으로 채움
                    if len(value) == 1:
                        data_dict[key] = value * max_length
                    else:
                        data_dict[key] = value + [None] * (max_length - len(value))
        return data_dict

    def process_xml_file(self, xml_file: str) -> Tuple[pd.DataFrame, str]:
        """
        단일 XML 파일을 처리하여 주요 정보를 추출하고 DataFrame으로 반환합니다.
        """
        result_dict = {'XMLFilePath': [xml_file], 'BioSampleName': ["NoName"]}

        # XML 파싱
        root = self.xml_parser.parse_xml(xml_file)
        if not root:
            logger.error(f"XML 파싱 실패: {xml_file}")
            return pd.DataFrame({'Error': ['XML 파싱 실패'], 'XMLFilePath': [xml_file]}), "Error"

        # BioSample 이름 추출
        biosample_name = self.xml_parser.extract_biosample_name(root)
        if biosample_name:
            result_dict['BioSampleName'] = [biosample_name]

        # 각 데이터 추출 함수와 결과에 사용할 접두사 정의
        extractors = [
            (lambda: self.xml_parser.extract_dataset_metadata(root), "SampleInfo_"),
            (lambda: self.xml_parser.extract_run_details(root), "RunDetails_"),
            (lambda: self.xml_parser.extract_wellinfo(root), "WellInfo_"),
            (lambda: self.xml_parser.extract_cell_statistics(xml_file), "CellStat_"),
            (lambda: self.xml_parser.extract_bam_path(root, xml_file), "BamFile_"),
        ]

        # 정의된 추출 함수들을 순차적으로 실행하여 데이터 수집
        for extract_func, prefix in extractors:
            try:
                df_extracted = extract_func()
                if df_extracted is not None and not df_extracted.empty:
                    for col in df_extracted.columns:
                        if len(df_extracted) >= 1:  # 첫 번째 행 데이터만 사용
                            key = f"{prefix}{col}"
                            # BamFile 관련 컬럼은 특별 처리 (BAMPath, BamFileUniqueId)
                            if prefix == "BamFile_" and col == "Path":
                                result_dict['BAMPath'] = [df_extracted[col].values[0]]
                            elif prefix == "BamFile_" and col == "UniqueId":
                                result_dict['BamFileUniqueId'] = [df_extracted[col].values[0]]
                            else:
                                result_dict[key] = [df_extracted[col].values[0]]
            except Exception as e:
                logger.warning(f"{prefix} 데이터 추출 중 오류 발생 ({xml_file}): {e}")

        # 데이터 정리 및 유효성 검사 (길이 맞추기 등)
        result_dict = self.validate_and_clean_data(result_dict)

        # 최종 DataFrame 생성
        try:
            return pd.DataFrame(result_dict), biosample_name
        except Exception as e:
            logger.error(f"DataFrame 생성 중 오류 ({xml_file}): {e}")
            return pd.DataFrame({'Error': ['처리 실패'], 'XMLFilePath': [xml_file]}), "Error"

    @staticmethod
    def clean_dataframe(df: pd.DataFrame) -> pd.DataFrame:
        """
        생성된 DataFrame에서 불필요한 열을 제거하고 열 순서를 조정합니다.
        """
        # DEFAULT_CONFIG에 정의된 불필요한 열 제거
        columns_to_drop = [col for col in DEFAULT_CONFIG['columns_to_drop'] if col in df.columns]
        if columns_to_drop:
            df = df.drop(columns_to_drop, axis=1)

        # XMLFilePath 열을 맨 뒤로 이동
        if 'XMLFilePath' in df.columns:
            columns_ordered = [col for col in df.columns if col != 'XMLFilePath']
            columns_ordered.append('XMLFilePath')
            df = df[columns_ordered]

        return df

    @staticmethod
    def extract_salesforce_data() -> pd.DataFrame:
        """
        세일즈포스 엑셀 파일에서 랩 실행 로그 데이터를 읽어와 전처리합니다.
        BioSampleName, RunName, ServiceId, ServiceName을 추출하고 정제합니다.
        """
        try:
            # Sample Name : runlog demulti 시트 M 열
            # pooling ID : runlog demulti 시트 B 열
            # Sample ID : runlog demulti 시트 E 열
            # run name : runlog demulti 시트 J 열
            lab_run_log = pd.read_excel(SALESFORCE_EXCEL_PATH, sheet_name="Revio Demulti")
            lab_run_log = lab_run_log[['Sample ID', 'pooling ID', 'Sample Name', 'run name', '서비스 ID', '서비스명']]
            lab_run_log.columns = ['Sample ID', 'RunDetails_WellSampleName_tmp', 'BioSampleName', 'RunDetails_RunName', 'ServiceId', 'ServiceName']
            
            # BioSampleName 양쪽 공백 제거 후 비어 있다면 Sample Name column 값으로 채움
            lab_run_log['BioSampleName'] = lab_run_log['BioSampleName'].str.strip()
            lab_run_log['BioSampleName'] = lab_run_log['BioSampleName'].replace('', pd.NA)
            lab_run_log['BioSampleName'] = lab_run_log['BioSampleName'].fillna(lab_run_log['Sample ID'])
            
            # Sample ID 열 제거
            #lab_run_log.drop(columns=['Sample ID'], inplace=True, errors='ignore')
            
            # BioSampleName이 비어있는 행 제거
            #lab_run_log = lab_run_log[lab_run_log['BioSampleName'].notnull()]
            #if 'BioSampleName' in lab_run_log.columns:
            #    lab_run_log = lab_run_log[lab_run_log['BioSampleName'] != '']

            # 병합된 셀로 인해 비어있는 RunName을 이전 값으로 채움 (forward fill)
            lab_run_log.fillna(method='ffill', inplace=True)

            # BioSampleName을 문자열로 변환
            lab_run_log['BioSampleName'] = lab_run_log['BioSampleName'].astype(str)
            lab_run_log['Sample ID'] = lab_run_log['Sample ID'].astype(str)

            # BioSampleName의 양쪽 공백 제거
            lab_run_log['BioSampleName'] = lab_run_log['BioSampleName'].str.strip()
            lab_run_log['Sample ID'] = lab_run_log['Sample ID'].str.strip()

            # RunDetails_RunName 열에서 공백 제거
            lab_run_log['RunDetails_RunName'] = lab_run_log['RunDetails_RunName'].str.replace(' ', '')

            # 중복 제거 (BioSampleName, RunDetails_RunName, ServiceId 기준, 마지막 값 유지)
            lab_run_log.drop_duplicates(subset=['BioSampleName', 'RunDetails_RunName', 'ServiceId'], keep='last', inplace=True)
            return lab_run_log
        except FileNotFoundError:
            logger.error(f"세일즈포스 엑셀 파일을 찾을 수 없습니다: {SALESFORCE_EXCEL_PATH}")
            return pd.DataFrame()  # 빈 DataFrame 반환
        except Exception as e:
            logger.error(f"세일즈포스 데이터 추출 중 오류 발생: {e}")
            return pd.DataFrame()  # 빈 DataFrame 반환

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
        """
        ISSUE_SAMPLE_PATH에 정의된 경로에서 이슈 샘플 정보를 찾아 DataFrame으로 반환합니다.
        각 샘플 폴더에서 XML을 파싱하여 BAM 경로와 NumRecords를 추출합니다.
        """
        issue_sample_dir = ISSUE_SAMPLE_PATH
        if not os.path.exists(issue_sample_dir):
            logger.warning(f"이슈 샘플 경로를 찾을 수 없음: {issue_sample_dir}")
            return pd.DataFrame()  # 빈 DataFrame 반환

        issue_samples_folders = [d for d in os.listdir(issue_sample_dir) if os.path.isdir(os.path.join(issue_sample_dir, d))]
        issue_data_list = []

        logger.info(f"이슈 샘플 폴더 {len(issue_samples_folders)}개 검색 시작...")

        for sample_name in issue_samples_folders:
            xml_file_name = f"{sample_name}.consensusreadset.xml"
            xml_file_full_path = os.path.join(issue_sample_dir, sample_name, xml_file_name)

            bam_path_val = None
            num_records_val = -1  # 기본값: 알 수 없음 또는 에러
            total_length_val = -1  # 기본값: 알 수 없음 또는 에러
            started = None  # ExternalResource의 CreatedAt 속성을 저장할 변수

            if os.path.exists(xml_file_full_path):
                root = self.xml_parser.parse_xml(xml_file_full_path)
                if root:
                    # BAM 경로 추출
                    bam_path_df = self.xml_parser.extract_bam_path(root, xml_file_full_path)
                    if bam_path_df is not None and not bam_path_df.empty and 'Path' in bam_path_df.columns:
                        bam_path_val = bam_path_df['Path'].iloc[0]

                    # NumRecords 추출
                    dataset_meta_df = self.xml_parser.extract_dataset_metadata(root)  # format_number_with_commas 적용됨
                    if dataset_meta_df is not None and not dataset_meta_df.empty and 'NumRecords' in dataset_meta_df.columns:
                        num_records_str = dataset_meta_df['NumRecords'].iloc[0]
                        try:
                            num_records_val = int(str(num_records_str).replace(',', ''))  # 쉼표 제거 후 정수 변환
                        except ValueError:
                            logger.warning(f"이슈 샘플 {sample_name}의 NumRecords ('{num_records_str}')를 정수로 변환하는데 실패.")
                            num_records_val = -1  # 변환 실패 시 기본값
                    
                    # TotalLength 추출
                    total_length_df = self.xml_parser.extract_dataset_metadata(root)
                    if total_length_df is not None and not total_length_df.empty and 'TotalLength' in total_length_df.columns:
                        total_length_str = total_length_df['TotalLength'].iloc[0]
                        try:
                            total_length_val = int(str(total_length_str).replace(',', ''))  # 쉼표 제거 후 정수 변환
                        except ValueError:
                            logger.warning(f"이슈 샘플 {sample_name}의 TotalLength ('{total_length_str}')를 정수로 변환하는데 실패.")
                            total_length_val = -1
                    
                    # ExternalResource의 CreatedAt 속성 추출
                    external_resource = root.find('.//pbbase:ExternalResource', self.xml_parser.namespaces)
                    if external_resource is not None and 'CreatedAt' in external_resource.attrib:
                        started = external_resource.attrib.get('CreatedAt')
                        logger.info(f"샘플 {sample_name}의 CreatedAt 추출: {started}")

                else:
                    logger.warning(f"이슈 샘플 XML 파싱 실패: {xml_file_full_path}")
            else:
                logger.warning(f"이슈 샘플 XML 파일을 찾을 수 없음: {xml_file_full_path}")

            issue_data_list.append({
                'BioSampleName': sample_name,
                'XMLFilePath': xml_file_full_path if os.path.exists(xml_file_full_path) else None,
                'BAMPath': bam_path_val,
                'SampleInfo_NumRecords': num_records_val,
                'SampleInfo_TotalLength': total_length_val,
                'RunDetails_WhenStarted': started,
                'RunDetails_WhenCreated': started  
            })

        if not issue_data_list:
            logger.info("처리할 이슈 샘플이 없습니다.")
            return pd.DataFrame()

        issue_df = pd.DataFrame(issue_data_list)
        logger.info(f"이슈 샘플 {len(issue_df)}개 정보 처리 완료.")
        return issue_df


# ==============================================================================
# 리포트 생성 클래스
# ==============================================================================

class ReportGenerator:
    """XML 파싱 결과를 바탕으로 리포트를 생성하는 클래스"""
    
    def __init__(self):
        self.data_processor = DataProcessor()
    
    def process_xml_list(self, xml_files: List[str], output_file: str) -> None:
        """
        XML 파일 목록을 처리하고, 결과를 통합하여 지정된 파일에 저장합니다.
        세일즈포스 데이터와 병합, 특정 조건에 따른 데이터 필터링 및 정렬을 수행합니다.
        """
        result_dfs = []
        processed_files = []
        failed_files = []
        re_demux_df_final = pd.DataFrame()  # 재처리된 이슈 샘플 저장용

        total_files = len(xml_files)
        logger.info(f"총 {total_files}개 XML 파일 처리 시작")

        # 각 XML 파일 처리
        for i, xml_file in enumerate(xml_files):
            try:
                if os.path.exists(xml_file):
                    logger.info(f"파일 처리 중 ({i+1}/{total_files}): {xml_file}")
                    df_processed, _ = self.data_processor.process_xml_file(xml_file)
                    result_dfs.append(df_processed)
                    processed_files.append(xml_file)
                else:
                    logger.warning(f"파일을 찾을 수 없음 ({i+1}/{total_files}): {xml_file}")
                    failed_files.append((xml_file, "파일 없음"))
            except Exception as e:
                logger.error(f"파일 처리 중 오류 발생 ({i+1}/{total_files}): {xml_file}, 오류: {e}")
                failed_files.append((xml_file, str(e)))
                continue  # 다음 파일 처리

        if not result_dfs:
            logger.warning("처리할 수 있는 XML 파일이 없습니다. 결과 파일이 생성되지 않습니다.")
            return

        # 모든 추출된 데이터프레임 결합 및 기본 정리
        combined_df = pd.concat(result_dfs, ignore_index=True)
        combined_df = self.data_processor.clean_dataframe(combined_df)


        # 1차 세일즈포스 데이터와 병합 ( [BioSampleName, RunDetails_RunName] 기준)
        logger.info("1차 세일즈포스 데이터 추출 및 병합 시작...")
        logger.info("세일즈포스 데이터에서 BioSampleName과 RunDetails_RunName 기준으로 병합 시작...")
        # 세일즈포스 데이터 추출
        salesforce_data = self.data_processor.extract_salesforce_data()
        if not salesforce_data.empty:
            if 'RunDetails_RunName' in combined_df.columns:
                combined_df['RunDetails_RunName'] = combined_df['RunDetails_RunName'].dropna().astype(str).str.strip()  # 공백 제거 및 문자열 변환
                combined_df = pd.merge(combined_df, salesforce_data, on=['BioSampleName', 'RunDetails_RunName'], how='left')
            else:
                logger.warning("RunDetails_RunName 컬럼이 combined_df에 없어 세일즈포스 데이터 병합을 건너뜁니다.")
        else:
            logger.warning("세일즈포스 데이터가 비어있어 병합을 건너뜁니다.")

        # 2차 세일즈포스 데이터와 병합 ( [Sample ID, RunDetails_RunName] 기준)
        # salesforce_data copy 후 기존 BioSampleName 컬럼 제거, Sample Id를 BioSampleName으로 변경
        salesforce_data_tmp = salesforce_data.copy()  # 원본 데이터 보호를 위해 복사
        salesforce_data_tmp.drop(columns=['BioSampleName'], inplace=True, errors='ignore') 
        salesforce_data_tmp.rename(columns={'Sample ID': 'BioSampleName'}, inplace=True, errors='ignore')
        
        logger.info("2차 세일즈포스 데이터 추출 및 병합 시작...")
        logger.info("세일즈포스 데이터에서 BioSampleName(엑셀에서 Sample ID열)과 RunDetails_RunName 기준으로 비어있는 ServiceId 채우기 시작...")
        if 'ServiceId' in combined_df.columns and 'BioSampleName' in combined_df.columns:
            # ServiceId가 비어있는 행 필터링
            empty_serviceid_df = combined_df[combined_df['ServiceId'].isna() | (combined_df['ServiceId'] == '')]
            if not empty_serviceid_df.empty:
                logger.info(f"ServiceId가 비어있는 행 {len(empty_serviceid_df)}개 발견. 세일즈포스 데이터로 채우기 시작.")
                
                # BioSampleName, RunDetails_RunName을 기준으로 salesforce_data에서 ServiceId 채우기
                for index, row in empty_serviceid_df.iterrows():
                    bio_sample_name = row['BioSampleName']
                    run_details_run_name = row.get('RunDetails_RunName', '')

                    # salesforce_data에서 해당 BioSampleName과 RunDetails_RunName에 해당하는 행 찾기
                    matching_row = salesforce_data_tmp[
                        (salesforce_data_tmp['BioSampleName'] == bio_sample_name) &
                        (salesforce_data_tmp['RunDetails_RunName'] == run_details_run_name)
                    ]

                    if not matching_row.empty:
                        # ServiceId 채우기
                        combined_df.at[index, 'ServiceId'] = matching_row['ServiceId'].values[0]
                        combined_df.at[index, 'ServiceName'] = matching_row['ServiceName'].values[0]  # ServiceName도 채움
                        logger.info(f"행 {index}의 ServiceId를 {matching_row['ServiceId'].values[0]}로 채움.")
                    else:
                        logger.warning(f"행 {index}의 BioSampleName '{bio_sample_name}'와 RunDetails_RunName '{run_details_run_name}'에 해당하는 ServiceId를 찾을 수 없음.")
            else:
                logger.info("ServiceId가 비어있는 행이 없습니다. 세일즈포스 데이터로 채우기를 건너뜁니다.")
        else:
            logger.warning("combined_df에 ServiceId 또는 BioSampleName 컬럼이 없습니다. 세일즈포스 데이터로 채우기를 건너뜁니다.")
        
        # 3차 세일즈포스 데이터 병합 ( [BioSampleName, RunDetails_WellSampleName_tmp] 기준)
        # RunDetails_WellSampleName_tmp 컬럼 생성 (RunDetails_WellSampleName에서 '-Cell' 앞부분만 추출)
        combined_df['RunDetails_WellSampleName_tmp'] = combined_df['RunDetails_WellSampleName'].str.split('-Cell').str[0]  # RunDetails_WellSampleName에서 '-Cell' 앞부분만 추출
        
        logger.info("3차 세일즈포스 데이터 추출 및 병합 시작...")
        logger.info("세일즈포스 데이터에서 BioSampleName(엑셀에서 Sample ID열)과 RunDetails_WellSampleName_tmp 기준으로 비어있는 ServiceId 채우기 시작...")
        if 'ServiceId' in combined_df.columns and 'BioSampleName' in combined_df.columns:
            # ServiceId가 비어있는 행 필터링
            empty_serviceid_df = combined_df[combined_df['ServiceId'].isna() | (combined_df['ServiceId'] == '')]

            if not empty_serviceid_df.empty:
                logger.info(f"ServiceId가 비어있는 행 {len(empty_serviceid_df)}개 발견. 세일즈포스 데이터로 채우기 시작.")
                
                # BioSampleName, 을 기준으로 salesforce_data에서 ServiceId 채우기
                for index, row in empty_serviceid_df.iterrows():
                    bio_sample_name = row['BioSampleName']
                    run_details_wellsample_name = row.get('RunDetails_WellSampleName_tmp', '')

                    # salesforce_data에서 해당 BioSampleName과 RunDetails_WellSampleName_tmp에 해당하는 행 찾기
                    matching_row = salesforce_data_tmp[
                        (salesforce_data_tmp['BioSampleName'] == bio_sample_name) &
                        (salesforce_data_tmp['RunDetails_WellSampleName_tmp'] == run_details_wellsample_name)
                    ]

                    if not matching_row.empty:
                        # ServiceId 채우기
                        combined_df.at[index, 'ServiceId'] = matching_row['ServiceId'].values[0]
                        combined_df.at[index, 'ServiceName'] = matching_row['ServiceName'].values[0]  # ServiceName도 채움
                        logger.info(f"행 {index}의 ServiceId를 {matching_row['ServiceId'].values[0]}로 채움.")
                    else:
                        logger.warning(f"행 {index}의 BioSampleName '{bio_sample_name}'와 RunDetails_WellSampleName_tmp '{run_details_wellsample_name}'에 해당하는 ServiceId를 찾을 수 없음.")
            else:
                logger.info("ServiceId가 비어있는 행이 없습니다. 세일즈포스 데이터로 채우기를 건너뜁니다.")
        else:
            logger.warning("combined_df에 ServiceId 또는 BioSampleName 컬럼이 없습니다. 세일즈포스 데이터로 채우기를 건너뜁니다.")

    
        combined_df.drop(columns=['Sample ID', 'RunDetails_WellSampleName_tmp'], inplace=True, errors='ignore')

        # BioSampleName이 비어있는 행 제거
        combined_df = combined_df[combined_df['BioSampleName'].notna()]
        combined_df = combined_df[combined_df['BioSampleName'] != "NoName"]  # "NoName"도 제거

        # SampleInfo_NumRecords 필터링 (1000 미만 제거)
        rows_to_remove_for_logging = pd.DataFrame()
        if 'SampleInfo_NumRecords' in combined_df.columns:
            # 콤마 제거 및 정수로 변환 (오류 발생 시 0으로 대체)
            combined_df['SampleInfo_NumRecords_Numeric'] = combined_df['SampleInfo_NumRecords'].astype(str).str.replace(',', '')
            combined_df['SampleInfo_NumRecords_Numeric'] = pd.to_numeric(combined_df['SampleInfo_NumRecords_Numeric'], errors='coerce').fillna(0).astype(int)

            # 1000 미만인 행 식별 (로그 기록용)
            rows_to_remove_for_logging = combined_df[combined_df['SampleInfo_NumRecords_Numeric'] < 1000].drop_duplicates(subset=['BioSampleName', 'BAMPath'])
            if not rows_to_remove_for_logging.empty:
                logger.info(f"SampleInfo_NumRecords가 1000 미만인 {len(rows_to_remove_for_logging)}개의 행이 제거 대상입니다.")
                for _, row in rows_to_remove_for_logging.iterrows():
                    logger.info(f"제거 대상 행: BioSampleName={row.get('BioSampleName', 'N/A')}, SampleInfo_NumRecords={row.get('SampleInfo_NumRecords', 'N/A')}, BAMPath={row.get('BAMPath', 'N/A')}")

            # 1000 이상인 행만 유지
            combined_df = combined_df[combined_df['SampleInfo_NumRecords_Numeric'] >= 1000].copy()  # .copy() 추가
            # SampleInfo_NumRecords를 다시 쉼표 포맷으로 변경
            if not combined_df.empty:
                combined_df['SampleInfo_NumRecords'] = combined_df['SampleInfo_NumRecords_Numeric'].apply(format_number_with_commas)
            combined_df.drop(columns=['SampleInfo_NumRecords_Numeric'], inplace=True, errors='ignore')

            # 이슈 샘플 처리 (제거 대상 중 이슈 샘플 경로에 있는 경우 복원 및 정보 업데이트)
            issue_df = self.data_processor.process_issue_xml()
            if not issue_df.empty and not rows_to_remove_for_logging.empty:
                # 제거 대상이었던 행들 중에서, BioSampleName이 issue_df에 있는 행들을 선택
                re_demux_df = rows_to_remove_for_logging[rows_to_remove_for_logging['BioSampleName'].isin(issue_df['BioSampleName'])].copy()

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
                    re_demux_df_final = re_demux_df.copy()  # 최종 저장용으로 복사
                    combined_df = pd.concat([combined_df, re_demux_df_final], ignore_index=True)

        # ServiceId와 ServiceName 열을 앞쪽으로 이동 (2, 3번째 위치)
        if 'ServiceId' in combined_df.columns and 'ServiceName' in combined_df.columns:
            cols = list(combined_df.columns)
            # 컬럼이 존재하면 이동, 없으면 로그만 남김
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
            combined_df = combined_df[cols]

        # 'RunDetails_WhenCreated' 열 기준으로 시간 순 정렬 (오류 발생 시 무시)
        if 'RunDetails_WhenCreated' in combined_df.columns:
            combined_df['RunDetails_WhenCreated'] = pd.to_datetime(combined_df['RunDetails_WhenCreated'], errors='coerce')
            combined_df = combined_df.sort_values(by='RunDetails_WhenCreated')

        # 최종 중복 제거 (모든 컬럼 기준, 마지막 값 유지)
        combined_df = combined_df.drop_duplicates(keep='last')

        # 첫 번째 column에 행 순서대로 번호 부여
        combined_df.insert(0, 'Number', range(1, len(combined_df) + 1))

        # ServiceName column 옆 column에 ServiceEName 추가
        ename_df = self.data_processor.salesforce_ENAME(combined_df)
        
        # ename_df의 EName column을 combined_df의 ServiceId에 매핑하여 serviceENAME 이라는 새로운 열 추가하고 ServiceName 옆에 위치
        combined_df = pd.merge(combined_df, ename_df, on='ServiceId', how='left')
        combined_df.drop(columns=['EName'], inplace=True, errors='ignore')

        # ServiceEName 열을 ServiceName 옆으로 이동
        if 'ServiceEName' in combined_df.columns:
            cols = list(combined_df.columns)
            cols.remove('ServiceEName')
            # ServiceEName을 ServiceName 옆으로 이동
            service_name_index = cols.index('ServiceName') + 1
            cols.insert(service_name_index, 'ServiceEName')
            combined_df = combined_df[cols]
        else:
            logger.warning("ServiceEName 컬럼이 없어 위치를 조정할 수 없습니다.")
        
        # 결과 저장
        try:
            combined_df.to_csv(output_file, index=False, sep="\t")
            logger.info(f"통합 데이터가 {output_file}에 저장되었습니다. 총 {len(combined_df)}개의 레코드가 있습니다.")
        except Exception as e:
            logger.error(f"결과 파일 저장 중 오류 발생 ({output_file}): {e}")

        # 처리된 파일 및 실패한 파일 목록 저장
        base_name = os.path.splitext(output_file)[0]
        processed_files_path = f"{base_name}{PROCESSED_FILES_SUFFIX}"
        with open(processed_files_path, 'w') as f:
            for xml_file_path in processed_files:
                f.write(f"{xml_file_path}\n")
        logger.info(f"처리된 파일 목록이 {processed_files_path}에 저장되었습니다.")

        if failed_files:
            failed_files_path = f"{base_name}{FAILED_FILES_SUFFIX}"
            with open(failed_files_path, 'w') as f:
                for xml_file_path, reason in failed_files:
                    f.write(f"{xml_file_path}\t{reason}\n")
            logger.warning(f"{len(failed_files)}개 파일 처리 실패, 목록이 {failed_files_path}에 저장되었습니다.")

        # 요약 로그
        logger.info(f"SUMMARY :: TOTAL PROCESSED XMLs :: {len(processed_files)}")
        logger.info(f"SUMMARY :: FAILED XMLs          :: {len(failed_files)}")
        logger.info(f"SUMMARY :: ROWS REMOVED (<1000 reads initially) :: {len(rows_to_remove_for_logging)}")
        logger.info(f"SUMMARY :: ROWS RE-DEMUXED (issue samples) :: {len(re_demux_df_final)}")
        logger.info(f"SUMMARY :: FINAL RECORDS IN OUTPUT :: {len(combined_df)}")

    def process_single_xml(self, xml_file: str, output_file: str) -> Optional[str]:
        """단일 XML 파일을 처리하고, BioSample 이름을 포함한 파일명으로 결과를 저장합니다."""
        if not os.path.exists(xml_file):
            logger.error(f"입력 XML 파일을 찾을 수 없음: {xml_file}")
            return None

        logger.info(f"단일 XML 파일 처리 중: {xml_file}")
        result_df, biosample_name = self.data_processor.process_xml_file(xml_file)

        if "Error" in biosample_name or result_df.empty:
            logger.error(f"단일 XML 파일 처리 실패: {xml_file}")
            return None

        result_df = self.data_processor.clean_dataframe(result_df)

        # 출력 파일명 결정 (BioSampleName 기반)
        base_output_name = os.path.splitext(output_file)[0]
        file_extension = os.path.splitext(output_file)[1] or ".tsv"  # 확장자 없으면 .tsv 사용
        
        # BioSampleName이 유효한 경우 파일명에 추가
        if biosample_name not in ["NoName", "Error", None, ""]:
            output_path = f"{base_output_name}_{biosample_name}_report{file_extension}"
        else:
            output_path = output_file  # 기본 출력 파일명 사용

        # 결과 저장
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
    """
    메인 실행 함수. 명령줄 옵션을 파싱하고, 옵션에 따라 XML 처리 방식을 결정합니다.
    - XMLLIST 파일이 제공되면 해당 목록의 XML들을 처리합니다.
    - REVIOXML 파일이 제공되면 해당 단일 XML 파일을 처리합니다.
    - 아무 옵션도 없으면 기본 경로에서 모든 XML 파일을 찾아 처리합니다.
    """
    options = get_options()
    exit_code = 0  # 기본 성공 코드
    report_generator = ReportGenerator()

    try:
        if options.xmllist and options.xmllist != "None":
            # XML 목록 파일 처리
            if not os.path.exists(options.xmllist):
                logger.error(f"XML 목록 파일을 찾을 수 없음: {options.xmllist}")
                return 1  # 오류 코드 반환
            try:
                with open(options.xmllist, 'r') as f:
                    xml_files_from_list = [line.strip() for line in f if line.strip() and line.strip().endswith('.xml')]
                if not xml_files_from_list:
                    logger.error(f"XML 목록 파일이 비어 있거나 유효한 XML 파일이 없음: {options.xmllist}")
                    return 1
                logger.info(f"XML 목록 파일에서 {len(xml_files_from_list)}개의 파일을 찾았습니다.")
                report_generator.process_xml_list(xml_files_from_list, options.output)
            except Exception as e:
                logger.error(f"XML 목록 파일 처리 중 오류 발생: {e}")
                return 1

        elif options.revioxml and options.revioxml != "None":
            # 단일 XML 파일 처리
            if not os.path.exists(options.revioxml):
                logger.error(f"XML 파일을 찾을 수 없음: {options.revioxml}")
                return 1
            report_generator.process_single_xml(options.revioxml, options.output)

        else:
            # 기본 경로에서 모든 XML 파일 처리
            logger.info(f"지정된 XML 파일 또는 목록이 없어, 기본 경로({ALL_XML_GLOB_PATTERN})에서 모든 XML 파일을 처리합니다.")
            all_xmls_found = glob.glob(ALL_XML_GLOB_PATTERN)
            # 제외 패턴에 해당하는 파일 필터링
            xml_files_to_process = [f for f in all_xmls_found if not EXCLUDE_XML_PATTERN_REGEX.search(f)]
            if not xml_files_to_process:
                logger.warning("처리할 XML 파일이 기본 경로에 없습니다.")
            else:
                logger.info(f"기본 경로에서 {len(xml_files_to_process)}개의 XML 파일을 처리합니다.")
                report_generator.process_xml_list(xml_files_to_process, options.output)

    except Exception as e:
        logger.error(f"메인 실행 중 예상치 못한 오류가 발생했습니다: {e}", exc_info=True)  # 상세 스택 트레이스 로깅
        exit_code = 1  # 오류 코드 설정
    
    return exit_code


if __name__ == "__main__":
    # 스크립트 실행 시작 로깅
    logger.info("==============================================================")
    logger.info("Revio Sequencing Report Generation 스크립트 시작")
    logger.info(f"실행 인수: {' '.join(sys.argv)}")
    logger.info("==============================================================")
    
    script_exit_code = main()
    
    # 스크립트 실행 종료 로깅
    logger.info("==============================================================")
    if script_exit_code == 0:
        logger.info("Revio Sequencing Report Generation 스크립트 성공적으로 완료")
    else:
        logger.error(f"Revio Sequencing Report Generation 스크립트 오류로 종료 (Exit Code: {script_exit_code})")
    logger.info("==============================================================")
    
    sys.exit(script_exit_code)