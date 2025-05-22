from simple_salesforce import Salesforce
import argparse, sys
from tabulate import tabulate
import json
import yaml
import pandas as pd
import tabulate


# def full_data():
#     merged_df = pd.read_csv("~/bin/sequel_full_db.tsv", sep="\t")
#     return merged_df


def _get_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        '-s',
        dest='K_NAME',
        help="Korean Name or Client info")
    parser.add_argument(
        '-n',
        dest='E_NAME',
        help="English Name or Client info")
    parser.add_argument(
        '-e',
        dest='ExID',
        help="Experiment ID")
    parser.add_argument(
        '-T',
        dest='Table',
        action='store_true',
        help="Print table format in terminal")
    parser.add_argument(
        '-y',
        dest='Yaml',
        action='store_true',
        help="Print yaml format")
    parser.add_argument(
        '-j',
        dest='Json',
        action='store_true',
        help="Print json format(default)")
    parser.add_argument(
        '-c',
        dest='csv',
        action='store_true',
        help="Save as CSV file")
    parser.add_argument(
        '-t',
        dest='tsv',
        action='store_true',
        help="Save as TSV file")


    if len(sys.argv)==0:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    Service = args.K_NAME
    Name = args.E_NAME
    ExID = args.ExID
    Table = args.Table
    Yaml = args.Yaml
    Json = args.Json
    CSV = args.csv
    TSV = args.tsv
    return(Service, Name, ExID, Table, Yaml, Json, CSV, TSV)

Service, Name, ExID, Table, Yaml, Json, CSV, TSV  = _get_args()

sf = Salesforce(
    username='bioinfo@dnalink.com',
    password='dnabio$191!',
    security_token='egCwveNX3jPiUslfYBx4qfH8T'
)

# 사용자 지정 객체의 이름 설정 (예: 'Experiment__c')
custom_object_name = 'Experiment__c'

# 필드 이름 설정 (레퍼런스 필드 이름, 예: 'Contact__c')
reference_field_name = 'Contact__c' #고객정보
reference_field_name2 = 'SalesUser__c' # 영업담당자 정보
reference_field_name3 = 'EName__c'
Direction_field_name = 'AnMemo__c'
service_field_name = 'ServiceName__c'

# SOQL 쿼리 설정 (검색어를 이용한 데이터 검색, Contact__c 필드 포함)

if Service:
    ServiceID_keyword = Service
    soql_query = f"SELECT {service_field_name}, {reference_field_name3}, {reference_field_name}, {reference_field_name2}, {Direction_field_name} FROM {custom_object_name} WHERE ServiceName__c LIKE '%{ServiceID_keyword}%'"# LIMIT 10"
elif Name:
    client_name_keyword = Name
    soql_query = f"SELECT {service_field_name}, {reference_field_name3}, {reference_field_name}, {reference_field_name2}, {Direction_field_name} FROM {custom_object_name} WHERE EName__c LIKE '%{client_name_keyword}%'"# LIMIT 10"
elif ExID:
    Expriment_keyword = ExID
    soql_query = f"SELECT {service_field_name}, {reference_field_name3}, {reference_field_name}, {reference_field_name2}, {Direction_field_name} FROM {custom_object_name} WHERE Name LIKE '%{Expriment_keyword}%'"# LIMIT 10"
else:
    sys.exit("No options provided. Use -s for Korean ServiceID/Name or -n for English ServiceID/Name or -e for experiment ID.")

try:
    data_sum = []
    # Salesforce에 SOQL 쿼리를 보내서 결과를 가져오기
    query_result = sf.query_all(soql_query)
    records = query_result['records']

    Index=0
	# match service information and client, Sales info
    for record in query_result['records']:
        reference_value = record[reference_field_name]
        reference_value2 = record[reference_field_name2]

        contact_query = f"SELECT Id, Name, Email, Phone FROM Contact WHERE Id = '{reference_value}'"
        contact_result = sf.query_all(contact_query)

        contact_query2 = f"SELECT Id, Name, Email, Phone FROM Contact WHERE Id = '{reference_value2}'"
        contact_result2 = sf.query_all(contact_query2)

		# query result merge
        data_1=query_result['records'][Index]
        data_1['Contact__c'] = contact_result['records'][0]
        data_1['SalesUser__c'] = contact_result2['records'][0]
        data_sum.append(data_1)
        Index+=1

	#Print result data
    if Table == True:
        #Table_output = tabulate(data_sum)
        #Table_output = '\t'.join(data_sum) 리스트가 아니라서 불가. collection.OrderedDict 타입이다.
        print(data_sum)

    if Yaml == True:
        yaml_output = yaml.dump(data_sum, default_flow_style=False, allow_unicode=True)
        print(yaml_output)
        #with open("query_result.yaml", "w", encoding='utf-8') as yaml_file:
        #    yaml_file.write(yaml_output)
    else :
        json_output = json.dumps(data_sum, indent=4, ensure_ascii=False)
        print(json_output)
        #with open('query_result.json', 'w', encoding='utf-8') as json_file:
        #    json.dump(data_sum, json_file, indent=4, ensure_ascii=False)

except Exception as e:
    print(f"Error: {e}")
