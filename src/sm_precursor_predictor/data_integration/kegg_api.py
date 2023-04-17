import io

from Bio.KEGG import REST as kegg_api
import pandas as pd


class KeggApi:

    @staticmethod
    def to_df(result):
        return pd.read_table(io.StringIO(result), header=None)

    @staticmethod
    def get_list(database):
        list_of_pathways = kegg_api.kegg_list(database).read()
        return KeggApi.to_df(list_of_pathways)

    @staticmethod
    def get_info(database):
        info = kegg_api.kegg_info(database).read()
        return KeggApi.to_df(info)

    @staticmethod
    def get(query, **kwargs):
        result = kegg_api.kegg_get(query, **kwargs).read()
        return result

    @staticmethod
    def get_links(database1, database2):
        result = kegg_api.kegg_link(database1, database2).read()
        return result

    @staticmethod
    def find(database, entry):
        result = kegg_api.kegg_find(database, entry).read()
        return result
