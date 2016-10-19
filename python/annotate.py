"""
Functions for the annotation workflows
"""


def get_hit_table_name_from_wildcards_db(wildcards, config, name_root='genes'):
    """
    Return the hit table name based on the db name using the db config info
    """
    db=wildcards.db
    db_type=config['dbs'][db].get('type','hmmer')
    if db_type=='hmmer':
        if 'frags' in config['dbs'][db]:
            template = "{name_root}.vs.{db}.tbl.dbatch"
        else:
            template = "{name_root}.vs.{db}.tbl"
    elif db_type=='lastdb':
        template = "{name_root}.vs.{db}.lastp"
    else:
        # Don't know what to do:
        raise Exception("Unknown database type for {}: {}".format(db,db_type))

    return template.format(name_root=name_root, **wildcards)


def get_db_assignment_params(wildcards, config):
    """
    return the params needed to turn hits from the given db (wildcards.db) into gene family assignments 
    using the assign_paths.py script
    """
    assign_type = config['dbs'][wildcards.db].get('assign_type','hitid').lower()
    if assign_type=='kegg':
        return '-p hitid -M kegg -m %s.kos' % (config['dbs'].get('KEGG',{'path':''})['path'])
    if assign_type=='pfam':
        return '-p pfam'
    return '-p hitid'


