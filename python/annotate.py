"""
Functions for the annotation workflows
"""


def get_db_types(config):
    gene_family_dbs = []
    for db in config['dbs']:
        db_type = config['dbs'][db].get('type','gene')
        if db_type == 'tax':
            taxdb=db
        elif db_type == 'gene':
            gene_family_dbs.append(db)
    # return, but make sure we found a tax db
    try:
        return (gene_family_dbs, taxdb)
    except NameError:
        raise Exception("No taxonomic db found. One of your configured annotation"
                        " databases must have 'istaxdb' set to True!")


def get_hit_table_name_from_wildcards_db(wildcards, config):
    """
    Return the hit table name based on the db name using the db config info
    """
    db=wildcards.db
    db_format=config['dbs'][db].get('format','hmmer')
    if db_format=='hmmer':
        if 'frags' in config['dbs'][db]:
            template = "{name_root}.vs.{db}.tbl.dbatch"
        else:
            template = "{name_root}.vs.{db}.tbl"
    elif db_format=='lastdb':
        template = "{name_root}.vs.{db}.lastp"
    else:
        # Don't know what to do:
        raise Exception("Unknown database format for {}: {}".format(db,db_format))

    name_root = config['annotation_hit_table_map']\
                        .get(wildcards.annotation_prefix,
                             wildcards.annotation_prefix)

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


