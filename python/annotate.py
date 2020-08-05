"""
Functions for the annotation workflows
"""
import os

def get_db_types(config):
    """ loop over DBs and:
         * identify all gene family dbs
         * identify all taxonomic dbs
    """
    gene_family_dbs = []
    taxdbs = []
    for db in config['dbs']:
        db_type = config['dbs'][db].get('type', 'gene')
        if db_type[0:3] == 'tax':
            taxdbs.append(db)
        elif db_type == 'gene':
            gene_family_dbs.append(db)
    return (gene_family_dbs, taxdbs)

def get_last_alg(dbformat, extension):
    """
    right now looks for last db type (lastp or lastn) and extension (faa or
    not) and returns lastp, lastx, or lastn.

    Support for other dbs can be added on request.
    """
    if dbformat == 'lastp':
        if extension == 'faa':
            search_alg = 'lastp'
        else:
            search_alg = 'lastx'
    elif dbformat == 'lastn':
        if extension == 'faa':
            raise Exception("I'm sorry, I don't know how to search for faa "
                            "sequences in a lastp database!")
        else:
            search_alg = 'lastn'
    elif dbformat == 'bwadb':
        search_alg = 'bwa.sam'
    elif dbformat == 'dmnd':
        if extension == 'faa':
            search_alg = 'dmndp'
        else:
            search_alg = 'dmndx'
    else:
        raise Exception(("I'm sorry, but the database format '{}' is not yet "
                         "supported").format(dbformat))
    return search_alg


def get_db_dot_fmt_strings(db_list, config, query_extension='fasta'):
    """
    Return a list of strings that are "{db}.{format}". Where db is the name of the database and format is the extension generated by the search (eg lastx, or tbl). There is a special case for fragmented HMM dbs where we need to add ".dbatch" to the format.
    """
    strings=[]
    for d in db_list:
        db_data = config['dbs'][d]
        format = db_data.get('format','tbl')
        if format.startswith('last'):
            format = get_last_alg(format, query_extension)
        if 'frags' in db_data and int(db_data['frags'])>1:
            format = format + ".dbatch"
        strings.append("{}.{}".format(d,format))
    return strings


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
        return '-p hitid -C first -M kegg -m %s.kos' % (config['dbs'].get('KEGG',{'path':''})['path'])
    if assign_type=='pfam':
        return '-p pfam -C all'
    return '-C first -p hitid'


def get_db_frag(config, db, N):
    full_hmm = config['dbs'][db]['path']
    n_frags = config['dbs'][db]['frags']
    template = get_db_frag_template(full_hmm, n_frags)
    return template.format(N=int(N))

def get_db_frag_template(full_hmm, n_frags):
    n_frag_digits = 3  # unless we can get the hard coded 000 out of the rules
    #n_frag_digits = len(str(n_frags))
    hmm_dir = os.path.dirname(full_hmm)
    frag_dir = os.path.join(hmm_dir, 'frag_{}'.format(n_frags))
    hmm_base, hmm_ext = os.path.splitext(os.path.basename(full_hmm))
    template = "{}{}{}.{{N:0{}d}}{}".format(frag_dir, os.path.sep,
                                            hmm_base, n_frag_digits,
                                            hmm_ext)
    return template

def get_db_frags(full_hmm, n_frags):
    template = get_db_frag_template(full_hmm, n_frags)
    return [template.format(N) for N in range(1, n_frags+1)]
