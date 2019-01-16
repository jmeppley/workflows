"""
functions to aid in dowloading remote files
"""
import re, os
from python.common import apply_defaults
from snakemake import logger
from snakemake.remote.SFTP import RemoteProvider as sftp_rp
from snakemake.io import ancient

url_rexp = re.compile(r'^([A-Za-z]+)://(.+)$')
mnt_rexp = re.compile(r'^(/mnt/([a-z]+ine)/.+)$')

# cache redentials for known hosts
providers = {}
def get_provider(protocol, host, config):
    """
    Creates a unique RemoteProvider object for each (protocol, host) pair.

    Connection options taken from:
    config['remote'][protocol][host]

    With defaults from:
    config['remote'][protocol]['defaults']

    protocol must be one of SFTP, SCP, HTTP, HTTPS, FTP, but case can be
    different. Case is used when keying options, though, so be consistent in
    your configuration
    """
   
    provider_key = (protocol, host)
    protocol = protocol.upper()
    if provider_key not in providers:
        remote_defaults = config.get('remote', {}) \
                                .get(protocol, {}) \
                                .get('defaults', {})
        remote_options = config.get('remote', {}) \
                                .get(protocol, {}) \
                                .get(host, {})
        apply_defaults(remote_options, remote_defaults)
        logger.debug("Creating RemoteProvider for {protocol} "
                     "using:\n{remote_options}".format(**vars()))
        if protocol in ['SFTP', 'SCP']:
            providers[provider_key] = sftp_rp(**remote_options)
        elif protocol in ['HTTP', 'HTTPS']:
            from snakemake.remote.HTTP import RemoteProvider as http_rp
            providers[provider_key] = http_rp(**remote_options)
        elif protocol == 'FTP':
            from snakemake.remote.FTP import RemoteProvider as ftp_rp
            providers[provider_key] = ftp_rp(**remote_options)
        else:
            raise Exception("Remote protocol {} not yet supported"
                            .format(provider_key[0]))
    return providers[provider_key]

def path_up_to_wildcard(full_path):
    """ If a given path has a wildcard placeholder ( eg {sample} ), 
    return the last directory before that point """
    path_fragment = full_path.split('{')[0]
    if path_fragment.endswith(os.path.pathsep):
        return path_fragment[:-1]
    return os.path.dirname(path_fragment)

def infer_provider(source, config, glob=False):
    """
    Checks the file path to see if we need a remote provider
    """
    try:
        # is it an explicit url
        # (EG: SFTP://lys.soest.hawaii.edu/mnt/lysine/...)
        m = url_rexp.search(source)
        if m is not None:
            logger.debug("EXPLICIT URL")
            # replace source file with a remote object
            protocol, source_path = m.groups()
            host = re.search(r'^([^/]+)/', source_path).group(1)
            return get_provider(protocol, host, config), source_path

        # special case: delong lab remote mounts
        m = mnt_rexp.search(source)
        if m is not None and not os.path.exists(path_up_to_wildcard(source) if glob else source):
            logger.debug("INFERRED URL")
            protocol = 'SFTP'
            config.setdefault('remote', {}) \
                  .setdefault(protocol, {}) \
                  .setdefault('defaults', {'username': 'readonly',
                                          })
            source_path, server_name = m.groups()
            host = server_name + config.get('remote', {}) \
                                       .get('domain', '.soest.hawaii.edu')
            source_path = host + source_path
            #return get_provider(protocol, host), "".join([host, source])
            return get_provider(protocol, host, config), source_path

    except Exception as e:
        print("Error in remote check: " + repr(e))
        raise e

    return None, source

def remote_wrapper(source, config, glob=False):
    """
    if file is a remote url 
         ( or a missing netowrk mount )
    return remote provider object for downloading
    or return wildcard lists if glob=True (using glob_wildcards)
    """
    provider, source = infer_provider(source, config, glob=glob)
    logger.debug("provider: {}\nsource: {}".format(provider, source))
    if glob:
        if provider is None:
            return glob_wildcards(source)
        else:
            return provider.glob_wildcards(source)
    else:
        if provider is None:
            return source
        else:
            # return provider.remote(source)
            # Hack to bypass bug in SFTP
            host, path = source.split("/", 1)
            full_path = "/" + path
            if isinstance(provider, sftp_rp):
                # download with rsync
                local_path = ancient("rsync/" + source)
                config.setdefault('download_map', {})[source] = \
                        {'host': host, 'remote_path': full_path,
                         'user': provider.kwargs['username']}
            else:
                # use remote()
                local_path = 'remote/' + source
                config.setdefault('download_map', {})[source] = \
                        {'remote': provider.remote(source)}
            return local_path

