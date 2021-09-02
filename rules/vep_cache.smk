rule get_vep_cache:
    output:
        directory(f"{VEP_CACHE}")
    params:
        species=get_params('annotate','species'),
        build=get_params('annotate','build'),
        release=get_params('annotate','release')
    threads:
        get_resource('default', 'threads')
    resources:
        mem=get_resource('default', 'mem'),
        walltime=get_resource('default', 'walltime')
    log:
        f"{LOGDIR}/get_vep_cache/cache.log"
    cache: True  # save space and time with between workflow caching (see docs)
    wrapper:
        "0.77.0/bio/vep/cache"