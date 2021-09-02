rule get_vep_cache:
    output:
        directory(f"{VEP_CACHE}")
    params:
        species=get_params('annotate','species'),
        build=get_params('annotate','build'),
        release=get_params('annotate','release')
    log:
        f"{LOGDIR}/get_vep_cache/cache.log"
    cache: True  # save space and time with between workflow caching (see docs)
    wrapper:
        "0.77.0/bio/vep/cache"