import sys
import logging
import ConfigParser
import logging

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

def loadconfig(args):
    cfg = args['config']
    config = ConfigParser.ConfigParser()
    config.readfp(open(cfg))
    log.info("Loading config %s", cfg)
    for k,v in args.items():
        if config.has_option('default',k):
            log.info("overriding %s: %s"%(k,v))
            config.set('default',k,v)
    #for k,v in config.items():

    arncfg = ConfigParser.ConfigParser()
    arncfg.readfp(open(config.get('default','arnconf')))
    log.info("Loading arn config %s", config.get('default','arnconf'))
    for k,v in args.items():
        if arncfg.has_option('default',k):
            log.info("overriding %s: %s"%(k,v))
            arncfg.set('default',k,v)

    return config, arncfg

def parsecmd():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='The main configuration file \
                            must be provided. Any parameter can be overwritten \
                            using --pname pvalue')
    parser.add_argument("config", metavar="FILE",
                        help="main configuration file")

    try:
        cfg =  ConfigParser.ConfigParser()
        cfg.readfp(open(sys.argv[-1]))
        for k,v in cfg.items('default'):
            parser.add_argument("--"+k, dest=k,
                            help = "Overwrites the configuration file value")
    except IOError:
        sys.argv.append('-h')

    parser.add_argument("-i", "--initdm", dest="initdm", metavar="NUM",
                        help="number of dm events to intialize population \
                        (controls genome size if initialization is random)")
    parser.add_argument("-o", "--overlapgenes", dest="overlapgenes",
                        metavar="BOOL",
                        help="if True genes and binding sites may overlap (but not promotors)")
    args = vars(parser.parse_args())
    return dict(filter(lambda kv: kv[1] != None,
                  args.items()))
    #options = filter(lambda k,v: v != None,
     #                options.items())
    #return options, args