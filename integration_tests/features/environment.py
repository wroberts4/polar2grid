import sys, os


def before_all(context):
    if context.config.userdata['datapath'].startswith(os.sep):
        context.data_path = context.config.userdata['datapath']
    else:
        context.data_path = os.path.join(os.getcwd(), context.config.userdata['datapath'])
    p2g_home = os.path.expandvars('$POLAR2GRID_HOME')
    context.p2g_path = os.path.join(p2g_home, 'bin')


def after_all(context):
    if context.failed:
        sys.exit(1)
