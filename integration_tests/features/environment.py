import sys, os


def before_all(context):
    if context.config.userdata['datapath'].startswith(os.sep):
        context.data_path = context.config.userdata['datapath']
    else:
        context.data_path = os.path.join(os.getcwd(), context.config.userdata['datapath'])
    context.path_to_scripts = os.path.expandvars('$PATH_TO_SCRIPTS')


def after_all(context):
    if context.failed:
        sys.exit(1)
