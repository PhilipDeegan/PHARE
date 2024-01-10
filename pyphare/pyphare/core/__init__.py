#
#
#
#
#


def _print(obj):
    for attr in dir(obj):
        if hasattr(obj, attr):
            print("obj.%s = %s" % (attr, getattr(obj, attr)))
