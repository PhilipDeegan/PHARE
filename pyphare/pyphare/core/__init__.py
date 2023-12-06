
def dbg_print_(*args):
    print(*args)

dbg_print = dbg_print_ if __debug__ else lambda *args: ...
