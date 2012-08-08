from yt.mods import *

def main():
    import enzo
    import message_passing
    print dir(message_passing)
    pf = EnzoStaticOutputInMemory()
    total_cells = sum(g.ActiveDimensions.prod() for g in pf.h.grids)
    message_passing.test_event(total_cells)
